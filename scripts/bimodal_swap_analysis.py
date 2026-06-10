#!/usr/bin/env python3
# bimodal_swap_analysis.py — Estimate MPKI impact of swapping the bimodal
# fallback inside a TAGE monitor run with a different bimodal predictor.
#
# Inputs:
#   --monitor DIR      Monitor run dir (e.g. full_out/monitor_mb2_all/) with
#                      per-trace <name>_summary.txt files.
#   --swap NAME=DIR    Path to a non-monitor full-eval dir (.out files from
#                      run_all / quick_eval) for the candidate bimodal.
#                      Repeatable, e.g. --swap bimA=full_out/bimA
#                                       --swap bimB=full_out/bimB
#
# For each trace it extracts:
#   - From monitor summary: instructions, total mispredicts, MPKI baseline,
#     bimodal provisioned branches B, bimodal correct C  → bimodal_misp = B-C.
#   - From swap .out (CSV): condbr, mispredicts → swap accuracy on this trace.
#
# Approximation: treats the swap predictor's overall conditional-branch
# accuracy as a proxy for the accuracy it would achieve on the FB-routed
# subset.  Computes:
#   new_misp = total_misp - bimodal_misp + B * (1 - swap_acc)
#   new_MPKI = new_misp / instructions * 1000
#   delta_MPKI = new_MPKI - baseline_MPKI

import argparse
import csv
import re
import sys
from pathlib import Path


CBP_CSV_FIELDS = [
    "name", "ninstr", "nbranch", "ncondbr", "npred", "extra_cycles",
    "p1p2_disagree", "p1p2_at_block_end", "mispredictions",
    "p1_latency", "p2_latency", "epi_fJ",
]


def parse_monitor_summary(path: Path):
    """Return dict with instr, branches, mispred, mpki, bim_prov, bim_corr."""
    text = path.read_text()
    m = re.search(
        r"Instructions:\s+(\d+)\s+Branches:\s+(\d+)\s+Mispredictions:\s+(\d+)"
        r"[^M]*MPKI:\s+([\d.]+)",
        text,
    )
    if not m:
        return None
    instr, branches, mispred, mpki = int(m[1]), int(m[2]), int(m[3]), float(m[4])

    # FB per-branch accuracy comes from "Provider Source Breakdown":
    #   "No tag match (FB only):    1329596 (28.2%)  acc=63.1%"
    #   "Sec-tag rejected (FB):     1295999 (27.5%)  acc=56.9%"
    fb_rows = re.findall(
        r"(?:No tag match \(FB only\)|Sec-tag rejected \(FB\)):\s+(\d+)\s+\([\d.]+%\)\s+acc=([\d.]+)%",
        text,
    )
    if len(fb_rows) < 2:
        return None
    bim_prov = 0
    bim_corr = 0
    for cnt_s, acc_s in fb_rows:
        cnt = int(cnt_s)
        acc = float(acc_s) / 100.0
        bim_prov += cnt
        bim_corr += round(cnt * acc)

    # Real FB-caused block mispredicts from "Block-mispredict Cause Attribution":
    #   "No tag match (FB only):    1234 (8.4%)"
    #   "Sec-tag rejected (FB):     5678 (38.5%)"
    # FB-caused total is what's actually attributable to the bimodal — what we
    # need for swap-impact calculation (not the phantom-inflated per-branch
    # routing counts).
    fb_caused = None
    cm = re.search(
        r"Block-mispredict Cause Attribution:.*?FB-caused total:\s+(\d+)",
        text, re.DOTALL,
    )
    if cm:
        fb_caused = int(cm[1])

    return dict(
        instr=instr, branches=branches, mispred=mispred, mpki=mpki,
        bim_prov=bim_prov, bim_corr=bim_corr,
        fb_caused=fb_caused,
    )


def parse_cbp_csv(path: Path):
    """Read a single-line .out csv from run/quick_eval."""
    with path.open() as f:
        row = next(csv.reader(f))
    if len(row) < len(CBP_CSV_FIELDS):
        return None
    rec = dict(zip(CBP_CSV_FIELDS, row))
    return dict(
        instr=int(rec["ninstr"]),
        condbr=int(rec["ncondbr"]),
        mispred=int(rec["mispredictions"]),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--monitor", required=True, type=Path,
                    help="Monitor run dir with <trace>_summary.txt files")
    ap.add_argument("--baseline", type=Path, default=None,
                    help="Non-monitor full-eval dir for the same predictor "
                         "(.out files). Used as the reference MPKI for delta "
                         "computation; monitor MPKI shown alongside.")
    ap.add_argument("--swap", action="append", default=[],
                    metavar="NAME=DIR",
                    help="Candidate bimodal: name=path/to/full_out_dir (repeatable)")
    ap.add_argument("--sort", choices=("name", "delta", "mpki"), default="name")
    args = ap.parse_args()

    swaps = {}
    for spec in args.swap:
        if "=" not in spec:
            sys.exit(f"--swap requires NAME=DIR (got {spec!r})")
        name, d = spec.split("=", 1)
        swaps[name] = Path(d)
    if not swaps:
        sys.exit("need at least one --swap NAME=DIR")

    rows = []
    for summary in sorted(args.monitor.glob("*_summary.txt")):
        trace = summary.name[: -len("_summary.txt")]
        mon = parse_monitor_summary(summary)
        if mon is None:
            continue

        row = {"trace": trace, **mon,
               "bim_prov_pct": 100.0 * mon["bim_prov"] / mon["branches"]
                                if mon["branches"] else 0.0,
               "bim_acc": (mon["bim_corr"] / mon["bim_prov"])
                          if mon["bim_prov"] else 0.0}

        # Non-monitor baseline MPKI (option A: reference for deltas)
        row["base_mpki"] = None
        if args.baseline is not None:
            bp = args.baseline / f"{trace}.out"
            if bp.exists():
                b = parse_cbp_csv(bp)
                if b is not None and b["instr"]:
                    row["base_mpki"] = b["mispred"] / b["instr"] * 1000.0
        ref_mpki = row["base_mpki"] if row["base_mpki"] is not None else mon["mpki"]
        row["ref_mpki"] = ref_mpki

        for name, d in swaps.items():
            out = d / f"{trace}.out"
            if not out.exists():
                row[f"{name}_mpki"] = None
                row[f"{name}_delta"] = None
                continue
            sw = parse_cbp_csv(out)
            if sw is None or sw["condbr"] == 0:
                row[f"{name}_mpki"] = None
                row[f"{name}_delta"] = None
                continue
            swap_acc = 1.0 - sw["mispred"] / sw["condbr"]
            # Use FB-caused block mispredicts (real attribution) when the
            # monitor produced them; otherwise fall back to per-branch math
            # (which over-estimates the swap delta for low-MPKI traces).
            bim_acc = (mon["bim_corr"] / mon["bim_prov"]
                       if mon["bim_prov"] else 0.0)
            if mon.get("fb_caused") is not None and (1.0 - bim_acc) > 1e-9:
                # Scale FB-caused mispredicts by the accuracy ratio.
                fb_caused = mon["fb_caused"]
                ratio = (1.0 - swap_acc) / (1.0 - bim_acc)
                new_fb_caused = fb_caused * ratio
                new_total_misp = mon["mispred"] - fb_caused + new_fb_caused
            else:
                # Legacy / fallback path: per-branch substitution.
                bim_misp = mon["bim_prov"] - mon["bim_corr"]
                new_bim_misp = mon["bim_prov"] * (1.0 - swap_acc)
                new_total_misp = mon["mispred"] - bim_misp + new_bim_misp
            new_mpki = new_total_misp / mon["instr"] * 1000.0
            row[f"{name}_acc"] = swap_acc
            row[f"{name}_mpki"] = new_mpki
            row[f"{name}_delta"] = new_mpki - ref_mpki
        rows.append(row)

    if args.sort == "name":
        rows.sort(key=lambda r: r["trace"])
    elif args.sort == "mpki":
        rows.sort(key=lambda r: -r["mpki"])
    elif args.sort == "delta":
        first = next(iter(swaps))
        rows.sort(key=lambda r: (r.get(f"{first}_delta") or 0.0))

    # Column layout: trace name left-aligned wide; numerics right-aligned
    trace_w = max([len(r["trace"]) for r in rows] + [len("trace")])

    fixed_cols = [
        ("trace",     trace_w, "<"),
        ("MPKI_mon",       9, ">"),
        ("MPKI_base",      9, ">"),
        ("bim_prov%",      9, ">"),
        ("bim_acc%",       8, ">"),
    ]
    swap_cols = []
    for name in swaps:
        swap_cols += [
            (f"{name}_acc%",   max(8, len(name) + 5), ">"),
            (f"{name}_MPKI",   max(9, len(name) + 5), ">"),
            (f"{name}_dMPKI",  max(9, len(name) + 6), ">"),
        ]
    cols = fixed_cols + swap_cols

    def fmt_row(values):
        return "  ".join(f"{v:{align}{w}}" for v, (_, w, align) in zip(values, cols))

    header_vals = [name for name, _, _ in cols]
    print(fmt_row(header_vals))
    print("  ".join("-" * w for _, w, _ in cols))

    tot_instr = 0
    tot_misp_mon = 0           # from monitor mispredictions
    tot_misp_base = 0.0        # from non-monitor .out, summed only where present
    tot_instr_base = 0
    tot_new_misp = {name: 0.0 for name in swaps}

    for r in rows:
        vals = [
            r["trace"],
            f"{r['mpki']:.3f}",
            "--" if r["base_mpki"] is None else f"{r['base_mpki']:.3f}",
            f"{r['bim_prov_pct']:.1f}",
            f"{100*r['bim_acc']:.1f}",
        ]
        for name in swaps:
            acc = r.get(f"{name}_acc")
            mpki = r.get(f"{name}_mpki")
            delta = r.get(f"{name}_delta")
            vals += [
                "--" if acc is None else f"{100*acc:.1f}",
                "--" if mpki is None else f"{mpki:.3f}",
                "--" if delta is None else f"{delta:+.3f}",
            ]
        print(fmt_row(vals))

        tot_instr += r["instr"]
        tot_misp_mon += r["mispred"]
        if r["base_mpki"] is not None:
            tot_misp_base += r["base_mpki"] * r["instr"] / 1000.0
            tot_instr_base += r["instr"]
        for name in swaps:
            mpki = r.get(f"{name}_mpki")
            if mpki is None:
                continue
            tot_new_misp[name] += mpki * r["instr"] / 1000.0

    if tot_instr:
        mon_mpki = tot_misp_mon / tot_instr * 1000.0
        base_mpki = (tot_misp_base / tot_instr_base * 1000.0
                     if tot_instr_base else None)
        ref_mpki = base_mpki if base_mpki is not None else mon_mpki
        print()
        print(f"AGG ({len(rows)} traces, {tot_instr} instr):")
        print(f"  monitor MPKI:      {mon_mpki:.3f}")
        if base_mpki is not None:
            print(f"  non-monitor MPKI:  {base_mpki:.3f}  "
                  f"(delta vs monitor {base_mpki - mon_mpki:+.3f})")
        ref_label = "non-monitor" if base_mpki is not None else "monitor"
        for name in swaps:
            new_mpki = tot_new_misp[name] / tot_instr * 1000.0
            print(f"  {name} swap MPKI:    {new_mpki:.3f}  "
                  f"(delta vs {ref_label} {new_mpki - ref_mpki:+.3f})")


if __name__ == "__main__":
    main()
