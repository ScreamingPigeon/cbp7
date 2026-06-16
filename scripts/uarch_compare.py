#!/usr/bin/env python3
"""
Compute uarch-metric distributions across the v5 monitor full sweep, then
compare the eight outlier traces against the corpus distribution.
"""

import csv
import os
import re
import statistics
import sys
from pathlib import Path

CORPUS_DIR = Path("full_out/monitor_all_v5")
OUTLIERS = [
    "web_99", "int_117", "snappy-unit_431", "505-mcf-1_45923",
    "gmsh-1.30841_0", "abc-3.3561_0", "gap-sssp-usroad.287_0",
    "541-leela-1_110755",
]


def grab_re(s, pat, group=1, conv=float, default=None):
    m = re.search(pat, s, re.DOTALL)
    if not m:
        return default
    try:
        return conv(m.group(group))
    except (ValueError, IndexError):
        return default


def extract(path):
    if not path.exists():
        return None
    s = path.read_text()
    d = {}

    # Top-line
    d["mpki"]            = grab_re(s, r"MPKI:\s+([\d.]+)")
    d["instr"]           = grab_re(s, r"Instructions:\s+(\d+)", conv=int)
    d["branches"]        = grab_re(s, r"Branches:\s+(\d+)", conv=int)
    d["mispred"]         = grab_re(s, r"Mispredictions:\s+(\d+)", conv=int)
    d["misp_pct"]        = grab_re(s, r"Mispredictions:\s+\d+\s+\(([\d.]+)%\)")
    d["extra_cycle_pct"] = grab_re(s, r"Extra cycles:\s+\d+\s+\(([\d.]+)%\)")
    d["true_blocks_pct"] = grab_re(s, r"True blocks:\s+\d+\s+\(([\d.]+)%\)")

    # Block structure
    d["avg_instr_per_block"]   = grab_re(s, r"Avg instr/block:\s+([\d.]+)")
    d["avg_branches_per_block"] = grab_re(s, r"Avg branches/block:\s+([\d.]+)")
    d["avg_nb_misp_blocks"]    = grab_re(
        s, r"Avg num_branch in misp blocks:\s+([\d.]+)")
    d["avg_nb_correct_blocks"] = grab_re(
        s, r"Avg num_branch in misp blocks:\s+[\d.]+\s+\(vs\s+([\d.]+)")

    # Branch density (derived)
    if d.get("instr") and d.get("branches"):
        d["br_per_kinstr"] = 1000.0 * d["branches"] / d["instr"]

    # Provider Distribution: Bimodal + T0..T14
    prov = re.search(r"Provider Distribution:(.*?)\n\n", s, re.DOTALL)
    if prov:
        ptext = prov.group(1)
        bim = re.search(r"Bimodal\s+\|\s+(\d+)\s+\|\s+(\d+)\s+\|\s+([\d.]+)%", ptext)
        if bim:
            d["bim_prov"] = int(bim.group(1))
            d["bim_acc"]  = float(bim.group(3))
            if d.get("branches"):
                d["bim_share_pct"] = 100.0 * d["bim_prov"] / d["branches"]
        tage_rows = re.findall(
            r"T(\d+)\s+\|\s+(\d+)\s+\|\s+(\d+)\s+\|\s+([\d.]+)%", ptext)
        if tage_rows:
            tage_total_prov    = sum(int(r[1]) for r in tage_rows)
            tage_total_correct = sum(int(r[2]) for r in tage_rows)
            d["tage_share_pct"] = 100.0 * tage_total_prov / d["branches"] if d.get("branches") else None
            d["tage_acc"] = 100.0 * tage_total_correct / tage_total_prov if tage_total_prov else None
            # T0
            t0 = [r for r in tage_rows if r[0] == "0"]
            if t0:
                d["t0_share_pct"] = 100.0 * int(t0[0][1]) / d["branches"] if d.get("branches") else None
                d["t0_acc"]       = float(t0[0][3])
            # Worst-accuracy table among those that provide >0.5% of branches
            sig = [r for r in tage_rows
                   if d.get("branches") and 100.0 * int(r[1]) / d["branches"] > 0.5]
            if sig:
                worst = min(sig, key=lambda r: float(r[3]))
                d["worst_table_acc"] = float(worst[3])
                d["worst_table_idx"] = int(worst[0])

    # Overall TAGE accuracy line
    d["overall_tage_acc"] = grab_re(s, r"Overall TAGE accuracy:\s+([\d.]+)%")

    # Provider Source Breakdown (per-branch routing & per-source accuracy)
    src = re.search(r"Provider Source Breakdown:(.*?)\n\n", s, re.DOTALL)
    if src:
        st = src.group(1)
        for label, key in [
            ("No tag match \\(FB only\\)",       "src_no_tag"),
            ("Sec-tag rejected \\(FB\\)",        "src_sec_rej"),
            ("TAGE provider, meta→alt",          "src_meta_alt"),
            ("TAGE provider, meta→pri",          "src_meta_pri"),
            ("TAGE provider, no meta",           "src_no_meta"),
        ]:
            m = re.search(rf"{label}:\s+(\d+)\s+\(([\d.]+)%\)\s+acc=([\d.]+)%", st)
            if m:
                d[f"{key}_pct"] = float(m.group(2))
                d[f"{key}_acc"] = float(m.group(3))

    # Sec-tag counterfactual
    d["sec_tag_cf_delta_pct"] = grab_re(
        s,
        r"Sec-tag counterfactual.*?Delta \(TAGE - FB\):.*?\(([+-][\d.]+)%\)")

    # Meta override
    d["meta_override_acc"] = grab_re(s, r"Meta Override:.*?Correct:\s+\d+\s+\(([\d.]+)%\)")

    # Block-mispredict Cause Attribution (per-block)
    cause = re.search(r"Block-mispredict Cause Attribution:(.*?)\n\n", s, re.DOTALL)
    if cause:
        ct = cause.group(1)
        for label, key in [
            (r"No tag match \(FB only\)",       "cause_no_tag"),
            (r"Sec-tag rejected \(FB\)",        "cause_sec_rej"),
            (r"TAGE provider, meta->alt",       "cause_meta_alt"),
            (r"TAGE provider, meta->pri",       "cause_meta_pri"),
            (r"TAGE provider, no meta",         "cause_no_meta"),
            (r"Phantom slot \(>=num_br\)",      "cause_phantom"),
            (r"Target miss \(BTB/next-PC\)",    "cause_target"),
            (r"FB-caused total",                "cause_fb_total"),
        ]:
            m = re.search(rf"{label}:\s+\d+\s+\(([\d.]+)%\)", ct)
            if m:
                d[key + "_pct"] = float(m.group(1))

    # Allocation
    alloc = re.search(r"Allocation:(.*?)\n  Provider\s+\|", s, re.DOTALL)
    if alloc:
        at = alloc.group(1)
        d["alloc_success_pct"] = grab_re(at, r"Success:\s+\d+\s+\(([\d.]+)%\)")
        att = grab_re(at, r"Attempts:\s+(\d+)", conv=int)
        blocked = grab_re(at, r"Blocked:\s+(\d+)", conv=int)
        sibling = grab_re(at, r"Sibling skips:\s+(\d+)", conv=int)
        if att:
            if blocked is not None:
                d["alloc_blocked_pct"] = 100.0 * blocked / att
            if sibling is not None:
                d["alloc_sibling_skip_pct"] = 100.0 * sibling / att
    d["unique_branch_pcs"] = grab_re(s, r"Unique branch PCs:\s+(\d+)", conv=int)

    # Bias distribution
    d["pcs_easy_bias"] = grab_re(s, r"Easy \(>=90% biased\):\s+(\d+)", conv=int)
    d["pcs_hard_bias"] = grab_re(s, r"Hard \(<60% biased\):\s+(\d+)", conv=int)

    # Tag collision rate
    d["tag_collision_pct"] = grab_re(
        s, r"Tag Collision Rate.*?Total checks:\s+\d+\s+Collisions:\s+\d+\s+\(([\d.]+)%\)")

    # Misprediction Direction vs Boundary
    d["misp_dir_pct"]      = grab_re(s, r"Direction\s+\(wrong.*?\(([\d.]+)%\)")
    d["misp_boundary_pct"] = grab_re(s, r"Boundary\s+\(.*?\(([\d.]+)%\)")

    return d


def dist(values, percentiles=(10, 25, 50, 75, 90, 99)):
    v = sorted(x for x in values if x is not None)
    if not v:
        return None
    n = len(v)
    return {
        "n": n,
        "mean":   statistics.mean(v),
        "median": v[n // 2],
        **{f"p{p}": v[min(n - 1, int(p * n / 100))] for p in percentiles},
        "min": v[0],
        "max": v[-1],
    }


def rank_of(values, target):
    if target is None:
        return None
    v = sorted(x for x in values if x is not None)
    if not v:
        return None
    below = sum(1 for x in v if x < target)
    return 100.0 * below / len(v)


def main():
    summaries = list(CORPUS_DIR.glob("*_summary.txt"))
    if not summaries:
        sys.exit(f"No summary files in {CORPUS_DIR}")
    data = {}
    for p in summaries:
        t = p.name.replace("_summary.txt", "")
        d = extract(p)
        if d:
            data[t] = d
    print(f"Loaded {len(data)} traces from {CORPUS_DIR}")

    # All metrics gathered
    all_keys = sorted({k for d in data.values() for k in d.keys()})
    # Skip raw count fields for distribution display
    skip_dist = {"instr", "branches", "mispred", "bim_prov"}
    metrics = [k for k in all_keys if k not in skip_dist]

    print()
    print("CORPUS-WIDE DISTRIBUTION")
    print(f'{"metric":<28} {"n":>4} {"mean":>9} {"p10":>9} {"p25":>9} {"p50":>9} {"p75":>9} {"p90":>9} {"p99":>9}')
    print("-" * 110)
    for m in metrics:
        vals = [d.get(m) for d in data.values()]
        st = dist(vals)
        if st is None:
            continue
        print(f'{m:<28} {st["n"]:>4} {st["mean"]:>9.2f} {st["p10"]:>9.2f} {st["p25"]:>9.2f} {st["p50"]:>9.2f} {st["p75"]:>9.2f} {st["p90"]:>9.2f} {st["p99"]:>9.2f}')

    # Per-metric per-outlier with corpus percentile
    print()
    print("OUTLIER VALUES AND CORPUS PERCENTILE")
    header = f'{"metric":<28}' + "".join(f"  {t[:20]:>22}" for t in OUTLIERS)
    print(header)
    print("-" * len(header))
    for m in metrics:
        vals = [d.get(m) for d in data.values()]
        st = dist(vals)
        if st is None:
            continue
        row = [f"{m:<28}"]
        for t in OUTLIERS:
            v = data.get(t, {}).get(m)
            if v is None:
                row.append("                    --")
            else:
                pct = rank_of(vals, v)
                row.append(f"{v:>10.2f} (p{pct:>4.1f}) ")
        print("  ".join(row))


if __name__ == "__main__":
    main()
