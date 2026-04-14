#!/usr/bin/env python3
"""
Fanout audit for HARCOM code.

Counts how many times each val/reg/arr variable is read (rvalue) vs written (lvalue),
and checks against any .fanout(hard<N>{}) declarations.

Usage:
    python3 scripts/fanout_audit.py predictors/custom/TageAhead.hpp
    python3 scripts/fanout_audit.py predictors/custom/TageAhead.hpp --function predict1
"""

import re
import sys
import argparse
from collections import defaultdict
from pathlib import Path


# Matches HARCOM type declarations:
#   val<N> name, val<N, T> name, reg<N> name, arr<T, N> name, auto name = ...
DECL_RE = re.compile(
    r"""(?:
        (?:val|reg|arr)\s*<[^>]+>\s+   # val<N>, reg<N>, arr<T,N>
        |auto\s+                        # auto (often deduced to val/reg)
    )
    (\w+)                               # capture variable name
    """,
    re.VERBOSE,
)

# Matches .fanout(hard<N>{}) or .fanout(hard<N>{...})
FANOUT_RE = re.compile(r"(\w+)\.fanout\(hard<(\d+)>\{[^}]*\}\)")

# Matches arr fanout: name.fanout(hard<N>{})
ARR_FANOUT_RE = re.compile(r"(\w+)\.fanout\(hard<(\d+)>\{[^}]*\}\)")

# C++ identifier
IDENT_RE = re.compile(r"\b([a-zA-Z_]\w*)\b")

# Matches lvalue assignment: name = expr (but not ==)
ASSIGN_RE = re.compile(r"\b(\w+)\s*(?:\[.*?\])?\s*=[^=]")

# Lines that are declarations (to skip counting the declared name as rvalue)
DECL_LINE_RE = re.compile(
    r"^\s*(?:(?:val|reg|arr|auto)\s*(?:<[^>]+>)?\s+\w+|"
    r"(?:static\s+)?constexpr\s+)"
)


def extract_functions(source: str) -> dict[str, str]:
    """Extract function bodies by name (very rough brace matching)."""
    # Match function signatures like: val<1> predict1(...) {
    func_re = re.compile(
        r"(?:val<\d+>|void)\s+(\w+)\s*\([^)]*\)\s*\{", re.DOTALL
    )
    functions = {}
    for m in func_re.finditer(source):
        name = m.group(1)
        start = m.end() - 1  # the opening {
        depth = 0
        end = start
        for i in range(start, len(source)):
            if source[i] == "{":
                depth += 1
            elif source[i] == "}":
                depth -= 1
                if depth == 0:
                    end = i + 1
                    break
        functions[name] = source[start:end]
    return functions


def strip_comments(source: str) -> str:
    """Remove C/C++ comments."""
    # Remove // line comments
    source = re.sub(r"//.*?$", "", source, flags=re.MULTILINE)
    # Remove /* block comments */
    source = re.sub(r"/\*.*?\*/", "", source, flags=re.DOTALL)
    return source


def strip_strings(source: str) -> str:
    """Remove string literals."""
    return re.sub(r'"[^"]*"', '""', source)


def find_declarations(source: str) -> set[str]:
    """Find all variable names declared as val/reg/arr/auto."""
    names = set()
    for m in DECL_RE.finditer(source):
        names.add(m.group(1))

    # Also pick up C-array declarations: reg<N> name[SIZE]
    c_arr_re = re.compile(r"(?:val|reg)\s*<[^>]+>\s+(\w+)\s*\[")
    for m in c_arr_re.finditer(source):
        names.add(m.group(1))

    return names


def count_rw(source: str, tracked_names: set[str]) -> dict:
    """
    For each tracked variable, count:
      - writes: appears on LHS of = (not ==)
      - reads: appears anywhere else (as rvalue)
      - fanout: declared fanout value from .fanout(hard<N>{})
    """
    stats = {}
    for name in tracked_names:
        stats[name] = {"writes": 0, "reads": 0, "fanout": 0, "fanout_line": None}

    lines = source.split("\n")

    for lineno, line in enumerate(lines, 1):
        stripped = line.strip()
        if not stripped or stripped.startswith("//"):
            continue

        # Check for fanout declarations
        for m in FANOUT_RE.finditer(line):
            fname = m.group(1)
            fval = int(m.group(2))
            if fname in stats:
                stats[fname]["fanout"] += fval
                stats[fname]["fanout_line"] = lineno

        # Find all identifiers in this line
        idents_in_line = IDENT_RE.findall(line)

        # Check if this line has an assignment to a tracked variable
        for m in ASSIGN_RE.finditer(line):
            lhs = m.group(1)
            if lhs in stats:
                stats[lhs]["writes"] += 1

        # Count reads: every occurrence that isn't on the LHS of assignment,
        # isn't in a declaration, and isn't part of .fanout()
        # Simple heuristic: count all occurrences, subtract writes and declarations
        for name in tracked_names:
            # Count raw occurrences using word boundary
            occurrences = len(re.findall(rf"\b{re.escape(name)}\b", line))
            if occurrences == 0:
                continue

            # Subtract if this is a declaration line for this variable
            if re.search(
                rf"(?:val|reg|arr|auto)\s*(?:<[^>]+>)?\s+{re.escape(name)}\b",
                line,
            ):
                occurrences -= 1  # the declaration itself

            # Subtract writes (LHS of =)
            write_count = len(
                re.findall(
                    rf"\b{re.escape(name)}\b\s*(?:\[[^\]]*\])?\s*=[^=]", line
                )
            )
            reads = occurrences - write_count

            # Subtract fanout self-reference (name.fanout)
            fanout_self = len(
                re.findall(rf"\b{re.escape(name)}\b\.fanout", line)
            )
            reads -= fanout_self

            if reads > 0:
                stats[name]["reads"] += reads

    return stats


def audit(source: str, scope_name: str = "global"):
    """Run the audit on a source string."""
    clean = strip_strings(strip_comments(source))
    names = find_declarations(clean)

    if not names:
        print(f"  [{scope_name}] No tracked variables found.")
        return

    stats = count_rw(clean, names)

    # Filter to variables that are actually used
    used = {
        k: v
        for k, v in stats.items()
        if v["reads"] > 0 or v["writes"] > 0
    }

    if not used:
        print(f"  [{scope_name}] No tracked variable usage found.")
        return

    # Report
    print(f"\n{'='*70}")
    print(f"  Fanout audit: {scope_name}")
    print(f"{'='*70}")
    print(f"  {'Variable':<25} {'Reads':>6} {'Writes':>7} {'Fanout':>7}  Status")
    print(f"  {'-'*25} {'-'*6} {'-'*7} {'-'*7}  {'-'*15}")

    issues = []
    for name in sorted(used.keys()):
        s = used[name]
        reads = s["reads"]
        writes = s["writes"]
        fanout = s["fanout"]

        if reads <= 1:
            status = "ok (single read)"
        elif fanout == 0 and reads > 1:
            status = f"MISSING fanout({reads})"
            issues.append((name, reads, fanout))
        elif fanout < reads:
            status = f"LOW fanout({fanout} < {reads})"
            issues.append((name, reads, fanout))
        elif fanout > reads:
            status = f"HIGH fanout({fanout} > {reads})"
        else:
            status = "ok"

        print(f"  {name:<25} {reads:>6} {writes:>7} {fanout:>7}  {status}")

    if issues:
        print(f"\n  ** {len(issues)} issue(s) found **")
        for name, reads, fanout in issues:
            if fanout == 0:
                print(f"     {name}: needs .fanout(hard<{reads}>{{}})")
            else:
                print(
                    f"     {name}: has fanout({fanout}), needs fanout({reads})"
                )
    else:
        print(f"\n  All fanouts look correct.")


def main():
    parser = argparse.ArgumentParser(description="HARCOM fanout audit")
    parser.add_argument("file", help="C++ source file to audit")
    parser.add_argument(
        "--function", "-f", help="Audit a specific function only"
    )
    parser.add_argument(
        "--all-functions",
        "-a",
        action="store_true",
        help="Audit each function separately",
    )
    args = parser.parse_args()

    source = Path(args.file).read_text()

    if args.function:
        functions = extract_functions(source)
        if args.function not in functions:
            print(f"Function '{args.function}' not found.")
            print(f"Available: {', '.join(functions.keys())}")
            sys.exit(1)
        audit(functions[args.function], args.function)
    elif args.all_functions:
        functions = extract_functions(source)
        if not functions:
            print("No functions found.")
            sys.exit(1)
        for fname, fbody in functions.items():
            audit(fbody, fname)
    else:
        audit(source, Path(args.file).name)


if __name__ == "__main__":
    main()
