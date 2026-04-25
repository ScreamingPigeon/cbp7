#!/usr/bin/env python3
"""Static check: find named vals that are read without fanout() or fo1() declaration.

Under CHECK_FANOUT, every val must have either fanout(hard<N>) or fo1() called
before any lvalue read. This script finds likely violations by scanning for
named val/auto declarations that are subsequently used without a .fanout() or
.fo1() call.

Limitations:
  - Regex-based, not a real parser. May have false positives/negatives.
  - Does not track if constexpr branches (may report both sides).
  - Does not track scope nesting perfectly.
  - Ignores rvalue temporaries from operators (those are harder to track).
"""

import re
import sys
from collections import defaultdict


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "predictors/custom/TageAhead.hpp"
    lines = open(path).readlines()

    # Track declarations and their usage
    # decl_line -> var_name
    declarations = {}  # (name, line) -> type
    has_fanout = set()  # (name, line) set of vars that have fanout/fo1
    reads = defaultdict(list)  # (name, decl_line) -> [read_lines]

    # Patterns
    # Match: val<...> name = ...; or auto name = ...;
    decl_pat = re.compile(
        r'(?:val<[^>]+>|auto)\s+(\w+)\s*='
    )
    # Match: varname.fanout( or varname.fo1()
    fanout_pat = re.compile(r'(\w+)\.fanout\(')
    fo1_pat = re.compile(r'(\w+)\.fo1\(\)')
    # Match: varname used as lvalue (not after . or ->)
    # We look for the variable name used in expressions

    # First pass: collect all declarations
    scope_stack = []  # track brace depth roughly
    for i, line in enumerate(lines, 1):
        # Skip comments and preprocessor
        stripped = line.lstrip()
        if stripped.startswith('//') or stripped.startswith('#'):
            continue

        for m in decl_pat.finditer(line):
            name = m.group(1)
            # Skip common non-val names
            if name in ('constexpr', 'const', 'static', 'return', 'void',
                        'int', 'u64', 'i64', 'bool', 'size_t'):
                continue
            declarations[(name, i)] = line.strip()

    # Second pass: find fanout/fo1 calls
    for i, line in enumerate(lines, 1):
        stripped = line.lstrip()
        if stripped.startswith('//') or stripped.startswith('#'):
            continue
        for m in fanout_pat.finditer(line):
            has_fanout.add(m.group(1))
        for m in fo1_pat.finditer(line):
            has_fanout.add(m.group(1))

    # Find vars that are declared but never get fanout/fo1
    # AND are used after declaration (read at least once)
    missing = []
    for (name, decl_line), decl_text in sorted(declarations.items(), key=lambda x: x[0][1]):
        if name in has_fanout:
            continue
        # Check if this var is used after its declaration line
        # Look for the name appearing in subsequent lines (not as part of another word)
        use_pat = re.compile(r'(?<![.\w])' + re.escape(name) + r'(?!\w*\s*[=])')
        used = False
        for j in range(decl_line, min(decl_line + 80, len(lines) + 1)):
            line = lines[j - 1]
            # Skip the declaration line itself
            if j == decl_line:
                continue
            if use_pat.search(line):
                # Make sure it's not just a comment
                code_part = line.split('//')[0]
                if use_pat.search(code_part):
                    used = True
                    break

        if used:
            missing.append((decl_line, name, decl_text))

    if not missing:
        print("No missing fanout/fo1 declarations found.")
        return

    print(f"Found {len(missing)} vals without fanout/fo1 that are read:")
    print()
    for line, name, text in missing:
        print(f"  line {line}: {name}")
        print(f"    {text}")
        print()


if __name__ == "__main__":
    main()
