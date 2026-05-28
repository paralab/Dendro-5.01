#!/usr/bin/env bash
#
# Full-tree clang-format of all FIRST-PARTY C/C++ sources, excluding vendored
# code. Use this to (re)establish the formatting baseline. For incremental
# (changed-lines-only) formatting use scripts/format-changed.sh instead.
#
#     scripts/format-all.sh           # reformat the whole first-party tree
#     scripts/format-all.sh --check   # exit non-zero if anything is unformatted
#
# Keep VENDORED in sync with .clang-format-ignore.

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

# Vendored paths (anchored, regex) — never reformat these.
VENDORED='^(IO/zlib/|lib/BHaHAHA/|IO/vtk/include/(json\.hpp|json_fwd\.hpp|base\.h|cencode\.h|cdecode\.h)$|include/dollar\.hpp$)'

mapfile -t files < <(
    git ls-files '*.c' '*.h' '*.hpp' '*.cpp' '*.cc' '*.cxx' '*.tcc' '*.cu' '*.cuh' \
        | grep -vE "$VENDORED"
)

if [ "${1:-}" = "--check" ]; then
    rc=0
    for f in "${files[@]}"; do
        diff -q <(clang-format "$f") "$f" >/dev/null 2>&1 || { echo "unformatted: $f"; rc=1; }
    done
    [ "$rc" -eq 0 ] && echo "All ${#files[@]} first-party files are clang-format clean."
    exit "$rc"
fi

printf '%s\n' "${files[@]}" | xargs clang-format -i
echo "Reformatted ${#files[@]} first-party files (vendored excluded). Review: git diff"
