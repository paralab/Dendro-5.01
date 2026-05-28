#!/usr/bin/env bash
#
# Format only the lines changed relative to a base ref (default: master),
# using the repo .clang-format. Untouched lines are left alone.
#
# This is the jj-friendly entry point: jj does not run git hooks, so run this
# before squashing/pushing a branch.
#
#     scripts/format-changed.sh            # format vs master
#     scripts/format-changed.sh <base-ref> # format vs another base
#
# Review the result with: jj diff   (or git diff)

set -euo pipefail

base="${1:-master}"
exts="tcc,h,hpp,cpp,cc,c"

if ! command -v git-clang-format >/dev/null 2>&1; then
    echo "git-clang-format not found in PATH." >&2
    exit 1
fi

git clang-format "$base" --extensions "$exts"
echo "Formatted lines changed since '$base'. Review with: jj diff"
