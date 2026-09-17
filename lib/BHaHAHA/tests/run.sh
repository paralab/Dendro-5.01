#!/bin/sh
# Build and test the C library from both C and C++ callers.
set -eu
project_dir=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
test_dir=$(mktemp -d)
trap 'rm -rf "$test_dir"' EXIT HUP INT TERM
make -C "$project_dir" -j2
"${CXX:-c++}" -std=c++11 -Wall -Wextra -Werror -pedantic -I"$project_dir" \
  "$project_dir/tests/public_header.cpp" \
  -Wl,--whole-archive "$project_dir/libBHaHAHA.a" -Wl,--no-whole-archive \
  -lm -fopenmp -o "$test_dir/public_header"
"$test_dir/public_header"
"${CC:-cc}" -std=gnu99 -Wall -Wextra -Werror -I"$project_dir" \
  "$project_dir/tests/schwarzschild.c" "$project_dir/libBHaHAHA.a" \
  -lm -fopenmp -o "$test_dir/schwarzschild"
OMP_NUM_THREADS=2 "$test_dir/schwarzschild"
