#!/usr/bin/env bash

cd "$(dirname "$0")/../.."
source_files=($(git ls-files '*.R'))
grep -hE '\b(require|library)\(\w+\)' "${source_files[@]}" | \
    sed '/^[[:space:]]*#/d' | \
    sed -E 's/.*\(([[:alnum:]]*)\).*/\1/' | \
    sort -uf \
    > DEPENDS
