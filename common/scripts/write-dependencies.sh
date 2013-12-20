#!/usr/bin/env bash

cd "$(dirname "$0")/../.."
grep -hE '\b(require|library)\b' {common,chip,rna}/scripts/*.R | \
    sed '/^[[:space:]]*#/d' | \
    sed -E 's/.*\(([[:alnum:]]*)\).*/\1/' | \
    sort -uf \
    > DEPENDS
