#!/usr/bin/env bash

set -euo pipefail

fasta=$1

fields=$(head -1 "$fasta" | awk '{print $1}' | tr "_" "\n" | wc -l)

paste \
    <(grep ^">" "$fasta" \
        | tr -d ">" \
        | cut -f 1-$((fields - 2)) -d "_" \
    ) \
    <(grep ^">" "$fasta" \
        | cut -f $((fields - 1)),$fields -d "_" \
        | tr "_" "\t" \
    ) \
| bedtools sort \
| awk '{print $1 "\t" $2 - 1 "\t" $3}'
