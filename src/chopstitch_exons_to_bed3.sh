#!/usr/bin/env bash

set -euo pipefail

fasta=$1

grep ^'>' "$fasta" \
| tr -d '>' \
| tr '_' '\t' \
| bedtools sort
