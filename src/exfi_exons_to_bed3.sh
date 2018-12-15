#!/usr/bin/env bash

set -euo pipefail

fasta=$1

grep ^'>' "$fasta" \
| tr -d '>' \
| tr '-' '\t' \
| tr ':' '\t' \
| sort -k1,1 -k2,2n
