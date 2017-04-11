#!/usr/bin/env bash
set -euo pipefail

mkdir -p data/reads/

gzip \
    --decompress \
    --keep \
    data/reads/genome.fa.gz
samtools faidx data/reads/genome.fa

>&2 echo "Simulating dna1"

wgsim \
    -N 10000 \
    -1 100 \
    -2 100 \
    -S 1 \
    data/reads/genome.fa \
    >(pigz --best > data/reads/dna1_1.fq.gz) \
    >(pigz --best > data/reads/dna1_2.fq.gz) \
> >(pigz --best > data/reads/dna1.log.gz) &

>&2 echo "Simulating dna2"

wgsim \
    -N 10000 \
    -1 100 \
    -2 100 \
    -S 2 \
    data/reads/genome.fa \
    >(pigz --best > data/reads/dna2_1.fq.gz) \
    >(pigz --best > data/reads/dna2_2.fq.gz) \
> >(pigz --best > data/reads/dna2.log.gz) &

wait
