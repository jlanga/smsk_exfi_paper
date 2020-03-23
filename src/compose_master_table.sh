#!/usr/bin/env bash

set -euo pipefail

# Get BWA stats
parallel --tag \
    tail -n+2 {} \
::: results/bwa/*.stats.tsv \
> bwa_initial.stats.tsv

parallel --tag \
    tail -n+2 {} \
::: results/bwa_review/*.stats.tsv \
> bwa_review.stats.tsv


# Get pr
parallel --tag \
    grep -v ^type {} \
:::  results/pr/*.tsv > pr_initial.tsv

parallel --tag \
    grep -v ^type {} \
:::  results/pr_review/*.tsv > pr_review.tsv



# exfi
# get fprs
grep FPR results/exfi/*.bf.log \
| cut -f 1,3 -d : \
| tr -d ' ' \
| tr ':' '\t' \
| sed 's/.log//' \
> exfi.fpr.tsv

# get time and memory of build
parallel --tag \
    cat {} \
::: results/exfi/*.bf.bmk \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> exfi.build.bmk

# get time and memory of predict
parallel \
    --tag \
    cat {} \
::: results/exfi/*.gfa.bmk \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> exfi.predict.bmk

# Get insertion rates
parallel --tag \
    'grep transcriptome {} | tail -1 | cut -f 5' \
::: results/exfi/*.bf.log > exfi.ir.tsv




# Chopstitch
# Get BF sizes, hash functions, actual fprs
echo -e "file\tfpr1\tfpr2\tbf1_bitsize\tbf2_bitsize\thash_fn1\thash_fn2" > chopstitch_build_data_initial.tsv
for i in  results/chopstitch/*.bf.log ; do
    paste \
        <(echo "$i") \
        <(grep 'Primary BF actual fpr:' "$i" | grep -oE '([0-9.]+)') \
        <(grep 'Secondary BF actual fpr:' "$i" | grep -oE '([0-9.]+)') \
        <(grep 'Primary BF bits:' "$i" | grep -oE '([0-9]+)' | head -1) \
        <(grep 'Secondary BF bits:' "$i" | grep -oE '([0-9]+)' | head -1) \
        <(grep 'Primary BF hashes:' "$i" | grep -oE '([0-9.]+)') \
        <(grep 'Secondary BF hashes:' "$i" | grep -oE '([0-9.]+)')
done >> chopstitch_build_data_initial.tsv

echo -e "file\tfpr1\tfpr2\tbf1_bitsize\tbf2_bitsize\thash_fn1\thash_fn2" > chopstitch_build_data_review.tsv
for i in  results/chopstitch_review/*.bf.log ; do
    paste \
        <(echo "$i") \
        <(grep 'Primary BF actual fpr:' "$i" | grep -oE '([0-9.]+)') \
        <(grep 'Secondary BF actual fpr:' "$i" | grep -oE '([0-9.]+)') \
        <(grep 'Primary BF bits:' "$i" | grep -oE '([0-9]+)' | head -1) \
        <(grep 'Secondary BF bits:' "$i" | grep -oE '([0-9]+)' | head -1) \
        <(grep 'Primary BF hashes:' "$i" | grep -oE '([0-9.]+)') \
        <(grep 'Secondary BF hashes:' "$i" | grep -oE '([0-9.]+)')
done >> chopstitch_build_data_review.tsv



# Get Build time - initial round
parallel --tag \
    cat {} \
::: results/chopstitch/*.bf.bmk  \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> chopstitch_initial.bf.bmk

# Build time - reviewer
parallel --tag \
    cat {} \
::: results/chopstitch_review/*.bf.bmk  \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> chopstitch_review.bf.bmk

# Get predict time -- initial
parallel --tag \
    cat {} \
::: results/chopstitch/*.exons.bmk \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> chopstitch_initial.predict.bmk

# Get predict time -- review
parallel --tag \
    cat {} \
::: results/chopstitch_review/*.exons.bmk \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> chopstitch_review.predict.bmk



# Gmap
# Get Build time
parallel \
    --tag cat {} \
::: results/gmap/*.build.bmk \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> gmap.build.bmk

# Get predict time
parallel --tag \
    cat {} \
::: results/gmap/*.gff3.bmk \
| grep -v h:m:s \
| sort \
| cut -f 1,3,4 \
| sed 's/.bmk//' \
> gmap.predict.bmk



Rscript src/compose_master_table.R
rm pr_initial.tsv pr_review.tsv
rm exfi.fpr.tsv exfi.build.bmk exfi.predict.bmk exfi.ir.tsv
rm chopstitch_build_data_*.tsv chopstitch_*.bf.bmk chopstitch_*.predict.bmk
rm gmap.build.bmk gmap.predict.bmk
rm bwa_*.stats.tsv