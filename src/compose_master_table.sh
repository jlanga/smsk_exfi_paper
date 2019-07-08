#!/usr/bin/env bash

set -euo pipefail

# Get BWA stats
parallel --tag tail -n+2 {} ::: results/bwa/*.stats.tsv > bwa.stats.tsv


# Get pr
parallel --tag grep -v ^type {} ::: results/pr/*.tsv > pr_all.tsv



# exfi
# get fprs
grep FPR results/exfi/*.bf.log | cut -f 1,3 -d : | tr -d ' ' | tr ':' '\t' | sed 's/.log//'> exfi.fpr.tsv

# get time and memory of build
parallel --tag cat {} ::: results/exfi/*.bf.bmk | grep -v h:m:s | sort | cut -f 1,3,4 | sed 's/.bmk//' > exfi.build.bmk

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
parallel \
    --tag \
    'grep transcriptome {} | tail -1 | cut -f 5' \
::: results/exfi/*.bf.log > exfi.ir.tsv




# Chopstitch
# Get BF sizes, hash functions, actual fprs
echo -e "file\tfpr1\tfpr2\tbf1_bitsize\tbf2_bitsize\thash_fn1\thash_fn2" > chopstitch_build_data.tsv
for i in  results/chopstitch/*.bf.log; do
    paste \
        <(echo $i) \
        <(grep 'Primary BF actual fpr:' $i | grep -oE '([0-9.]+)') \
        <(grep 'Secondary BF actual fpr:' $i | grep -oE '([0-9.]+)') \
        <(grep 'Primary BF bits:' $i | grep -oE '([0-9]+)' | head -1) \
        <(grep 'Secondary BF bits:' $i | grep -oE '([0-9]+)' | head -1) \
        <(grep 'Primary BF hashes:' $i | grep -oE '([0-9.]+)') \
        <(grep 'Secondary BF hashes:' $i | grep -oE '([0-9.]+)')
done >> chopstitch_build_data.tsv

# Get Build time
parallel --tag cat {} ::: results/chopstitch/*.bf.bmk | grep -v h:m:s | sort | cut -f 1,3,4 | sed 's/.bmk//' > chopstitch.bf.bmk

# Get predict time
parallel --tag cat {} ::: results/chopstitch/*.exons.bmk | grep -v h:m:s | sort | cut -f 1,3,4 | sed 's/.bmk//' > chopstitch.predict.bmk



# Gmap
# Get Build time
parallel --tag cat {} ::: results/gmap/*.build.bmk | grep -v h:m:s | sort | cut -f 1,3,4 | sed 's/.bmk//' > gmap.build.bmk

# Get predict time
parallel --tag cat {} ::: results/gmap/*.gff3.bmk | grep -v h:m:s | sort | cut -f 1,3,4 | sed 's/.bmk//' > gmap.predict.bmk



Rscript src/compose_master_table.R
rm pr_all.tsv
rm exfi.fpr.tsv exfi.build.bmk exfi.predict.bmk exfi.ir.tsv
rm chopstitch_build_data.tsv chopstitch.bf.bmk
rm gmap.build.bmk gmap.predict.bmk
rm bwa.stats.tsv chopstitch.predict.bmk
