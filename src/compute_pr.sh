#!/usr/bin/env bash

set -euo pipefail

compute_tp(){
    predicted=$1
    observed=$2
    identity=$3
    bedtools intersect -a "$predicted" -b "$observed" -f "$identity" -r \
    | wc -l
}

compute_fp(){
    predicted=$1
    observed=$2
    identity=$3
    bedtools intersect -a "$predicted" -b "$observed" -f "$identity" -r -v \
    | wc -l
}

compute_fn(){
    predicted=$1
    observed=$2
    identity=$3
    bedtools intersect -a "$observed" -b "$predicted" -f "$identity" -r -v \
    | wc -l
}

compute_pr(){
    predicted=$1
    observed=$2
    identity=$3

    tp=$(compute_tp "$predicted" "$observed" "$identity")
    fp=$(compute_fp "$predicted" "$observed" "$identity")
    fn=$(compute_fn "$predicted" "$observed" "$identity")

    precision=$(echo "$tp / ($tp + $fp)" | bc -l)
    recall=$(echo "$tp / ($tp + $fn)" | bc -l)
    f1=$(echo "2 * $precision * $recall / ($precision + $recall)" | bc -l)
    echo -e "identity\\ttp\\tfp\\tfn\\tprecision\\trecall\\tf1"
    echo -e "$identity\\t$tp\\t$fp\\t$fn\\t$precision\\t$recall\\t$f1"

}

predicted=$1
observed=$2
identity=$3

compute_pr "$predicted" "$observed" "$identity"
