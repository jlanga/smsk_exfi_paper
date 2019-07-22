#!/usr/bin/env bash

compute_stats() {

    n=$(samtools view "$1" | cut -f 1 | sort -u | wc -l)
    unmapped=$(samtools view "$1" | awk '$3 == "*"' | cut -f 1 | sort -u | wc -l)
    mapped=$(samtools view "$1" | awk '$3 != "*"'| cut -f 1 | sort -u | wc -l)
    # perfect=$(samtools view "$1" | awk '$6 ~ /^[0-9]+M$/' | cut -f 1 | sort -u | wc -l)
    perfect=$(samtools view "$1" | awk '$6 ~ /^[0-9MID]+$/' | cut -f 1 | sort -u | wc -l)
    unique=$(samtools view "$1" | cut -f 1 | sort | uniq -c | awk '$1 == 1' | wc -l)
    multi=$(samtools view "$1" | cut -f 1 | sort | uniq -c | awk '$1 > 1' | wc -l)
    mapped_transcripts=$(samtools view "$1" | cut -f 1 | cut -f 1-5 -d _  | cut -f 1 -d : | uniq | sort -u | wc -l)
    transcripts=$(grep -c ^">" "$2")

    echo -e \
        "exons\\tunmapped\\tmapped\\tperfect_match\\tuniquely_mapped\\tmultimapped\\ttranscripts\\tmapped_transcripts"
    echo -e \
        "$n\\t$unmapped\\t$mapped\\t$perfect\\t$unique\\t$multi\\t$transcripts\\t$mapped_transcripts"

}



compute_stats "$1" "$2"
