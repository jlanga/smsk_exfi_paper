shell.prefix("set -euo pipefail;")
configfile: "config.yaml"

dna_pe = [
    sample_name for sample_name in config["samples"]
        if config["samples"][sample_name]["type"] == "PE" and 
        config["samples"][sample_name]["molecule"] == "dna"
]


snakefiles = "bin/snakefiles/"

include: snakefiles + "folders.py"
#include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "exons.py"
include: snakefiles + "bwa.py"



rule all:
    input:
        # raw
        expand(
            raw + "{sample}_{end}.fq.gz",
            sample = dna_pe,
            end = "1 2".split(" ")
        ),
        raw + "transcriptome.fa.fai",
        # exfi
        expand(
            exons + "k{kmer}_l{levels}_m{size}.bloom",
            kmer = config["exons"]["kmer"],
            levels = config["exons"]["levels"],
            size = config["exons"]["size"]
        ),
        exons + "exons_raw.fa",
        exons + "exons_filtered_by_length.fa",
        