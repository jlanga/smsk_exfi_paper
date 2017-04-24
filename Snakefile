shell.prefix("set -euo pipefail;")
configfile: "config.yaml"



# Some variables
dna_pe = [
    sample_name for sample_name in config["samples"]
        if config["samples"][sample_name]["type"] == "PE" and 
        config["samples"][sample_name]["molecule"] == "dna"
]
THREAD_LIMIT = 4


# Read subsnakefiles
snakefiles = "bin/snakefiles/"
include: snakefiles + "folders.py"
include: snakefiles + "generic.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "exons.py"
include: snakefiles + "bwa.py"
include: snakefiles + "dist.py"



rule all:
    input:
        ## raw
        # expand(
        #    raw + "{sample}_{end}.fq.gz",
        #    sample = dna_pe,
        #    end = "1 2".split(" ")
        #),
        #raw + "transcriptome.fa.fai",
        ## exfi
        #expand(
        #    exons + "k{kmer}_l{levels}_m{size}.bloom",
        #    kmer = config["exons"]["kmer"],
        #    levels = config["exons"]["levels"],
        #    size = config["exons"]["size"]
        #),
        #exons + "raw.fa.fai",
        #exons + "filtered_by_length.fa.fai",
        #exons + "filtered_by_extensibility.fa.fai"
        # bwa
        #bwa + "exome",
        #bwa + "transcriptome",
        #bwa + "genome",
        expand(
            bwa + "{input}_vs_{reference}.bam.bai",
            input = ["raw", "filtered_by_length", "filtered_by_extensibility"],
            reference = ["exome", "transcriptome", "genome"]
        ),
        ## dist
        dist + "exon_histogram.pdf",
        dist + "exon_density.pdf",
