shell.prefix("set -euo pipefail;")
configfile: "src/config.yaml"



# Some variables
dna_pe = [
    sample_name for sample_name in config["samples"]
        if config["samples"][sample_name]["type"] == "PE" and
        config["samples"][sample_name]["molecule"] == "dna"
]
THREAD_LIMIT = 64


# Read subsnakefiles
snakefiles = "src/snakefiles/"
include: snakefiles + "folders.py"
include: snakefiles + "generic.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "exfi.py"
include: snakefiles + "bwa_ref.py"
include: snakefiles + "dist.py"
include: snakefiles + "pr.py"



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
        #    exfi + "k{kmer}_l{levels}_m{size}.bloom",
        #    kmer = config["exons"]["kmer"],
        #    levels = config["exons"]["levels"],
        #    size = config["exons"]["size"]
        #),
        #exfi + "raw.fa.fai",
        #exfi + "filtered_by_length.fa.fai",
        #exfi + "filtered_by_extensibility.fa.fai"
        # bwa_ref
        #bwa + "exome",
        #bwa + "transcriptome",
        #bwa + "genome",
        # expand(
        #     bwa_ref + "{input}_vs_{reference}.bam.bai",
        #     input = ["exons"],
        #     reference = ["exome", "transcriptome", "genome"]
        # ),
        expand(
            bwa_ref + "report_{reference}.html",
            reference = ["exome", "transcriptome", "genome"]
        ),
        # # bwa_exons
        # expand(
        #     bwa_exfi + "{sample}_vs_{exons}.bam.bai",
        #     sample = dna_pe,
        #     exons = ["exons"],
        # ),
        ## dist
        dist + "exon_histogram.pdf",
        dist + "exon_density.pdf",
        ## pr
        pr + "pr.tsv"
