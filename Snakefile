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
include: snakefiles + "pr.py"
include: snakefiles + "bwa.py"



rule all:
    input:
        ## raw
        # expand(
        #    raw + "{sample}_{end}.fq.gz",
        #    sample = dna_pe,
        #    end = "1 2".split(" ")
        # ),
        # raw + "transcriptome.fa.fai",
        ## exfi
        # expand(
        #    exfi + "k{kmer}_l{levels}_m{size}.bloom",
        #    kmer = config["exfi"]["kmer"],
        #    levels = config["exfi"]["levels"],
        #    size = config["exfi"]["size"]
        # ),
        # exfi + "splice_graph.gfa",
        # exfi + "exons.fa",
        # exfi + "gapped_transcripts.fa",
        ## pr
        pr + "pr.tsv",
