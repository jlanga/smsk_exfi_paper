# pylint: disable=syntax-error

import pandas as pd
import yaml

from snakemake.utils import min_version
min_version("5.3")

shell.prefix("set -euo pipefail;")

params = yaml.load(open("params.yml", "r"))
features = yaml.load(open("features.yml", "r"))
samples = pd.read_table("samples.tsv")

singularity: "docker://continuumio/miniconda3:4.4.10"

# Some variables
SAMPLES =  samples["sample"].values.tolist()

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
        #    RAW + "{sample}_{end}.fq.gz",
        #    sample = dna_pe,
        #    end = "1 2".split(" ")
        # ),
        # RAW + "transcriptome.fa.fai",
        ## exfi
        # expand(
        #    EXFI + "k{kmer}_l{levels}_m{size}.bloom",
        #    kmer = config["exfi"]["kmer"],
        #    levels = config["exfi"]["levels"],
        #    size = config["exfi"]["size"]
        # ),
        # EXFI + "splice_graph.gfa",
        # EXFI + "exons.fa",
        # EXFI + "gapped_transcripts.fa",
        ## pr
        PR + "pr.tsv",
        BWA + "stats.tsv"
