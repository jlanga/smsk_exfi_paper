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
SPECIES = samples["sample"].unique().tolist()


# Read subsnakefiles
snakefiles = "src/snakefiles/"
include: snakefiles + "folders.py"
include: snakefiles + "generic.py"
include: snakefiles + "clean.py"
include: snakefiles + "raw.py"
include: snakefiles + "exfi.py"
include: snakefiles + "chopstitch.py"
include: snakefiles + "pr.py"
include: snakefiles + "bwa.py"
include: snakefiles + "gmap.py"
include: snakefiles + "tuxedo2.py"



rule all:
    input:
        # chopstitch fpr1 and same as exfi
        expand(
            BWA + "{sample}.k{kmer}.fpr{fpr}.stats.tsv",
            sample="char",
            kmer=25,
            fpr=[0.01, 0.0017]
        ),
        expand(
            PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv",
            sample="drer",
            kmer=25,
            fpr=[0.01, 0.0042],
            identity=["0.90", "0.95", "1.00"]
        ),
        # expand(  # Doesn't work
        #     PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv",
        #     sample="hsap",
        #     kmer=25,
        #     fpr=[0.01, 0.0068],
        #     identity=0.90
        # ),
        # exfi_baited vs unbaited vs chopstitch, varying mem and fpr
        # drer memory exfi baited and unbaited
        expand(
            PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.{identity}.tsv",
            kmer=25,
            sample="drer",
            levels=2,
            size=range(4, 60 + 1, 4),
            sampling=100,
            type=["unbaited", "baited"],
            identity=["0.90", "0.95", "1.00"]
        ),
        expand(
            PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv",
            sample="drer",
            kmer=25,
            fpr=[0.01, 0.02, 0.03, 0.04, 0.05, "0.10", 0.15, "0.20"],
            identity=["0.90", "0.95", "1.00"]
        ),
        # # drer varying k
        expand(
            PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.{identity}.tsv",
            sample="drer",
            kmer=range(21, 35 + 1, 2),
            levels=2,
            size=60,
            sampling=100,
            type="baited",
            identity=["0.90", "0.95", "1.00"]
        ),
        # # drer varying sampling
        expand(
            PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.{identity}.tsv",
            sample="drer",
            kmer=25,
            levels=2,
            size=60,
            sampling=range(10,100 + 1, 10),
            type="baited",
            identity=["0.90", "0.95", "1.00"]
        ),
        # # drer varying identity
        expand(
            PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.{identity}.tsv",
            sample="drer",
            kmer=25,
            levels=2,
            size=60,
            sampling=100,
            type="baited",
            identity=[x / 10 for x in range(1, 10 + 1)]
        ),
        # char mappings: exfi_baited, exfi_unbaited chopstitch with same fpr as
        # unbaited
        expand(
             BWA + \
                "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.stats.tsv",
            sample="char",
            kmer=25,
            levels=2,
            size=60,
            sampling=100,
            type=["baited", "unbaited"]
        ),
        expand(
            BWA + "{sample}.k{kmer}.fpr{fpr}.stats.tsv",
            sample="char",
            kmer=25,
            fpr=0.0017
        ),
        # # sample cases for exfi m60 and exfi frugal (4G)
        expand(
            PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.{identity}.tsv",
            sample=["drer", "hsap"],
            kmer=25,
            levels=2,
            size=[4, 60],
            sampling=100,
            type="baited",
            identity=["0.90", "0.95", "1.00"]
        ),
        expand(
             BWA + \
                "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.stats.tsv",
            sample="char",
            kmer=25,
            levels=2,
            size=4,
            sampling=100,
            type=["baited"]
        ),
        expand(
            PR + "{sample}.gmap.{identity}.tsv",
            sample=["drer", "hsap"],
            identity=["0.90", "0.95", "1.00"]
        )
