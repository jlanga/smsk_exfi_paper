# pylint: disable=syntax-error

import pandas as pd
import yaml

from snakemake.utils import min_version
min_version("5.3")

shell.prefix("set -euo pipefail;")

params = yaml.load(open("params.yml", "r"))
features = yaml.load(open("features.yml", "r"))
samples = pd.read_csv("samples.tsv", sep='\t')

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
include: snakefiles + "snp.py"
include: snakefiles + "bowtie2.py"


rule all:
    input:
        # chopstitch fpr1 and same as exfi
        expand(
            BWA + "char.k25.fpr{fpr}.chopstitch.stats.tsv",
            fpr=[0.01, 0.0017]
        ),
        expand(
            PR + "drer.k25.fpr{fpr}.chopstitch.0.95.tsv",
            fpr=[0.01, 0.0042],
        ),
        expand(  # Doesn't work below 0.15
            PR + "hsap.k25.fpr{fpr}.chopstitch.0.95.tsv",
            fpr=[0.15, 0.2, 0.3, 0.4],
        ),
        expand(
            BWA + "hsap.k25.fpr{fpr}.chopstitch.stats.tsv",
            fpr=[0.15, 0.2, 0.3, 0.4]
        ),
        # exfi_baited vs unbaited, varying mem and fpr
        # drer memory exfi baited and unbaited
        expand(
            PR + "drer.k25.l2.m{size}.100.{type}.exfi.0.95.tsv",
            size=range(4, 60 + 1, 4),
            type=["unbaited", "baited"],
        ),
        expand(
            BWA + "drer.k25.l2.m{size}.100.{type}.exfi.stats.tsv",
            size=range(4, 60 + 1, 4),
            type=["unbaited", "baited"],
        ),
        expand(
            PR + "drer.k25.fpr{fpr}.chopstitch.0.95.tsv",
            fpr=[0.01, 0.02, 0.03, 0.04, 0.05, "0.10", 0.15, "0.20"]
        ),
        expand(
            BWA + "drer.k25.fpr{fpr}.chopstitch.stats.tsv",
            fpr=[0.01, 0.02, 0.03, 0.04, 0.05, "0.10", 0.15, "0.20"]
        ),
        # # drer varying k
        expand(
            PR + "drer.k{kmer}.l2.m60.100.baited.exfi.0.95.tsv",
            kmer=range(21, 65 + 1, 2),
        ),
        expand(
            BWA + "drer.k{kmer}.l2.m60.100.baited.exfi.stats.tsv",
            kmer=range(21, 65 + 1, 2),
        ),
        # # drer varying sampling
        expand(
            PR + "drer.k25.l2.m60.{sampling}.baited.exfi.0.95.tsv",
            sampling=range(10,100 + 1, 10),
        ),
        expand(
            BWA + "drer.k25.l2.m60.{sampling}.baited.exfi.stats.tsv",
            sampling=range(10,100 + 1, 10),
        ),
        # # drer varying identity
        expand(
            PR + "drer.k25.l2.m60.100.baited.exfi.{identity}.tsv",
            identity=[x / 10 for x in range(1, 10 + 1)]
        ),
        # char mappings: exfi_baited, exfi_unbaited chopstitch with same fpr as
        # unbaited
        expand(
            BWA + "char.k25.l2.m60.100.{type}.exfi.stats.tsv",
            type=["baited", "unbaited"]
        ),
        BWA + "char.k25.fpr0.0017.chopstitch.stats.tsv",
        # # sample cases for exfi m60 and exfi frugal (4G)
        expand(
            PR + "{sample}.k25.l2.m{size}.100.baited.exfi.0.95.tsv",
            sample=["drer", "hsap"],
            size=[4, 60],
        ),
        expand(
            BWA + "{sample}.k25.l2.m{size}.100.baited.exfi.stats.tsv",
            sample=["drer", "hsap"],
            size=[4, 60],
        ),
        expand(
            PR + "{sample}.k25.l2.m{size}.100.baited.exfi.0.95.tsv",
            sample=["drer", "hsap"],
            size=[4, 60],
        ),
        expand(
            BWA + "{sample}.k25.l2.m4.100.baited.exfi.stats.tsv",
            sample=["char", "drer", "hsap", "ssal"],
        ),
        expand(
            PR + "{sample}.gmap.{identity}.tsv",
            sample=["drer", "hsap"],
            identity=0.95
        ),
        expand(
            BWA + '{sample}.gmap.stats.tsv',
            sample = ['drer', 'hsap']
        ),
        PR + "ssal.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "ssal.k25.fpr0.01.chopstitch.0.95.tsv",
        PR + "ssal.gmap.0.95.tsv",
        BWA + "drer.k25.fpr0.0042.chopstitch.stats.tsv",
        BWA + "drer.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "hsap.k25.fpr0.15.chopstitch.stats.tsv",
        # EXFI + "ttin.k25.l2.m4.100.baited.gfa",
        BWA + 'char.gmap.stats.tsv',
        BWA + "char.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "char.k25.fpr0.01.chopstitch.stats.tsv",
        # BWA + 'haxy.gmap.stats.tsv',
        # BWA + "haxy.k25.l2.m4.100.baited.exfi.stats.tsv",
        # BWA + "haxy.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "ssal.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "ssal.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + 'ssal.gmap.stats.tsv',
        #BWA + "amex.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "amex.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "amex.k25.l2.m60.100.baited.exfi.stats.tsv",
        #BWA + 'amex.gmap.stats.tsv',
        #BWA + "plam.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "plam.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "plam.k25.l2.m60.100.baited.exfi.stats.tsv",
        BWA + 'plam.gmap.stats.tsv',
        PR + "char_ref.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "char_ref.k25.fpr0.01.chopstitch.0.95.tsv",
        PR + "char_ref.gmap.0.95.tsv",
        BWA + "char_ref.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "char_ref.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + 'char_ref.gmap.stats.tsv',
