# pylint: disable=syntax-error

import pandas as pd
import yaml

from snakemake.utils import min_version
min_version("5.3")

shell.prefix("set -euo pipefail;")

params = yaml.safe_load(open("params.yml", "r"))
features = yaml.safe_load(open("features.yml", "r"))
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
include: snakefiles + "chopstitch_review.smk"


# rule exfi_varying_pr_identity:
#     """Just check divergences in bedtools intersect minimum identity"""
#     expand(
#         PR + "drer.k25.l2.m60.100.baited.exfi.{identity}.tsv",
#         identity=[x / 10 for x in range(1, 10 + 1)]
#     )


rule exfi_baited_vs_unbaited:
    """In the zebrafish dataset, find out the effect of baiting vs using the 
    full dataset in the Bloom filters. Memory takes values in 4, 8, 12, ..., 64
    Gigabytes.

    For figure 3.
    """
    input:
        expand(
            PR + "drer.k25.l2.m{size}.100.{type}.exfi.0.95.tsv",
            size=range(4, 60 + 1, 4),
            type=["unbaited", "baited"],
        ),
        expand(
            BWA + "drer.k25.l2.m{size}.100.{type}.exfi.stats.tsv",
            size=range(4, 60 + 1, 4),
            type=["unbaited", "baited"],
        )

rule exfi_optimal_k:
    """In the zebrafish dataset, find optimal k value in terms of mapping and 
    bedtools metrics. Check values for low and high memory usage. k takes values
    in 21, 23, 25, ..., 65.

    For figure 4.
    """
    input:
        expand(
            PR + "drer.k{kmer}.l2.m{m}.100.baited.exfi.0.95.tsv",
            kmer=range(21, 65 + 1, 2),
            m=[4, 60]
        ),
        expand(
            BWA + "drer.k{kmer}.l2.m{m}.100.baited.exfi.stats.tsv",
            kmer=range(21, 65 + 1, 2),
            m=[4, 60]
        )


rule exfi_optimal_depth:
    """
    In the zebrafish dataset, find out the optimal experiment depth when k=25, 
    high and low memory, and 10% increments of depth in zebrafish (~6x each)

    For figure 5.
    """
    input:
        expand(
            PR + "drer.k25.l2.m{m}.{sampling}.baited.exfi.0.95.tsv",
            sampling=range(10,100 + 1, 10),
            m=[4, 60]
        ),
        expand(
            BWA + "drer.k25.l2.m{m}.{sampling}.baited.exfi.stats.tsv",
            sampling=range(10,100 + 1, 10),
            m=[4, 60]
        )


rule comparison_drer:
    """
    Compare how CS, EXFI and GMAP behave on zebrafish varying the memory/FPR
    values. Parameters:
    - EXFI: k=25, m=4G, l=2
    - CS_k25: k=25, FPR
    - GMAP: default

    Note: some values 
    For figure 6.
    """
    input:
        



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

        # Unknown why
        expand(
            PR + "drer.k25.fpr{fpr}.chopstitch.0.95.tsv",
            fpr=[0.01, 0.02, 0.03, 0.04, 0.05, "0.10", 0.15, "0.20"]
        ),
        expand(
            BWA + "drer.k25.fpr{fpr}.chopstitch.stats.tsv",
            fpr=[0.01, 0.02, 0.03, 0.04, 0.05, "0.10", 0.15, "0.20"]
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


rule comparison_species:
    """
    Comparison between the 6 species.
    Parameters:
        - EXFI:
            - k=25, l=2, m=4G in all species
            - k=25, l=2, m=60G in megagenomes, since 4 was too poor
        - CS:
            - k=25, FPR1=FPR2=0.01 in drer, hsap, char,  
            - k=25, FPR1=FPR2=FPR(exfi_4G) 
        - GMAP:
            - default in every case
            - gmapl in amex and plam because they are too big
    
    For table 2 (original)
    """
    input:
        # -- with reference genome --
        # drer
        PR + "drer.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "drer.k25.fpr0.01.chopstitch.0.95.tsv",
        PR + "drer.gmap.0.95.tsv",
        BWA + "drer.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "drer.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "drer.gmap.stats.tsv",
        # hsap
        PR + "hsap.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "hsap.k25.fpr0.15.chopstitch.0.95.tsv",
        PR + "hsap.gmap.0.95.tsv",
        BWA + "hsap.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "hsap.k25.fpr0.15.chopstitch.stats.tsv",
        BWA + "hsap.gmap.statstsv"
        # ssal
        PR + "ssal.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "ssal.k25.fpr0.01.chopstitch.0.95.tsv",
        PR + "ssal.gmap.0.95.tsv",
        BWA + "ssal.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "ssal.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + 'ssal.gmap.stats.tsv',
        # char_ref
        PR + "char_ref.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "char_ref.k25.fpr0.01.chopstitch.0.95.tsv",
        PR + "char_ref.gmap.0.95.tsv",
        BWA + "char_ref.k25.fpr0.01.chopstitch.stats.tsv",
        BWA + "char_ref.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + 'char_ref.gmap.stats.tsv',
        # -- no reference genome --
        # char
        BWA + 'char.gmap.stats.tsv',
        BWA + "char.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "char.k25.fpr0.01.chopstitch.stats.tsv",
        # amex
        # BWA + "amex.k25.fpr0.01.chopstitch.stats.tsv",  # std::mem_alloc()
        BWA + "amex.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "amex.k25.l2.m60.100.baited.exfi.stats.tsv",
        # BWA + 'amex.gmap.stats.tsv',  # Loads the index, maps, but after a while occupies all RAM
        # plam
        # BWA + "plam.k25.fpr0.01.chopstitch.stats.tsv",  # std::mem_alloc()
        BWA + "plam.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "plam.k25.l2.m60.100.baited.exfi.stats.tsv",
        BWA + 'plam.gmap.stats.tsv',


rule cs_varying_fprs:
    """
    Part of the first review. 
    Run chopstitch varying FPR1 and FPR2 in the six species. Commented cases
    don't work in 64GB server.
    Also, run chopstitch with default k=50, since authors says it is better to
    predict microexons and perform correction.
    """
    input:
        # -- with reference genome --
        # # Zebrafish case
        expand(  # Drer fpr1 big, fpr2 in 0.01 0.05
            PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv",
            sample=["drer"],
            kmer=50,
            fpr1="0.01 0.05 0.10 0.15 0.16 0.17 0.18 0.19".split(),
            fpr2="0.01 0.05".split(),
            identity="0.95"
        ),
        expand(  # Drer k=50 and multiple frp1 and fpr2
            BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
            sample=["drer"],
            kmer=50,
            fpr1="0.01 0.05 0.10 0.15".split(),
            fpr2="0.01 0.05".split(),
            identity="0.95"
        ),
        # # Human case
        # expand(
        #     PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv",
        #     sample=["hsap"],
        #     kmer=50,
        #     fpr1="0.01 0.05 0.10 0.15 0.16 0.17 0.18 0.19".split(),
        #     fpr2="0.01 0.05".split(),
        #     identity="0.95"
        # ),
        # expand(
        #     BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
        #     sample=["hsap"],
        #     kmer=50,
        #     fpr1="0.01 0.05 0.10 0.15".split(),
        #     fpr2="0.01 0.05",
        #     identity="0.95"
        # ),
        expand(  # Fishes 
            BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
            sample="char char_ref ssal".split(),
            kmer=50,
            fpr1=0.01,
            fpr2=0.01
        ),
        expand( # Fishes with reference
            PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv",
            sample="char_ref ssal".split(),
            kmer=50,
            fpr1=0.01,
            fpr2=0.01,
            identity=0.95
        ),
        # expand(  # Megagenomes - Doesn't work < 64 Gb
        #     BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
        #     sample="plam amex".split(),
        #     kmer=50,
        #     fpr1="0.01 0.05 0.10 0.15".split(),
        #     fpr2="0.01 0.05".split()
        # ),
        expand(  # Megagenomes
            BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
            sample="plam amex".split(),
            kmer=50,
            fpr1="0.20".split(),
            fpr2="0.01 0.05".split()
        ),
        expand(  # Human, high fpr1, low fpr2
            PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv",
            sample="hsap".split(),
            kmer=50,
            fpr1="0.16 0.18 0.20".split(),
            fpr2="0.01 0.05".split(),
            identity=0.95
        )



rule comparison_species_version2:
    input:
        # -- with reference genome --
        # drer
        PR + "drer.k25.fpr0.01.chopstitch.0.95.tsv",
        PR_R + "drer.k_50.fpr1_0.01.fpr2_0.01.chopstitch.0.95.tsv",
        PR + "drer.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "drer.gmap.0.95.tsv",
        BWA + "drer.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "drer.k25.fpr0.01.chopstitch.stats.tsv",
        BWA_R + "drer.k_25.fpr1_0.01.fpr2_0.01.chopstitch.stats.tsv",
        BWA + "drer.gmap.stats.tsv",
        # hsap
        PR + "hsap.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "hsap.k25.fpr0.15.chopstitch.0.95.tsv",
        # expand(  # None of these work
        #     PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv",
        #     sample=["hsap"],
        #     kmer=50,
        #     fpr1="0.01 0.05 0.10 0.15 0.16 0.17 0.18 0.19".split(),
        #     fpr2="0.01 0.05".split(),
        #     identity="0.95"
        # ),
        PR + "hsap.gmap.0.95.tsv",
        BWA + "hsap.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "hsap.k25.fpr0.15.chopstitch.stats.tsv",
        # expand(  # None of these work
        #     BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
        #     sample=["hsap"],
        #     kmer=50,
        #     fpr1="0.01 0.05 0.10 0.15".split(),
        #     fpr2="0.01 0.05",
        #     identity="0.95"
        # ),
        BWA + "hsap.gmap.statstsv"
        # ssal
        PR + "ssal.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "ssal.k25.fpr0.01.chopstitch.0.95.tsv",
        PR_R + "ssal.k_50.fpr1_0.01.fpr2_0.02.chopstitch.stats.tsv",
        PR + "ssal.gmap.0.95.tsv",
        BWA + "ssal.k25.fpr0.01.chopstitch.stats.tsv",
        BWA_R + "ssal.k_50.fpr1_0.01.fpr2_0.01.chopstitch.stats.tsv",
        BWA + "ssal.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + 'ssal.gmap.stats.tsv',
        # char_ref
        PR + "char_ref.k25.l2.m4.100.baited.exfi.0.95.tsv",
        PR + "char_ref.k25.fpr0.01.chopstitch.0.95.tsv",
        PR_R + "char_ref.k_25.fpr1_0.01.fpr2_0.01.chopstitch.0.95.tsv",
        PR + "char_ref.gmap.0.95.tsv",
        BWA + "char_ref.k25.fpr0.01.chopstitch.stats.tsv",
        BWA_R + "char_ref.k_50.fpr1_0.01.fpr2_0.01.chopstitch.stats.tsv",
        BWA + "char_ref.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + 'char_ref.gmap.stats.tsv',
        # -- no reference genome --
        # char
        BWA + 'char.gmap.stats.tsv',
        BWA + "char.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA_R + "char.k_50.fpr1_0.01.fpr2_0.01.chopstitch.stats.tsv",
        BWA + "char.k25.fpr0.01.chopstitch.stats.tsv",
        # amex
        # BWA + "amex.k25.fpr0.01.chopstitch.stats.tsv",  # std::mem_alloc()
        BWA + "amex.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "amex.k25.l2.m60.100.baited.exfi.stats.tsv",
        # BWA + 'amex.gmap.stats.tsv',  # Loads the index, maps, but after a while clogs the RAM
        # plam
        # BWA + "plam.k25.fpr0.01.chopstitch.stats.tsv",  # std::mem_alloc()
        expand(  # Megagenomes
            BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv",
            sample="plam amex".split(),
            kmer=50,
            fpr1="0.20".split(),
            fpr2="0.01 0.05".split()
        )
        BWA + "plam.k25.l2.m4.100.baited.exfi.stats.tsv",
        BWA + "plam.k25.l2.m60.100.baited.exfi.stats.tsv",
        BWA + 'plam.gmap.stats.tsv',         

rule all:
    input:
        rules.exfi_baited_vs_unbaited.input,
        rules.exfi_optimal_k.input,
        rules.exfi_optimal_depth.input,
        rules.comparison_species_version2.input