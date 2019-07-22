def get_reads(wildcards):
    sample = wildcards.sample
    library = wildcards.library
    forward, reverse = (
        samples
        [(samples["sample"] == sample and samples["library"] == library)]
        [["forward", "reverse"]]
        .values
        .tolist()[0]
    )
    return forward, reverse


rule raw_link_pe_sample:
    input:
        get_reads
    output:
        forward = RAW + "{sample}_{library}_1.fq.gz",
        reverse = RAW + "{sample}_{library}_2.fq.gz"
    log:
        RAW + "link_dna_pe_{sample}_{library}.log"
    benchmark:
        RAW + "link_dna_pe_{sample}_{library}.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input[0]}) "
            "{output.forward} 2> {log}; "
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input[1]}) "
            "{output.reverse} 2>> {log}"


def get_transcriptome(wildcards):
    return features[wildcards.sample]["transcriptome"]


def get_genome(wildcards):
    return features[wildcards.sample]["genome"]


def get_annotation(wildcards):
    return features[wildcards.sample]["annotation"]


rule raw_link_transcriptome:
    input:
        get_transcriptome
    output:
        RAW + "{sample}.rna.fa"
    shell:
        "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule raw_link_genome:
    input:
        get_genome
    output:
        RAW + "{sample}.dna.fa"
    shell:
        "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule raw_link_annotation:
    input:
        get_annotation
    output:
        RAW + "{sample}.gff3"
    shell:
        "ln --symbolic $(readlink --canonicalize {input}) {output}"


rule raw_reference:
    input:
        expand(
            RAW + "{sample}.{ending}",
            sample=SPECIES,
            ending=["dna.fa", "rna.fa", "gff3"]
        )
