def get_reads(wildcards):
    sample = wildcards.sample
    forward, reverse = (
        samples
        [(samples["sample"] == sample)]
        [["forward", "reverse"]]
        .values
        .tolist()[0]
    )
    return forward, reverse


rule raw_link_pe_sample:
    input:
        get_reads
    output:
        forward = RAW + "{sample}_1.fq.gz",
        reverse = RAW + "{sample}_2.fq.gz"
    log:
        RAW + "link_dna_pe_{sample}.log"
    benchmark:
        RAW + "link_dna_pe_{sample}.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input[0]}) "
            "{output.forward} 2> {log}; "
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input[0]}) "
            "{output.reverse} 2>> {log}"



rule raw_link_assembly:
    input:
        fasta = features["assembly"]
    output:
        fasta = RAW + "assembly.fa"
    log:
        RAW + "link_assembly.log"
    benchmark:
        RAW + "link_assembly.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"



# rule raw_link_exome:
#     input:
#         fasta= config["reference"]["exome"]
#     output:
#         fasta= RAW + "exome.fa"
#     log:
#         RAW + "link_exome.log"
#     benchmark:
#         RAW + "link_exome.json"
#     shell:
#         "ln "
#             "--symbolic "
#             "$(readlink --canonicalize {input.fasta}) "
#             "{output.fasta} 2> {log}"



# rule raw_reduce_exome:
#     input:
#         fasta = RAW + "exome.fa"
#     output:
#         fasta = RAW + "exome_reduced.fa"
#     shell:
#         "reduce_exons "
#             "--input-fasta {input.fasta} "
#             "--output-fasta {output.fasta}"



rule raw_link_transcriptome:
    input:
        fasta= features["reference"]["transcriptome"]
    output:
        fasta= RAW + "transcriptome.fa"
    log:
        RAW + "link_transcriptome.log"
    benchmark:
        RAW + "link_transcriptome.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"



rule raw_link_genome:
    input:
        fasta= features["reference"]["genome"]
    output:
        fasta= RAW + "genome.fa"
    log:
        RAW + "link_genome.log"
    benchmark:
        RAW + "link_genome.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"


rule raw_link_annotation:
    input:
        gff3_gz = features["reference"]["annotation"]
    output:
        gff3_gz = RAW + "annotation.gff3.gz"
    log:
        RAW + "link_annotation.log"
    benchmark:
        RAW + "link_annotation.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.gff3_gz}) "
            "{output.gff3_gz} "
        "2> {log}"
