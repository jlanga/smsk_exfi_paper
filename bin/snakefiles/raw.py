rule raw_link_pe_sample:
    output:
        forward= raw + "{sample}_1.fq.gz",
        reverse= raw + "{sample}_2.fq.gz"
    params:
        forward= lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reverse= lambda wildcards: config["samples"][wildcards.sample]["reverse"]
    log:
        raw + "link_dna_pe_{sample}.log"
    benchmark:
        raw + "link_dna_pe_{sample}.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {params.forward}) "
            "{output.forward} 2> {log}; "
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {params.reverse}) "
            "{output.reverse} 2>> {log}"



rule raw_link_assembly:
    output:
        fasta= raw + "transcriptome.fa"
    params:
        fasta= config["assembly"]
    log:
        raw + "link_assembly.log"
    benchmark:
        raw + "link_assembly.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {params.fasta}) "
            "{output.fasta} 2> {log}"



# rule raw_index_assembly:
#     input:
#         fasta= raw + "transcriptome.fa"
#     output:
#         fai= raw + "transcriptome.fa.fai"
#     log:
#         raw + "index_assembly.log"
#     benchmark:
#         raw + "index_assembly.json"
#     shell:
#         "samtools faidx {input.fasta} 2> {log}"


rule raw_link_exome:
    input:
        fasta = config["reference"]["exome"]
    output:
        fasta = raw + "exome.fa"
    shell: 
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta}"