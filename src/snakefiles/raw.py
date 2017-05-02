rule raw_link_pe_sample:
    input:
        forward= lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reverse= lambda wildcards: config["samples"][wildcards.sample]["reverse"]
    output:
        forward= raw + "{sample}_1.fq.gz",
        reverse= raw + "{sample}_2.fq.gz"
    log:
        raw + "link_dna_pe_{sample}.log"
    benchmark:
        raw + "link_dna_pe_{sample}.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.forward}) "
            "{output.forward} 2> {log}; "
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.reverse}) "
            "{output.reverse} 2>> {log}"



rule raw_link_assembly:
    input:
        fasta= config["assembly"]
    output:
        fasta= raw + "assembly.fa"
    log:
        raw + "link_assembly.log"
    benchmark:
        raw + "link_assembly.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"



rule raw_link_exome:
    input:
        fasta= config["reference"]["exome"]
    output:
        fasta= raw + "exome.fa"
    log:
        raw + "link_exome.log"
    benchmark:
        raw + "link_exome.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"



rule raw_link_transcriptome:
    input:
        fasta= config["reference"]["transcriptome"]
    output:
        fasta= raw + "transcriptome.fa"
    log:
        raw + "link_transcriptome.log"
    benchmark:
        raw + "link_transcriptome.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"



rule raw_link_genome:
    input:
        fasta= config["reference"]["genome"]
    output:
        fasta= raw + "genome.fa"
    log:
        raw + "link_genome.log"
    benchmark:
        raw + "link_genome.json"
    shell:
        "ln "
            "--symbolic "
            "$(readlink --canonicalize {input.fasta}) "
            "{output.fasta} 2> {log}"