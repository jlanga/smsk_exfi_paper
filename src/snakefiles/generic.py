rule index_bam:
    input:
        bam = "{filename}.bam"
    output:
        bai = "{filename}.bam.bai"
    conda:
        "generic.yml"
    shell:
        "samtools index {input.bam}"


rule index_fa:
    input:
        fasta = "{filename}.fa"
    output:
        fai = "{filename}.fa.fai"
    conda:
        "generic.yml"
    shell:
        "samtools faidx {input.fasta}"



rule index_fasta:
    input:
        fasta = "{filename}.fasta"
    output:
        fai = "{filename}.fasta.fai"
    conda:
        "generic.yml"
    shell:
        "samtools faidx {input.fasta}"
