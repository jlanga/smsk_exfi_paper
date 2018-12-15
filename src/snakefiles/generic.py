rule index_bam:
    input: "{filename}"
    output: "{filename}.bai"
    conda: "generic.yml"
    shell: "samtools index {input}"


rule index_fasta:
    input: "{filename}"
    output: "{filename}.fai"
    conda: "generic.yml"
    shell: "samtools faidx {input}"
