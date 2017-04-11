rule bwa_index:
    input:
        fasta = lambda wildcards: config["reference"][wildcards.reference]
    output:
        mock = touch(bwa + "{reference}"),
        other_files = expand(
            bwa + "{reference}.{extension}",
            extension = "ann bwt pac sa".split(),
            reference = "{reference}"
        )
    threads:
        1
    log:
        bwa + "index_{reference}.log"
    benchmark:
        bwa + "index_{reference}.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log}"





rule bwa_align:
    input:
        fasta = exons + "{exon_file}.fa",
        reference = bwa + "{reference}"
    output:
        bam = bwa + "{exon_file}_vs_{reference}.bam"
    threads:
        THREAD_LIMIT
    log:
        bwa + "{exon_file}_vs_{reference}.log"
    benchmark:
        bwa + "{exon_file}_vs_{reference}.json"
    shell:
        "(bwa mem "
            "-t {threads} "
            "{input.reference} "
            "{input.fasta} "
        "| samtools view "
            "-Shu "
        "| samtools sort "
            "-l 9 "
            "--output-fmt BAM "
            "-@ {threads} "
        "> {output.bam}) "
        "2> {log}"