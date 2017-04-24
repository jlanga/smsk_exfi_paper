rule bwa_index:
    input:
        fasta = lambda wildcards: config["reference"][wildcards.reference]
    output:
        mock = protected(
            touch(bwa + "{reference}")
        ),
        other_files = protected(
            expand(
                bwa + "{reference}.{extension}",
                extension = "ann bwt pac sa".split(),
                reference = "{reference}"
            )
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
        "2> {log} 1>&2"



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



rule bwa_stats:
    input:
        bam = bwa + "{exon_file}_vs_{reference}.bam"
    output:
        stats= bwa + "{exon_file}_vs_{reference}.stats"
    shell:
        "samtools stats {input.bam} > {output.stats}"



rule bwa_report:
    input:
        stats = expand(
            bwa + "{exon_file}_vs_{reference}.stats",
            exon_file = ["raw", "filtered_by_length", "filtered_by_extensibility"],
            reference = "{reference}"
        )
    output:
        list_stats = temp(
            bwa + "list_{reference}.tsv"
        ),
        report = bwa + "report_{reference}.html"
    params:
        name = bwa + "report_{reference}"
    log:
        bwa + "report_{reference}.log"
    benchmark:  
        bwa + "report_{reference}.json"
    shell:
        "(ls -1 {input.stats} > {output.list_stats}; "
        "multiqc "
            "--filename {params.name} "
            "--file-list {output.list_stats}) "
        "2> {log}"