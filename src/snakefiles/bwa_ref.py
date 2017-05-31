rule bwa_ref_index:
    input:
        fasta = raw + "{reference}.fa"
    output:
        mock = protected(
            touch(bwa_ref+ "{reference}")
        ),
        other_files = protected(
            expand(
                bwa_ref+ "{reference}.{extension}",
                extension = "ann bwt pac sa".split(),
                reference = "{reference}"
            )
        )
    log:
        bwa_ref+ "index_{reference}.log"
    benchmark:
        bwa_ref+ "index_{reference}.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log} 1>&2"



rule bwa_ref_align:
    input:
        fasta = exons + "{exon_file}.fa",
        reference = bwa_ref+ "{reference}"
    output:
        bam = bwa_ref+ "{exon_file}_vs_{reference}.bam"
    threads:
        THREAD_LIMIT
    log:
        bwa_ref+ "{exon_file}_vs_{reference}.log"
    benchmark:
        bwa_ref+ "{exon_file}_vs_{reference}.json"
    shell:
        "(bwa bwasw "
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



rule bwa_ref_stats:
    input:
        bam = bwa_ref + "{exon_file}_vs_{reference}.bam"
    output:
        stats= bwa_ref + "{exon_file}_vs_{reference}.stats"
    shell:
        "bamtools stats -in {input.bam} > {output.stats}"



rule bwa_ref_report:
    input:
        stats = expand(
            bwa_ref+ "{exon_file}_vs_{reference}.stats",
            exon_file = ["raw", "filtered_by_length", "filtered_by_extensibility"],
            reference = "{reference}"
        )
    output:
        list_stats = temp(
            bwa_ref+ "list_{reference}.tsv"
        ),
        report = bwa_ref+ "report_{reference}.html"
    params:
        name = bwa_ref+ "report_{reference}"
    log:
        bwa_ref+ "report_{reference}.log"
    benchmark:
        bwa_ref+ "report_{reference}.json"
    shell:
        "(ls -1 {input.stats} > {output.list_stats}; "
        "multiqc "
            "--filename {params.name} "
            "--file-list {output.list_stats}) "
        "2> {log}"
