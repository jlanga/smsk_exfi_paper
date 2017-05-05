rule bwa_exons_index:
    input:
        fasta = exons + "{exons}.fa"
    output:
        mock = protected(
            touch(bwa_exons + "{exons}")
        ),
        other_files = protected(
            expand(
                bwa_exons + "{exons}.{extension}",
                extension = "ann bwt pac sa".split(),
                exons = "{exons}"
            )
        )
    log:
        bwa_exons + "index_{exons}.log"
    benchmark:
        bwa_exons + "index_{exons}.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log} 1>&2"



rule bwa_exons_align:
    input:
        forward = lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reverse = lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reference = bwa_exons + "{exons}",
    output:
        bam = bwa_exons+ "{sample}_vs_{exons}.bam"
    threads:
        THREAD_LIMIT
    log:
        bwa_exons+ "{sample}_vs_{exons}.log"
    benchmark:
        bwa_exons+ "{sample}_vs_{exons}.json"
    shell:
        "(bwa mem "
            "-t {threads} "
            "{input.reference} "
            "{input.forward} "
            "{input.reverse} "
        "| samtools view "
            "-Shu "
        "| samtools sort "
            "-l 9 "
            "--output-fmt BAM "
            "-@ {threads} "
        "> {output.bam}) "
        "2> {log}"



rule bwa_exons_stats:
    input:
        bam = bwa_exons + "{sample}_vs_{exons}.bam"
    output:
        stats= bwa_exons + "{sample}_vs_{exons}.stats"
    shell:
        "bamtools stats -in {input.bam} > {output.stats}"



rule bwa_exons_report:
    input:
        stats = expand(
            bwa_exons+ "{sample}_vs_{exons}.stats",
            exons = ["raw", "filtered_by_length", "filtered_by_extensibility"],
            sample = dna_pe
        )
    output:
        list_stats = temp(
            bwa_exons+ "list_{exons}.tsv"
        ),
        report = bwa_exons+ "report_{exons}.html"
    params:
        name = bwa_exons+ "report_{exons}"
    log:
        bwa_exons+ "report_{exons}.log"
    benchmark:  
        bwa_exons+ "report_{exons}.json"
    shell:
        "(ls -1 {input.stats} > {output.list_stats}; "
        "multiqc "
            "--filename {params.name} "
            "--file-list {output.list_stats}) "
        "2> {log}"