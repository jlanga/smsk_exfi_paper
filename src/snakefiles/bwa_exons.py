rule bwa_exfi_index:
    input:
        fasta = exfi + "{exons}.fa"
    output:
        mock = protected(
            touch(bwa_exfi + "{exons}")
        ),
        other_files = protected(
            expand(
                bwa_exfi + "{exons}.{extension}",
                extension = "ann bwt pac sa".split(),
                exons = "{exons}"
            )
        )
    log:
        bwa_exfi + "index_{exons}.log"
    benchmark:
        bwa_exfi + "index_{exons}.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log} 1>&2"



rule bwa_exfi_align:
    input:
        forward = lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reverse = lambda wildcards: config["samples"][wildcards.sample]["forward"],
        reference = bwa_exfi + "{exons}",
    output:
        bam = bwa_exfi+ "{sample}_vs_{exons}.bam"
    threads:
        THREAD_LIMIT
    log:
        bwa_exfi+ "{sample}_vs_{exons}.log"
    benchmark:
        bwa_exfi+ "{sample}_vs_{exons}.json"
    shell:
        "(bwa mem "
            "-t {threads} "
            "{input.reference} "
            "{input.forward} "
            "{input.reverse} "
        "| samtools view "
            "-F4 -h "
        "| samtools sort "
            "-l 9 "
            "--output-fmt BAM "
            "-@ {threads} "
        "> {output.bam}) "
        "2> {log}"



rule bwa_exfi_stats:
    input:
        bam = bwa_exfi + "{sample}_vs_{exons}.bam"
    output:
        stats= bwa_exfi + "{sample}_vs_{exons}.stats"
    shell:
        "bamtools stats -in {input.bam} > {output.stats}"



rule bwa_exfi_report:
    input:
        stats = expand(
            bwa_exfi+ "{sample}_vs_{exons}.stats",
            exons = ["exons"],
            sample = dna_pe
        )
    output:
        list_stats = temp(
            bwa_exfi+ "list_{exons}.tsv"
        ),
        report = bwa_exfi+ "report_{exons}.html"
    params:
        name = bwa_exfi+ "report_{exons}"
    log:
        bwa_exfi+ "report_{exons}.log"
    benchmark:
        bwa_exfi+ "report_{exons}.json"
    shell:
        "(ls -1 {input.stats} > {output.list_stats}; "
        "multiqc "
            "--filename {params.name} "
            "--file-list {output.list_stats}) "
        "2> {log}"
