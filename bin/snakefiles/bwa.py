rule bwa_index_exome:
    input:
        fasta = config["reference"]["exome"]
    output:
        mock = touch(bwa + "exome"),
        other_files = expand(
            bwa + "exome.{extension}",
            extension = "ann bwt pac sa".split()
        )
    threads:
        1
    log:
        bwa + "index_exome.log"
    benchmark:
        bwa + "index_exome.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log}"



rule bwa_index_transcriptome:
    input:
        fasta = config["reference"]["transcriptome"]
    output:
        mock = touch(bwa + "transcriptome"),
        other_files = expand(
            bwa + "transcriptome.{extension}",
            extension = "amb ann bwt pac sa".split()
        )
    threads:
        1
    log:
        bwa + "index_transcriptome.log"
    benchmark:
        bwa + "index_transcriptome.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log}"



rule bwa_index_genome:
    input:
        fasta = config["reference"]["genome"]
    output:
        mock = touch(bwa + "genome"),
        other_files = expand(
            bwa + "genome.{extension}",
            extension = "ann bwt pac sa".split()
        )
    threads:
        1
    log:
        bwa + "index_genome.log"
    benchmark:
        bwa + "index_genome.json"
    shell:
        "bwa index "
            "-p {output.mock} "
            "{input.fasta} "
        "2> {log}"



rule align:
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