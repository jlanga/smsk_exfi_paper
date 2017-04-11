rule bwa_index_exome:
    input:
        fasta = config["reference"]["exons"]
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



rule bwa_raw_vs_exome:
    input:
        fasta = exons + "raw.fa",
        reference = bwa + "exome"
    output:
        bam = bwa + "raw_vs_exome.bam"
    threads:
        24
    log:
        bwa + "raw_vs_exome.log"
    benchmark:
        bwa + "raw_vs_exome.json"
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



rule bwa_raw_vs_transcriptome:
    input:
        fasta = exons + "raw.fa",
        reference = bwa + "transcriptome"
    output:
        bam = bwa + "raw_vs_transcriptome.bam"
    threads:
        24
    log:
        bwa + "raw_vs_transcriptome.log"
    benchmark:
        bwa + "raw_vs_transcriptome.json"
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



rule bwa_raw_vs_genome:
    input:
        fasta = exons + "raw.fa",
        reference = bwa + "genome"
    output:
        bam = bwa + "raw_vs_genome.bam"
    threads:
        24
    log:
        bwa + "raw_vs_genome.log"
    benchmark:
        bwa + "raw_vs_genome.json"
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



rule bwa_filtered_vs_exome:
    input:
        fasta = exons + "filtered.fa",
        reference = bwa + "exome"
    output:
        bam = bwa + "filtered_vs_exome.bam"
    threads:
        24
    log:
        bwa + "filtered_vs_exome.log"
    benchmark:
        bwa + "filtered_vs_exome.json"
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



rule bwa_filtered_vs_transcriptome:
    input:
        fasta = exons + "filtered.fa",
        reference = bwa + "transcriptome"
    output:
        bam = bwa + "filtered_vs_transcriptome.bam"
    threads:
        24
    log:
        bwa + "filtered_vs_transcriptome.log"
    benchmark:
        bwa + "filtered_vs_transcriptome.json"
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



rule bwa_filtered_vs_genome:
    input:
        fasta = exons + "filtered.fa",
        reference = bwa + "genome"
    output:
        bam = bwa + "filtered_vs_genome.bam"
    threads:
        24
    log:
        bwa + "filtered_vs_genome.log"
    benchmark:
        bwa + "filtered_vs_genome.json"
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