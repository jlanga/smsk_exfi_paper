rule bwa_index:
    input:
        reference = RAW + "{sample}.dna.fa"
    output:
        expand(
            BWA + "{sample}.{suffix}",
            sample="{sample}",
            suffix = "amb ann bwt pac sa".split()
        ),
        mock = touch(BWA + "{sample}.bwa")
    log:
        BWA + "{sample}_index.log"
    benchmark:
        BWA + "{sample}_index.bmk"
    params:
        prefix = BWA + "{sample}"
    conda: "bwa.yml"
    shell:
        "bwa index "
            "-p {params.prefix} "
            "{input.reference} "
        "2> {log} 1>&2"



rule bwa_map_exfi:
    input:
        index = BWA + "{sample}.bwa",
        exons = EXFI + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.fa"
    output:
        bam = BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam"
    log:
        BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam.log"
    benchmark:
        BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam.bmk"
    conda: "bwa.yml"
    threads: 4
    params:
        index = BWA + "{sample}"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "-D 0.1 "
            "{params.index} "
            "{input.exons} "
        "| samtools sort "
            "-l 9 "
            "-@ {threads} "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"



rule bwa_map_chopstitch:
    input:
        exons = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.exons.fa",
        index = BWA + "{sample}.bwa"
    output:
        bam = BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bam"
    log:
        BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bam.log"
    benchmark:
         BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bam.bmk"
    conda: "bwa.yml"
    threads: 4
    params:
        index = BWA + "{sample}"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "-D 0.1 "
            "{params.index} "
            "{input.exons} "
        "| samtools sort "
            "-@ {threads} "
            "-l 9 "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"


rule bwa_map_gmap:
    input:
        exons = GMAP + "{sample}.exons.fa",
        index = BWA + "{sample}.bwa"
    output:
        bam = BWA + "{sample}.gmap.bam"
    log:
        BWA + "{sample}.gmap.bam.log"
    benchmark:
         BWA + "{sample}.gmap.bam.bmk"
    conda: "bwa.yml"
    threads: 4
    params:
        index = BWA + "{sample}"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "-D 0.1 "
            "{params.index} "
            "{input.exons} "
        "| samtools sort "
            "-l 9 "
            "-@ {threads} "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"



rule bwa_stats_exfi:
    input:
        bam = BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam",
        fa = RAW + "{sample}.rna.fa"
    output: BWA + \
        "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.stats.tsv"
    log: BWA + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.stats.log"
    benchmark: BWA + \
        "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.stats.bmk"
    conda: "bwa.yml"
    shell:  'bash src/compute_bwa_stats.sh {input} > {output} 2> {log}'


rule bwa_stats_chopstitch:
    input:
        bam = BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bam",
        fa = RAW + "{sample}.rna.fa"
    output: BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.stats.tsv"
    log: BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.stats.log"
    benchmark: BWA + "{sample}.k{kmer}.fpr{fpr}.chopstitch.stats.bmk"
    conda: 'bwa.yml'
    shell: 'bash src/compute_bwa_stats.sh {input} > {output} 2> {log}'



rule bwa_stats_gmap:
    input:
        bam = BWA + "{sample}.gmap.bam",
        fa = RAW + "{sample}.rna.fa"
    output: BWA + "{sample}.gmap.stats.tsv"
    log: BWA + '{sample}.gmap.stats.tsv.log'
    benchmark: BWA + '{sample}.gmap.statsself.tsv.bmk'
    conda: 'bwa.yml'
    shell: 'bash src/compute_bwa_stats.sh {input} > {output} 2> {log}'
