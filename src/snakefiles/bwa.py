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
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.bam"
    log:
        BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.log"
    benchmark:
        BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.bmk"
    conda: "bwa.yml"
    params:
        index = BWA + "{sample}"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "-D 0.1 "
            "{params.index} "
            "{input.exons} "
        "| samtools view -Shub "
        "| samtools sort "
            "-l 9 "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"



rule bwa_map_chopstitch:
    input:
        exons = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.processedexons.fa",
        index = BWA + "{sample}.bwa"
    output:
        bam = BWA + "{sample}.k{kmer}.fpr{fpr}.processedexons.bam"
    log:
        BWA + "{sample}.k{kmer}.fpr{fpr}.processedexons.log"
    benchmark:
         BWA + "{sample}.k{kmer}.fpr{fpr}.processedexons.bmk"
    conda: "bwa.yml"
    params:
        index = BWA + "{sample}"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "-D 0.1 "
            "{params.index} "
            "{input.exons} "
        "| samtools view -Shub "
        "| samtools sort "
            "-l 9 "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"



rule bwa_stats_exfi:
    input:
        bam = BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.bam"
    output:
        BWA + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.stats.tsv"
    log:
        BWA + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.stats.log"
    benchmark:
        BWA + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.stats.bmk"
    conda: "bwa.yml"
    shell:
        """
        n=$(samtools view {input.bam} | cut -f 1 | sort -u | wc -l)
        unmapped=$(samtools view {input.bam} | awk '$3 == "*"' | wc -l)
        mapped=$(samtools view {input.bam} | awk '$3 != "*"'| cut -f 1 | sort -u | wc -l)
        perfect=$(samtools view {input.bam} | awk '$6 ~ /^[0-9]+M$/' | cut -f 1 | sort -u | wc -l)
        unique=$(samtools view {input.bam} | cut -f 1 | sort | uniq -c | awk '$1 == 1' | wc -l)
        multi=$(samtools view {input.bam} | cut -f 1 | sort | uniq -c | awk '$1 > 1' | wc -l)

        echo exons"\t"unmapped"\t"mapped"\t"perfect_match"\t"uniquely_mapped"\t"multimapped > {output}
        echo $n"\t"$unmapped"\t"$mapped"\t"$perfect"\t"$unique"\t"$multi >> {output}
        """



rule bwa_stats_chopstitch:
    input:
        bam = BWA + "{sample}.k{kmer}.fpr{fpr}.processedexons.bam"
    output:
        tsv = BWA + "{sample}.k{kmer}.fpr{fpr}.stats.tsv"
    log:
        BWA + "{sample}.k{kmer}.fpr{fpr}.stats.log"
    benchmark:
        BWA + "{sample}.k{kmer}.fpr{fpr}.stats.bmk"
    conda: "bwa.yml"
    shell:
        """
        n=$(samtools view {input.bam} | cut -f 1 | sort -u | wc -l)
        unmapped=$(samtools view {input.bam} | awk '$3 == "*"' | wc -l)
        mapped=$(samtools view {input.bam} | awk '$3 != "*"'| cut -f 1 | sort -u | wc -l)
        perfect=$(samtools view {input.bam} | awk '$6 ~ /^[0-9]+M$/' | cut -f 1 | sort -u | wc -l)
        unique=$(samtools view {input.bam} | cut -f 1 | sort | uniq -c | awk '$1 == 1' | wc -l)
        multi=$(samtools view {input.bam} | cut -f 1 | sort | uniq -c | awk '$1 > 1' | wc -l)

        echo exons"\t"unmapped"\t"mapped"\t"perfect_match"\t"uniquely_mapped"\t"multimapped > {output.tsv}
        echo $n"\t"$unmapped"\t"$mapped"\t"$perfect"\t"$unique"\t"$multi >> {output.tsv}
        """
