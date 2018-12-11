rule bwa_index:
    input:
        reference = RAW + "genome.fa"
    output:
        expand(
            BWA + "genome.{suffix}",
            suffix = "amb ann bwt pac sa".split()
        ),
        mock = touch(BWA + "genome")
    log:
        BWA + "index.log"
    benchmark:
        BWA + "index.bmk"
    params:
        prefix = BWA + "genome"
    conda: "bwa.yml"
    shell:
        "bwa index "
            "-p {params.prefix} "
            "{input.reference} "
        "2> {log} 1>&2"



rule bwa_map:
    input:
        exons = EXFI + "exons.fa",
        index = BWA + "genome"
    output:
        bam = BWA + "exons.bam"
    log:
        BWA + "map.log"
    benchmark:
        BWA + "index.bmk"
    conda: "bwa.yml"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "-D 0.1 "
            "{input.index} "
            "{input.exons} "
        "| samtools view -Shub "
        "| samtools sort "
            "-l 9 "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"



rule bwa_stats:
    input: bam = BWA + "exons.bam"
    output: tsv = BWA + "stats.tsv"
    log: BWA + "stats.log"
    benchmark: BWA + "stats.bmk"
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
