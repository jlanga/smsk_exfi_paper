rule snp_index:
    input:
        reference = EXFI + \
            "{prefix}.exons.fa"
    output:
        expand(
            SNP + "{prefix}.exons.fa.{suffix}",
            prefix="{prefix}",
            suffix = "amb ann bwt pac sa".split()
        ),
        mock = touch(
            SNP + "{prefix}.exons.fa.bwa"
        )
    log:
        SNP + "{prefix}.exons.fa.bwa.log"
    benchmark:
        SNP + "{prefix}.exons.fa.bwa.bmk"
    params:
        prefix = SNP + "{prefix}.exons.fa"
    conda: "snp.yml"
    shell:
        "bwa index "
            "-p {params.prefix} "
            "{input.reference} "
        "2> {log} 1>&2"


def get_forwards(wildcards):
    """Get the forward read files in wildcards for a species"""
    return samples[samples['sample'] == wildcards.sample].forward.values.tolist()

def get_reverses(wildcards):
    """Get the reverse read files in wildcards for a species"""
    return samples[samples['sample'] == wildcards.sample].reverse.values.tolist()

rule snp_bwa_map_to_exons:
    input:
        index = SNP + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.fa.bwa",
        forwards = get_forwards,
        reverses = get_reverses
    output:
        bam = SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam"
    log:
        SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam.log"
    benchmark:
        SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam.bmk"
    conda: "snp.yml"
    threads: 24
    params:
        index = SNP + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.fa"
    shell:
        "(bwa mem "
            "-t {threads} "
            "-M "
            "{params.index} "
            "<(gzip --decompress --stdout {input.forwards}) "
            "<(gzip --decompress --stdout {input.reverses}) "
        "| samtools view "
            "-q 20 "
            "-F 0x0004  `# read unmapped. Throw away` "
            "-F 0x0008  `# mate unmapped. Throw away` "
            "-u "
            "- "
        "| samtools sort "
            "-l 9 "
            "-@ 24 "
            "-o {output.bam} "
            "/dev/stdin ) "
        "2> {log} 1>&2"




rule snp_call:
    input:
        bam =  SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bam",
        reference = EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.fa"
    output:
        vcf = SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.vcf.gz"
    threads:
        4
    conda: "snp.yml"
    log:
        SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.vcf.gz.log"
    benchmark:
        SNP + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.vcf.gz.bmk"
    shell:
        "(bcftools mpileup "
            "-R {} "
            "-O u "
            "-f {input.reference} "
            "{input.bam} "
        "| bcftools call "
            "-vmO z "
            "-o {output.vcf}) "
        "2> {log} 1>&2"
