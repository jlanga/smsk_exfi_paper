rule tuxedo2_hisat2_build:
    """Build HISAT2 index"""
    input:
        fa = RAW + "{sample}.dna.fa"
    output:
        suffixes=protected(expand(
            TUXEDO2 + "{sample}.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split(),
            sample="{sample}"
        )),
        mock=protected(touch(TUXEDO2 + "{sample}.ok"))
    params:
        prefix = TUXEDO2 + "{sample}"
    threads:
        1000
    log:
        TUXEDO2 + "hisat2_build_{sample}.log"
    benchmark:
        TUXEDO2 + "hisat2_build_{sample}.bmk"
    conda:
        "tuxedo2.yml"
    shell:
        "hisat2-build "
            "-p {threads} "
            "{input.fa} "
            "{params.prefix} "
        "2> {log} 1>&2"



rule tuxedo2_hisat2_align:
    """Map reads with hisat2 and convert to bam on the fly"""
    input:
        mock = TUXEDO2 + "{sample}.ok",
        index_files=expand(
            TUXEDO2 + "{sample}.{extension}.ht2",
            extension="1 2 3 4 5 6 7 8".split(),
            sample="{sample}"
        ),
        transcriptome = RAW + "{sample}.rna.fa",
        reference=RAW + "{sample}.dna.fa"
    output:
        cram = protected(TUXEDO2 + "{sample}.cram")
    threads:
        4  # Samtools is the bottleneck
    log:
        TUXEDO2 + "hisat2_align_{sample}.log"
    benchmark:
        TUXEDO2 + "hisat2_align_{sample}.bmk"
    params:
        index_prefix = TUXEDO2 + "{sample}"
    conda:
        "tuxedo2.yml"
    shell:
        "(hisat2 "
            "--threads {threads} "
            "--dta "
            "-f "
            "-x {params.index_prefix} "
            "-U {input.transcriptome} "
            "-S /dev/stdout "
        "| samtools sort "
            "-@ {threads} "
            "-l 9 "
            "--output-fmt CRAM "
            "--reference {input.reference} "
            "-o {output.cram} "
            "- ) "
        "2> {log} 1>&2"




rule tuxedo2_stringtie_assemble:
    """assemble transcripts with stringtie"""
    input:
        cram = TUXEDO2 + "{sample}.cram",
        reference_fa = RAW + "{sample}.dna.fa"
    output:
        gtf = TUXEDO2 + "{sample}.gtf"
    params:
        label = "{sample}"
    log:
        TUXEDO2 + "stringtie_{sample}.log"
    benchmark:
        TUXEDO2 + "stringtie_{sample}.bmk"
    conda:
        "tuxedo2.yml"
    shell:
        "(samtools view "
            "-h "
            "--reference {input.reference_fa} "
            "{input.cram} "
        "| stringtie "
            "-p {threads} "
            "-o {output.gtf} "
            "-l {params.label} "
            "-f 0 "  # Minimum abundance
            "-m 30 "  # Minimum transcript length
            "-c 1 "  # Minimum coverage
            "-M 1 "
            "- )"
        "2> {log} 1>&2"
