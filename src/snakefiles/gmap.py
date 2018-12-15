rule gmap_build:
    input:
        RAW + "{sample}.dna.fa"
    output:
        touch(GMAP + "{sample}/gmap_build.ok")
    log:
        GMAP + "{sample}.build.log"
    benchmark:
        GMAP + "{sample}.build.bmk"
    params:
        sample = "{sample}",
        gmap_dir = GMAP
    conda:
        "gmap.yml"
    shell:
        'gmap_build '
            '--db {params.sample} '
            '--dir {params.gmap_dir} '
            '{input} '
        '2> {log} 1>&2'



rule gmap_align:
    input:
        rna = RAW + "{sample}.rna.fa",
        index = touch(GMAP + "{sample}/gmap_build.ok")
    output:
        GMAP + "{sample}.gff3"
    log:
        GMAP + "{sample}.gff3.log"
    benchmark:
        GMAP + "{sample}.gff3.bmk"
    params:
        sample = "{sample}",
        gmap_dir = GMAP
    conda:
        "gmap.yml"
    threads:
        32
    shell:
        'gmap '
            '--db {params.sample} '
            '--dir {params.gmap_dir} '
            '{input.rna} '
            '-f gff3_match_cdna '
            '--nthreads {threads} '
        '> {output} '
        '2> {log}'
