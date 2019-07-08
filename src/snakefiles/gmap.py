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
        index = GMAP + "{sample}/gmap_build.ok"
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
        24
    shell:
        'gmapl '
            '--db {params.sample} '
            '--dir {params.gmap_dir} '
            '{input.rna} '
            '-f gff3_match_cdna '
            '--nthreads {threads} '
        '> {output} '
        '2> {log}'



rule gmap_extract_exons:
    input:
        gff3 = GMAP + '{sample}.gff3',
        reference = RAW + '{sample}.rna.fa'
    output:
        GMAP + '{sample}.exons.fa'
    log:
        GMAP + '{sample}.exons.fa.log'
    benchmark:
        GMAP + '{sample}.exons.fa.bmk'
    conda: 'gmap.yml'
    shell:
        '(Rscript src/gff3_gmap_to_exon_bed3.R {input} '
        '| bedtools getfasta -fi {input.reference} -bed /dev/stdin -fo {output}'
        ') 2> {log}'
