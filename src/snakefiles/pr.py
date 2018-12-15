rule pr_gff3_to_exon_bed3:
    input:
        RAW + "{sample}.gff3"
    output:
        PR + "{sample}.true.bed"
    log:
        PR + "{sample}.true.log"
    benchmark:
        PR + "{sample}.true.bmk"
    conda:
        "pr.yml"
    shell:
        """
        (Rscript src/gff3_to_exon_bed3.R \
            {input} \
        | sort -k 1,1 -k2,2n \
        > {output}) \
        2> {log}
        """



rule pr_exfi_exons_to_bed3:
    input:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.fa"
    output:
        PR + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bed"
    log:
        PR + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.log"
    benchmark:
        PR + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bmk"
    conda:
        "pr.yml"
    shell:
        'bash src/exfi_exons_to_bed3.sh {input} '
        '| sort -k 1,1 -k2,2n '
        '> {output} 2> {log}'


rule pr_chopstitch_exons_to_bed3:
    input:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.processedexons.fa"
    output:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bed"
    log:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.log"
    benchmark:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bmk"
    conda:
        "pr.yml"
    shell:
        'bash src/chopstitch_exons_to_bed3.sh {input} '
        '| sort -k 1,1 -k2,2n '
        '> {output} 2> {log}'



rule pr_exfi:
    input:
        obs = PR + "{sample}.true.bed",
        pred = PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.bed"
    output:
        PR + \
         "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi.{identity}.tsv"
    params:
        identity = "{identity}"
    conda:
        "pr.yml"
    shell:
        'bash src/compute_pr.sh '
            '{input.pred} '
            '{input.obs} '
            '{params.identity} '
        '> {output}'



rule pr_chopstitch:
    input:
        obs = PR + "{sample}.true.bed",
        pred = PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bed"
    output:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv"
    params:
        identity = "{identity}"
    conda:
        "pr.yml"
    shell:
        'bash src/compute_pr.sh '
            '{input.pred} '
            '{input.obs} '
            '{params.identity} '
        '> {output}'



rule pr_gmap_gff3_to_bed:
    input:
        gff3 = GMAP + "{sample}.gff3"
    output:
        bed = PR + "{sample}.gmap.bed"
    conda:
        "pr.yml"
    log:
        PR + "{sample}.gmap.bed.log"
    benchmark:
        PR + "{sample}.gmap.bed.bmk"
    shell:
        'Rscript src/gff3_gmap_to_exon_bed3.R {input} '
        '| sort -k 1,1 -k2,2n '
        '> {output} 2> {log}'


rule pr_gmap:
    input:
        obs = PR + "{sample}.true.bed",
        pred = PR + "{sample}.gmap.bed"
    output:
        PR + "{sample}.gmap.{identity}.tsv"
    params:
        identity = "{identity}"
    conda:
        "pr.yml"
    shell:
        'bash src/compute_pr.sh '
            '{input.pred} '
            '{input.obs} '
            '{params.identity} '
        '> {output}'
