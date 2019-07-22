def get_gff_type(wildcards):
    return features[wildcards.sample]["reference_type"]


rule pr_exfi:
    input:
        obs = RAW + "{sample}.gff3",
        pred = EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gfa",
        fasta = RAW + "{sample}.rna.fa"
    output:
        PR + \
            "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi."
            "{identity}.tsv"
    params:
        identity = "{identity}",
        gff_type = get_gff_type
    conda:
        "pr.yml"
    log:
        PR + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi."
        "{identity}.tsv.log"
    benchmark:
        PR + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exfi."
        "{identity}.tsv.bmk"
    shell:
        'compare_to_gff3 '
            '--input-splice-graph {input.pred} '
            '--input-fasta {input.fasta} '
            '--input-gff3 {input.obs} '
            '--simmilarity-fraction {params.identity} '
            '--type-gff3 {params.gff_type} '
            '--verbose '
        '>{output} '
        '2>{log}'



rule pr_chopstitch:
    input:
        obs = RAW + "{sample}.gff3",
        pred = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bed",
        fasta = RAW + "{sample}.rna.fa"
    output:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv"
    params:
        identity = "{identity}",
        gff_type = get_gff_type
    conda:
        "pr.yml"
    log:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv.log"
    benchmark:
        PR + "{sample}.k{kmer}.fpr{fpr}.chopstitch.{identity}.tsv.bmk"
    shell:
        'compare_to_gff3 '
            '--input-splice-graph {input.pred} '
            '--input-fasta {input.fasta} '
            '--input-gff3 {input.obs} '
            '--simmilarity-fraction {params.identity} '
            '--type-gff3 {params.gff_type} '
            '--verbose '
        '>{output} '
        '2>{log}'


rule pr_gmap:
    input:
        obs = RAW + "{sample}.gff3",
        pred = GMAP + "{sample}.gff3",
        fasta = RAW + "{sample}.rna.fa"
    output:
        PR + "{sample}.gmap.{identity}.tsv"
    params:
        identity = "{identity}",
        gff_type = get_gff_type
    conda:
        "pr.yml"
    log:
        PR + "{sample}.gmap.{identity}.tsv.log"
    benchmark:
        PR + "{sample}.gmap.{identity}.tsv.bmk"
    shell:
        'compare_to_gff3 '
            '--input-splice-graph {input.pred} '
            '--input-gff3 {input.obs} '
            '--input-fasta {input.fasta} '
            '--simmilarity-fraction {params.identity} '
            '--type-gff3 {params.gff_type} '
            '--verbose '
        '>{output} '
        '2>{log}'
