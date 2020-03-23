def get_reads(wildcards):
    sample = wildcards.sample
    reads_nested = (
        samples
        [(samples["sample"] == sample)]
        [["forward", "reverse"]]
        .values
        .tolist()
    )

    reads_flat = [item for sublist in reads_nested for item in sublist]

    return reads_flat



rule chopstitch_review_build_bloom:
    input:
        get_reads
    output:
        bf = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.full.unbaited.bf",
        inf = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.full.unbaited.inf"
    params:
        kmer = "{kmer}",
        fpr1 = "{fpr1}",
        fpr2 = "{fpr2}",
        bf = "Bfilter.bf",
        inf = "Bfilter.inf"
    threads:
        24  # All
    log:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.full.unbaited.bf.log"
    benchmark:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.full.unbaited.bf.bmk"
    shell:
        """
        ./bin/CreateBloom \
            -t {threads} \
            -k {params.kmer} \
            --fpr1 {params.fpr1} \
            --fpr2 {params.fpr2} \
            {input} \
            2> {log}

        mv {params.bf} {output.bf}
        mv {params.inf} {output.inf}
        """



rule chopstitch_review_find_exons:
    input:
        fa = RAW + "{sample}.rna.fa",
        bf = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.full.unbaited.bf",
        inf = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.full.unbaited.inf"
    output:
        exons = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.exons.fa",
        processed_exons = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.processedexons.fa"
    params:
        bf = "Bfilter.bf",
        inf = "Bfilter.inf",
        exons = "exons.fa",
        processed_exons = "processed_exons.fa"
    shadow:
        "shallow"
    log:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.find.exons.log"
    benchmark:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.find.exons.bmk"
    shell:
        """
        ln -s {input.bf} {params.bf}
        ln -s {input.inf} {params.inf}

        ./bin/FindExons \
            --input-bloom {params.bf} \
            --lsplicesignals AG,TG,AC,GC,GG \
            --rsplicesignals GT,TT,AT \
            --allexons \
            {input.fa} \
        2> {log}

        mv {params.exons} {output.exons}
        mv {params.processed_exons} {output.processed_exons}
        """



rule pr_chopstitch_review_exons_to_bed3:
    input:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.exons.fa"
    output:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bed"
    log:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.log"
    benchmark:
        CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bmk"
    conda:
        "pr.yml"
    shell:
        'bash src/chopstitch_exons_to_bed3.sh {input} '
        '| sort -k 1,1 -k2,2n '
        '> {output} 2> {log}'


rule pr_chopstitch_review:
    input:
        obs = RAW + "{sample}.gff3",
        pred = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bed",
        fasta = RAW + "{sample}.rna.fa"
    output:
        PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv"
    params:
        identity = "{identity}",
        gff_type = get_gff_type
    conda:
        "pr.yml"
    log:
        PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv.log"
    benchmark:
        PR_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.{identity}.tsv.bmk"
    conda: "pr.yml  "
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


rule bwa_map_chopstitch_review:
    input:
        exons = CHOPSTITCH_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.exons.fa",
        index = BWA + "{sample}.bwa"
    output:
        bam = BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bam"
    log:
        BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bam.log"
    benchmark:
        BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bam.bmk"
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



rule bwa_stats_chopstitch_review:
    input:
        bam = BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.bam",
        fa = RAW + "{sample}.rna.fa"
    output: BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.tsv"
    log: BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.log"
    benchmark: BWA_R + "{sample}.k_{kmer}.fpr1_{fpr1}.fpr2_{fpr2}.chopstitch.stats.bmk"
    conda: 'bwa.yml'
    shell: 'bash src/compute_bwa_stats.sh {input} > {output} 2> {log}'