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



rule chopstitch_build_bloom:
    input:
        get_reads
    output:
        bf = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.full.unbaited.bf",
        inf = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.full.unbaited.inf"
    params:
        kmer = "{kmer}",
        fpr = "{fpr}",
        bf = "Bfilter.bf",
        inf = "Bfilter.inf"
    threads:
        24  # All
    log:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.full.unbaited.bf.log"
    benchmark:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.full.unbaited.bf.bmk"
    shell:
        """
        ./bin/CreateBloom \
            -t {threads} \
            -k {params.kmer} \
            --fpr1 {params.fpr} \
            --fpr2 {params.fpr} \
            {input} \
            2> {log}

        mv {params.bf} {output.bf}
        mv {params.inf} {output.inf}
        """



rule chopstitch_find_exons:
    input:
        fa = RAW + "{sample}.rna.fa",
        bf = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.full.unbaited.bf",
        inf = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.full.unbaited.inf"
    output:
        exons = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.exons.fa",
        processed_exons = CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.processedexons.fa"
    params:
        bf = "Bfilter.bf",
        inf = "Bfilter.inf",
        exons = "exons.fa",
        processed_exons = "processed_exons.fa"
    shadow:
        "shallow"
    log:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.find.exons.log"
    benchmark:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.find.exons.bmk"
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



rule pr_chopstitch_exons_to_bed3:
    input:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.exons.fa"
    output:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bed"
    log:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.chopstitch.log"
    benchmark:
        CHOPSTITCH + "{sample}.k{kmer}.fpr{fpr}.chopstitch.bmk"
    conda:
        "chopstitch.yml"
    shell:
        'bash src/chopstitch_exons_to_bed3.sh {input} '
        '| sort -k 1,1 -k2,2n '
        '> {output} 2> {log}'
