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



rule exfi_build_unbaited_bloom_filter:
    input:
        fq_gz = get_reads
    output:
        bloom_filter = protected(
            expand(
                EXFI + "{sample}.k{kmer}.l{levels}.m{size}.100.unbaited.bf",
                sample = "{sample}",
                kmer = "{kmer}",
                levels = "{levels}",
                size = "{size}"
            )
        )
    params:
        kmer   = "{kmer}",
        levels = "{levels}",
        size   = "{size}"
    threads:
        24
    log:
        EXFI + "bbf.{sample}.k{kmer}.l{levels}.m{size}.100.unbaited.log"
    benchmark:
        EXFI + "bbf.{sample}.k{kmer}.l{levels}.m{size}.100.unbaited.bmk"
    conda: "exfi.yml"
    shell:
        "abyss-bloom build "
            "--kmer {params.kmer} "
            "--bloom-size {params.size}G "
            "--threads {threads} "
            "--levels {params.levels} "
            "--verbose "
            "{output.bloom_filter} "
            "{input.fq_gz} "
        "2> {log}"



rule exfi_build_baited_bloom_filter_with_sampling:
    input:
        reads = get_reads,
        transcriptome = get_transcriptome
    output:
        bloom_filter = protected(
            expand(
                EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.baited.bf",
                sample="{sample}",
                kmer="{kmer}",
                levels="{levels}",
                size="{size}",
                sampling="{sampling}"
            )
        )
    params:
        sampling = "{sampling}",
        kmer = "{kmer}",
        levels = "{levels}",
        size = "{size}",
    threads:
        96  # Block everything
    log:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.baited.bf.log"
    benchmark:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.baited.bf.bmk"
    conda:
        "exfi.yml"
    shell:
        """
        tempfile=$(mktemp --dry-run)
        mkfifo $tempfile

        ( if [[ {params.sampling} -ge 100 ]]; then
            pigz -dc {input.reads} > $tempfile
        else
            seqtk sample \
                -s 1 \
                <(pigz -dc {input.reads}) \
                $(echo "{params.sampling} / 100" | bc -l) \
            > $tempfile
        fi ) &

        build_baited_bloom_filter \
            --input-fasta {input.transcriptome} \
            --kmer {params.kmer} \
            --bloom-size {params.size}G \
            --threads {threads} \
            --levels {params.levels} \
            --output-bloom {output.bloom_filter} \
            --verbose \
            $tempfile \
        2> {log}

        rm $tempfile
        """


rule exfi_build_splice_graph:
    input:
        transcriptome = RAW + "{sample}.rna.fa",
        fai = RAW + "{sample}.rna.fa.fai",
        bloom_filter = protected(
            expand(
                EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.bf",
                sample="{sample}",
                kmer="{kmer}",
                levels="{levels}",
                size="{size}",
                sampling="{sampling}",
                type="{type}"
            )
        )
    output:
        gfa = EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gfa"
    params:
        kmer = "{kmer}",
        max_fp_bases = params["exfi"]["max_fp_bases"],
        max_overlap = params["exfi"]["max_overlap"],
        max_gap_size =  params["exfi"]["max_gap_size"]
    threads: 32
    log:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gfa.log"
    benchmark:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gfa.bmk"
    conda: "exfi.yml"
    shell:
        "build_splice_graph "
            "--input-fasta {input.transcriptome} "
            "--input-bloom {input.bloom_filter} "
            "--kmer {params.kmer} "
            "--correct "
            "--polish "
            "--threads {threads} "
            "--max-overlap {params.max_overlap} "
            "--max-gap-size {params.max_gap_size} "
            "--max-fp-bases {params.max_fp_bases} "
            "--output-gfa {output.gfa} "
            "--verbose "
        "2> {log}"



rule exfi_gfa_to_exons:
    input:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gfa"
    output:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.fa"
    params:
        extra = params["exfi"]["gfa1_to_exons_extra"]
    log:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.log"
    benchmark:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.exons.bmk"
    conda: "exfi.yml"
    shell:
        "gfa1_to_exons "
            "--input-gfa {input} "
            "--output-fasta {output} "
            "--verbose "
            "{params.extra} "
        "2> {log}"



rule exfi_gfa_to_gapped_transcript:
    input:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gfa"
    output:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gapped.fa"
    params:
        extra = params["exfi"]["gfa1_to_gapped_transcript_extra"]
    log:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gapped.log"
    benchmark:
        EXFI + "{sample}.k{kmer}.l{levels}.m{size}.{sampling}.{type}.gapped.bmk"
    conda: "exfi.yml"
    shell:
        "gfa1_to_gapped_transcripts "
            "--input-gfa {input.gfa} "
            "--output-fasta {output.transcripts} "
            "--verbose "
            "{params.extra} "
        "2> {log}"
