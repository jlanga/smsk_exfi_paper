rule exfi_build_baited_bloom_filter:
    input:
        transcriptome = RAW + "assembly.fa",
        fq_gz = expand(
            RAW + "{sample}_{end}.fq.gz",
            sample = SAMPLES,
            end = "1 2".split()
        )
    output:
        bloom_filter = protected(
            expand(
                EXFI + "k{kmer}_l{levels}_m{size}.bloom",
                kmer   = params["exfi"]["kmer"],
                levels = params["exfi"]["levels"],
                size   = params["exfi"]["size"]
            )
        )
    params:
        kmer   = params["exfi"]["kmer"],
        levels = params["exfi"]["levels"],
        size   = params["exfi"]["size"]
    threads:
        32
    log:
        EXFI + "build_baited_bloom_filter.log"
    benchmark:
        EXFI + "build_baited_bloom_filter.json"
    conda: "exfi.yml"
    shell:
        "build_baited_bloom_filter "
            "--input-fasta {input.transcriptome} "
            "--kmer {params.kmer} "
            "--bloom-size={params.size} "
            "--threads={threads} "
            "--levels={params.levels} "
            "--output-bloom {output.bloom_filter} "
            "--verbose "
            "{input.fq_gz} "
        "2> {log}"



rule exfi_build_splice_graph:
    input:
        transcriptome = RAW + "assembly.fa",
        fai = RAW + "assembly.fa.fai",
        bloom_filter = expand(
            EXFI + "k{kmer}_l{levels}_m{size}.bloom",
            kmer   = params["exfi"]["kmer"],
            levels = params["exfi"]["levels"],
            size   = params["exfi"]["size"]
        )
    output:
        gfa = EXFI + "splice_graph.gfa"
    params:
        kmer = params["exfi"]["kmer"],
        max_fp_bases = params["exfi"]["max_fp_bases"],
        max_overlap = params["exfi"]["max_overlap"],
        max_gap_size =  params["exfi"]["max_gap_size"]
    threads: 32
    log: EXFI + "build_splice_graph.log"
    benchmark: EXFI + "build_splice_graph.json"
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
        gfa = EXFI + "splice_graph.gfa"
    output:
        exons = EXFI + "exons.fa"
    params:
        extra = params["exfi"]["gfa1_to_exons_extra"]
    log:
        EXFI + "gfa_to_exons.log"
    benchmark:
        EXFI + "gfa_to_exons.json"
    conda: "exfi.yml"
    shell:
        "gfa1_to_exons "
            "--input-gfa {input.gfa} "
            "--output-fasta {output.exons} "
            "--verbose "
            "{params.extra} "
        "2> {log}"



rule exfi_gfa_to_gapped_transcript:
    input:
        gfa = EXFI + "splice_graph.gfa"
    output:
        transcripts = EXFI + "gapped_transcripts.fa"
    params:
        extra = params["exfi"]["gfa1_to_gapped_transcript_extra"]
    log:
        EXFI + "gfa_to_gapped_transcripts.log"
    benchmark:
        EXFI + "gfa_to_gapped_transcripts.json"
    conda: "exfi.yml"
    shell:
        "gfa1_to_gapped_transcripts "
            "--input-gfa {input.gfa} "
            "--output-fasta {output.transcripts} "
            "--verbose "
            "{params.extra} "
        "2> {log}"
