rule exfi_build_baited_bloom_filter:
    input:
        transcriptome = raw + "assembly.fa",
        fq_gz = expand(
            raw + "{sample}_{end}.fq.gz",
            sample = dna_pe,
            end = "1 2".split()
        )
    output:
        bloom_filter = protected(
            expand(
                exfi + "k{kmer}_l{levels}_m{size}.bloom",
                kmer   = config["exfi"]["kmer"],
                levels = config["exfi"]["levels"],
                size   = config["exfi"]["size"]
            )
        )
    params:
        kmer   = config["exfi"]["kmer"],
        levels = config["exfi"]["levels"],
        size   = config["exfi"]["size"]
    threads:
        32
    log:
        exfi + "build_baited_bloom_filter.log"
    benchmark:
        exfi + "build_baited_bloom_filter.json"
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
        transcriptome = raw + "assembly.fa",
        fai = raw + "assembly.fa.fai",
        bloom_filter = expand(
            exfi + "k{kmer}_l{levels}_m{size}.bloom",
            kmer   = config["exfi"]["kmer"],
            levels = config["exfi"]["levels"],
            size   = config["exfi"]["size"]
        )
    output:
        gfa = exfi + "splice_graph.gfa"
    params:
        kmer = config["exfi"]["kmer"],
        max_fp_bases = config["exfi"]["max_fp_bases"],
        max_overlap = config["exfi"]["max_overlap"],
        max_gap_size =  config["exfi"]["max_gap_size"]
    threads: 32
    log: exfi + "build_splice_graph.log"
    benchmark: exfi + "build_splice_graph.json"
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
        gfa = exfi + "splice_graph.gfa"
    output:
        exons = exfi + "exons.fa"
    params:
        extra = config["exfi"]["gfa1_to_exons_extra"]
    log:
        exfi + "gfa_to_exons.log"
    benchmark:
        exfi + "gfa_to_exons.json"
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
        gfa = exfi + "splice_graph.gfa"
    output:
        transcripts = exfi + "gapped_transcripts.fa"
    params:
        extra = config["exfi"]["gfa1_to_gapped_transcript_extra"]
    log:
        exfi + "gfa_to_gapped_transcripts.log"
    benchmark:
        exfi + "gfa_to_gapped_transcripts.json"
    conda: "exfi.yml"
    shell:
        "gfa1_to_gapped_transcripts "
            "--input-gfa {input.gfa} "
            "--output-fasta {output.transcripts} "
            "--verbose "
            "{params.extra} "
        "2> {log}"
