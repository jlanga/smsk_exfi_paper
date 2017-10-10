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
    shell:
        "build_baited_bloom_filter "
            "--input-fasta {input.transcriptome} "
            "--kmer {params.kmer} "
            "--bloom-size={params.size} "
            "--threads={threads} "
            "--levels={params.levels} "
            "--output-bloom {output.bloom_filter} "
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
        max_fp_bases = config["exfi"]["max_fp_bases"]
    threads:
        1
    log:
        exfi + "build_splice_graph.log"
    benchmark:
        exfi + "build_splice_graph.json"
    shell:
        "build_splicegraph "
            "--input-fasta {input.transcriptome} "
            "--input-bloom {input.bloom_filter} "
            "--kmer {params.kmer} "
            "--max-fp-bases {params.max_fp_bases} "
            "--output-gfa {output.gfa} "
        "2> {log}"



rule exfi_gfa_to_exons:
    input:
        gfa = exfi + "splice_graph.gfa"
    output:
        exons = exfi + "exons.fa"
    threads:
        1
    params:
        extra = config["exfi"]["gfa_to_exons_extra"]
    log:
        exfi + "filter_by_length.log"
    benchmark:
        exfi + "filter_by_length.json"
    shell:
        "gfa_to_exons "
            "--input-gfa {input.gfa} "
            "--output-fasta {output.exons} "
            "{params.extra} "
        "2> {log}"



rule exfi_gfa_to_gapped_transcript:
    input:
        gfa = exfi + "splice_graph.fa"
    output:
        transcripts = exfi + "gapped_transcripts.fa"
    threads:
        1
    params:
        extra = config["exfi"]["gfa_to_gapped_transcript_extra"]
    log:
        exfi + "filter_by_length.log"
    benchmark:
        exfi + "filter_by_length.json"
    shell:
        "gfa_to_gapped_transcript "
            "--input-gfa {input.gfa} "
            "--output-fasta {output.transcripts} "
            "{params.extra} "
        "2> {log}"