rule exons_build_bloom_filter:
    input:
        fq_gz = expand(
            raw + "{sample}_{end}.fq.gz",
            sample = dna_pe,
            end = "1 2".split()
        )
    output:
        bloom_filter = protected(
                expand(
                exons + "k{kmer}_l{levels}_m{size}.bloom",
                kmer   = config["exons"]["kmer"],
                levels = config["exons"]["levels"],
                size   = config["exons"]["size"]
            )
        )
    params:
        kmer   = config["exons"]["kmer"],
        levels = config["exons"]["levels"],
        size   = config["exons"]["size"]
    threads:
        32
    log:
        exons + "build_bloom_filter.log"
    benchmark:
        exons + "build_bloom_filter.json"
    shell:
        "abyss-bloom build "
            "--kmer={params.kmer} "
            "--verbose "
            "--bloom-size={params.size} "
            "--threads={threads} "
            "--levels={params.levels} "
            "{output.bloom_filter} "
            "<(pigz --decompress --stdout {input.fq_gz}) "
        "2> {log}"



rule exons_find_exons:
    input:
        transcriptome = raw + "assembly.fa",
        fai = raw + "assembly.fa.fai",
        bloom_filter = expand(
            exons + "k{kmer}_l{levels}_m{size}.bloom",
            kmer   = config["exons"]["kmer"],
            levels = config["exons"]["levels"],
            size   = config["exons"]["size"]
        )
    output:
        exons = exons + "raw.fa"
    params:
        kmer = config["exons"]["kmer"]
    threads:
        1
    log:
        exons + "find_exons.log"
    benchmark:
        exons + "find_exons.json"
    shell:
        "find_exons "
            "--input-fasta {input.transcriptome} "
            "--input-bloom {input.bloom_filter} "
            "--kmer {params.kmer} "
            "--output-fasta {output.exons} "
        "2> {log}"



rule exons_filter_by_length:
    input:
        exons_raw = exons + "raw.fa"
    output:
        exons_out = exons + "filtered_by_length.fa"
    threads:
        1
    params:
        kmer = config["exons"]["kmer"],
        trim_left = config["exons"]["trim_left"],
        trim_right = config["exons"]["trim_right"]
    log:
        exons + "filter_by_length.log"
    benchmark:
        exons + "filter_by_length.json"
    shell:
        "filter_exons_by_length "
            "--input-fasta {input.exons_raw} "
            "--minimum-length {params.kmer} "
            "--trim-left {params.trim_left} "
            "--trim-right {params.trim_right} "
            "--output-fasta {output.exons_out} "
        "2> {log}"



rule exons_filter_by_extensibility:
    input:
        exons_raw = exons + "raw.fa",
        bloom_filter = expand(
            exons + "k{kmer}_l{levels}_m{size}.bloom",
            kmer   = config["exons"]["kmer"],
            levels = config["exons"]["levels"],
            size   = config["exons"]["size"]
        )
    output:
        exons_out = exons + "filtered_by_extensibility.fa"
    threads:
        1
    params:
        kmer = config["exons"]["kmer"],
    log:
        exons + "filter_by_extensibility.log"
    benchmark:
        exons + "filter_by_extensibility.json"
    shell:
        "filter_exons_by_extensibility "
            "--input-fasta {input.exons_raw} "
            "--input-bloom {input.bloom_filter} "
            "--kmer {params.kmer} "
            "--output-fasta {output.exons_out} "
        "2> {log}"