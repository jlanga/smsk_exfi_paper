rule pr_gff3_to_exon_bed:
    input:
        gff3_gz = raw + "annotation.gff3.gz"
    output:
        bed = pr + "true_exons.bed"
    log:
        pr + "gff3_to_exon_bed.log"
    benchmark:
        pr + "gff3_to_exon_bed.json"
    shell:
        "(Rscript src/gff3_to_exon_bed.R "
            "{input.gff3_gz} "
        "| bedtools sort "
        "> {output.bed}) "
        "2> {log}"



rule pr_exons_to_bed:
    input:
        fasta = exfi + "exons.fa"
    output:
        bed = pr + "predicted_exons.bed"
    log:
        pr + "exons_to_bed.log"
    benchmark:
        pr + "exons_to_bed.json"
    shell:
        '(grep ^">" {input.fasta} '
        '| tr -d ">" '
        '| tr "-" "\t" '
        '| tr ":" "\t" '
        '| bedtools sort '
        '> {output.bed}) '
        '2> {log}'



rule pr_true_positives:
    """Predicted exons in true"""
    input:
        true= pr + "true_exons.bed",
        pred= pr + "pred_exons.bed"
    output:
        bed = pr + "true_positives.bed"
    params:
        fraction_overlap = config["pr"]["fraction_overlap"]
    log:
        pr + "true_positives.log"
    benchmark:
        pr + "true_positives.json"
    shell:
        "bedtools intersect "
            "-a {input.pred} "
            "-b {input.true} "
            "-f {params.fraction_overlap} "
            "-r "
        "> {output.bed}"



rule pr_compute_false_positives:
    """Predicted exons not in true"""
    input:
        true= pr + "true_exons.bed",
        pred= pr + "pred_exons.bed"
    output:
        bed = pr + "false_positives.bed"
    params:
        fraction_overlap = config["pr"]["fraction_overlap"]
    log:
        pr + "true_positives.log"
    benchmark:
        pr + "true_positives.json"
    shell:
        "bedtools intersect "
            "-a {input.pred} "
            "-b {input.true} "
            "-f {params.fraction_overlap} "
            "-r "
            "-v "
        "> {output.bed}"



rule pr_compute_false_negatives:
    """True exons not in predicted"""
    input:
        true= pr + "true_exons.bed",
        pred= pr + "pred_exons.bed"
    output:
        bed = pr + "false_negatives.bed"
    params:
        fraction_overlap = config["pr"]["fraction_overlap"]
    log:
        pr + "true_positives.log"
    benchmark:
        pr + "true_positives.json"
    shell:
        "bedtools intersect "
            "-a {input.true} "
            "-b {input.pred} "
            "-f {params.fraction_overlap} "
            "-r "
            "-v "
        "> {output.bed}"


rule pr_precision_recall:
    input:
        tp = pr + "true_positives.bed",
        fp = pr + "false_positives.bed",
        fn = pr + "false_negatives.bed"
    output:
        tsv = pr + "pr.tsv"
    log:
        pr + "pr.log"
    benchmark:
        pr + "precision_recall.json"
    shell:
        '(tp=$(wc -l < {input.tp}); '
        'fp=$(wc -l < {input.fp}); '
        'fn=$(wc -l < {input.fn}); '
        'precision=$(echo "$tp / ($tp + $fp)" | bc -l); '
        'recall=$(echo "$tp / ($tp + $fn)" | bc -l); '
        'echo -e tp"\t"fp"\t"fn"\t"precision"\t"recall > {output.tsv}; '
        'echo -e $tp"\t"$fp"\t"$fn"\t"$precision"\t"$recall >> {output.tsv}) '
        '2> {log}'
