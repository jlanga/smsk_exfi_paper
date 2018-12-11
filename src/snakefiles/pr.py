rule pr_gff3_to_exon_bed:
    input:
        gff3_gz = RAW + "annotation.gff3.gz"
    output:
        bed = PR + "true_exons.bed"
    log:
        PR + "gff3_to_exon_bed.log"
    benchmark:
        PR + "gff3_to_exon_bed.json"
    conda:
        "pr.yml"
    shell:
        "(Rscript src/gff3_to_exon_bed.R "
            "{input.gff3_gz} "
        "| sort -k1,1 -k2,2n "
        "> {output.bed}) "
        "2> {log}"



rule pr_exons_to_bed:
    input:
        fasta = EXFI + "exons.fa"
    output:
        bed = PR + "pred_exons.bed"
    log:
        PR + "exons_to_bed.log"
    benchmark:
        PR + "exons_to_bed.json"
    conda:
        "pr.yml"
    shell:
        '(grep ^">" {input.fasta} '
        '| tr -d ">" '
        '| tr "-" "\t" '
        '| tr ":" "\t" '
        '| sort -k1,1 -k2,2n '
        '> {output.bed}) '
        '2> {log}'



rule pr_true_positives:
    """Predicted exons in true"""
    input:
        true= PR + "true_exons.bed",
        pred= PR + "pred_exons.bed"
    output:
        bed = PR + "true_positives.bed"
    params:
        fraction_overlap = params["pr"]["fraction_overlap"]
    log:
        PR + "true_positives.log"
    benchmark:
        PR + "true_positives.json"
    conda:
        "pr.yml"
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
        true= PR + "true_exons.bed",
        pred= PR + "pred_exons.bed"
    output:
        bed = PR + "false_positives.bed"
    params:
        fraction_overlap = params["pr"]["fraction_overlap"]
    log:
        PR + "true_positives.log"
    benchmark:
        PR + "true_positives.json"
    conda:
        "pr.yml"
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
        true= PR + "true_exons.bed",
        pred= PR + "pred_exons.bed"
    output:
        bed = PR + "false_negatives.bed"
    params:
        fraction_overlap = params["pr"]["fraction_overlap"]
    log:
        PR + "true_positives.log"
    benchmark:
        PR + "true_positives.json"
    conda:
        "pr.yml"
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
        tp = PR + "true_positives.bed",
        fp = PR + "false_positives.bed",
        fn = PR + "false_negatives.bed"
    output:
        tsv = PR + "pr.tsv"
    log:
        PR + "pr.log"
    benchmark:
        PR + "precision_recall.json"
    shell:
        '(tp=$(wc -l < {input.tp}); '
        'fp=$(wc -l < {input.fp}); '
        'fn=$(wc -l < {input.fn}); '
        'precision=$(echo "$tp / ($tp + $fp)" | bc -l); '
        'recall=$(echo "$tp / ($tp + $fn)" | bc -l); '
        'f1=$(echo "2 * $precision * $recall / ($precision + $recall)" | bc -l); '
        'echo -e tp"\t"fp"\t"fn"\t"precision"\t"recall"\t"f1 > {output.tsv}; '
        'echo -e $tp"\t"$fp"\t"$fn"\t"$precision"\t"$recall"\t"$f1 >> {output.tsv}) '
        '2> {log}'
