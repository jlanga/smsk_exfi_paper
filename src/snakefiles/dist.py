rule dist_plot:
    input:
        raw_fai = exfi + "exons.fa.fai",
        true_fai = raw + "exome.fa.fai"
    output:
        dist + "exon_histogram.pdf",
        dist + "exon_density.pdf"
    log:
        dist + "plot.log"
    benchmark:
        dist + "plot.json"
    shell:
        "Rscript src/plot_exon_length_distribution.R 2> {log}"
