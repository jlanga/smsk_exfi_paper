rule dist_plot:
    input:
        raw_fai = exons + "raw.fa.fai",
        length_fai = exons + "filtered_by_length.fa.fai",
        extensibility_fai = exons + "filtered_by_extensibility.fa.fai",
        true_fai = raw + "exome_reduced.fa.fai"
    output:
        dist + "exon_histogram.pdf",
        dist + "exon_density.pdf"
    log:
        dist + "plot.log"
    benchmark:
        dist + "plot.json"
    shell:
        "Rscript src/plot_exon_length_distribution.R 2> {log}"
