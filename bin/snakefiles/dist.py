rule dist_plot:
    input:
        raw_fai = exons + "raw.fa.fai",
        filtered_fai = exons + "filtered.fa.fai",
        true_fai = config["reference"]["exome"] + ".fai"
    output:
        dist + "exon_histogram.pdf",
        dist + "exon_density.pdf"
    log:
        dist + "plot.log"
    benchmark:
        dist + "plot.json"
    shell:
        "Rscript bin/plot_exon_length_distribution.R 2> {log}"