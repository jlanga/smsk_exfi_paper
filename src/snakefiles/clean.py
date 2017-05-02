rule clean:
    shell:
        "rm -r " + " ".join([
            dist, bwa, exons, raw
        ])
