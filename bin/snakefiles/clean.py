rule clean:
    shell:
        "rm -r " +
            dist + bwa + exons + raw
