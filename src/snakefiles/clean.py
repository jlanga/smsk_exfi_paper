rule clean:
    shell:
        "if [ -d results ]; then "
            "rm -ri results; "
        "fi"
