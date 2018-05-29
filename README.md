# exfi_validation: A Snakemake workflow to validate exfi's performance

[![Build Status](https://travis-ci.org/jlanga/exfi_validation.svg?branch=master)](https://travis-ci.org/jlanga/exfi_validation)

## 1. Description

Workflow to:

- Find putative exons (exfi, biobloomtools, abyss-bloom, bedtools)
- Map agaist different exonic, transcriptomic and genomic references (bwa + samtools)
- Check precision/recall values (bedtools)
- Make reports

## 2. First steps

Follow the contents of the `.travis.yml` file:

0. Install (ana|mini)conda

- [Anaconda](https://www.continuum.io/downloads)

- [miniconda](http://conda.pydata.org/miniconda.html)

1. Clone this repo and install

    ```sh
    git clone --recursive https://github.com/jlanga/exfi_validation.git
    cd exfi_validation/
    bash bin/install/from_src.sh
    ```

2. Execute the pipeline:

    ```sh
    snakemake --use-conda -j
    ```

Once you know it works, edit the `config.yaml` with paths to your data


## 3. File organization

The hierarchy of the folder is the one described in [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424):

```
smsk
├── bin: your binaries, scripts, installation and virtualenv related files.
├── data: raw data, hopefully links to backup data.
├── README.md
├── results: processed data.
|   ├── raw: links to raw data. For your security.
|   ├── exfi: products from the exfi pipeline: bloom filter, gfa with the splice graph, exons.fa with just the exons
|   ├── pr: precision/recall values
|   └── bwa: exons vs reference and some stats.
└── src: additional source code: snakefiles, env files, submodules to other repos.
```


## Bibliography

- [A Quick Guide to Organizing Computational Biology Projects](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1000424)

- [Snakemake—a scalable bioinformatics workflow engine](http://bioinformatics.oxfordjournals.org/content/28/19/2520)
