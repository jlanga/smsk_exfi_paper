#/usr/bin/env Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

infile <- args[1]
outfile <- file('/dev/stdout', 'w')


gff3_to_bed <- function(path){
    
    gff_columns <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attribute")
    
    raw <- path %>%
        read_tsv(
            col_names = gff_columns,
            col_types = c("ccciicccc"),
            comment = "#",
            progress = FALSE,
            na = "."
        ) %>%
        select(-c(seqid, source, score, phase)) %>%
        filter(type == "cDNA_match") %>%
        select(-type) %>%
        extract(
            col = attribute,
            into = "transcript_id",
            regex = 'Target=([A-Za-z0-9_.]+)',
            remove = TRUE,
            convert = TRUE
        )
    
    positive <- raw %>%
        filter(strand == "+") %>%
        select(-strand) %>%
        group_by(transcript_id) %>%
        arrange(transcript_id, start, end) %>%
        ungroup()
    
    negative <- raw %>%
        filter(strand == "-") %>%
        select(-strand) %>%
        group_by(transcript_id) %>%
        arrange(transcript_id, desc(start), desc(end)) %>%
        ungroup()
    
    rbind(positive, negative) %>%
        mutate(length = end - start + 1) %>%
        group_by(transcript_id) %>%
        mutate(
            transcript_end = cumsum(length),
            transcript_start = transcript_end - length
        ) %>%
        select(
            chrom = transcript_id,
            chromStart = transcript_start,
            chromEnd = transcript_end
        ) %>%
        ungroup() %>%
        distinct()
}

infile %>%
    gff3_to_bed() %>%
    format_tsv(col_names = FALSE) %>%
    cat()

