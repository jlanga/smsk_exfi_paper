library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

infile <- args[1]
outfile <- args[2]
gff3_to_bed <- function(path){
    
    path %>%
        read_tsv(
            col_names = c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attribute"),
            col_types = c("cccddcccc"),
            comment = "#",
            progress = FALSE,
            na = "."
        ) %>%
        select(-c(source, score, phase)) %>%
        filter(type == "exon") %>%
        select(-type) %>%
        extract(
            col = attribute,
            into = c("transcript_id", "exon_id", "rank"),
            regex = "Parent=transcript:(.+);Name=.*;exon_id=(.+);rank=(.+);version=.+",
            remove=TRUE,
            convert=TRUE
        ) %>%
        arrange(transcript_id, rank) %>%
        mutate(length = end - start + 1) %>% 
        group_by(transcript_id) %>%
        mutate(
            transcript_end = cumsum(length) %>% as.integer(),
            transcript_start = transcript_end - length %>% as.integer(),
        ) %>%
        select(
            chrom = transcript_id,
            chromStart = transcript_start,
            chromEnd = transcript_end,
            name = exon_id
        )
}

infile %>%
    gff3_to_bed() %>%
    format_tsv(col_names=FALSE) %>%
    cat()

