#!/usr/bin/env Rscript

library(tidyverse)

raw_exons <- "results/exons/raw.fa.fai" %>%
    read_tsv(
        file = .,
        col_names = FALSE
    ) %>%
    select(exon_id = X1, exon_length = X2) %>%
    mutate(label = "Raw exons")

length_exons <- "results/exons/filtered_by_length.fa.fai" %>%
    read_tsv(
        file = .,
        col_names = FALSE
    ) %>%
    select(exon_id = X1, exon_length = X2) %>%
    mutate(label = "Filtered  by length exons")

extensibility_exons <- "results/exons/filtered_by_extensibility.fa.fai" %>%
    read_tsv(
        file = .,
        col_names = FALSE
    ) %>%
    select(exon_id = X1, exon_length = X2) %>%
    mutate(label = "Filtered by extensibility exons")

true_exons <- "results/raw/exome_reduced.fa.fai"  %>%
    read_tsv(
        file = .,
        col_names = FALSE
    ) %>%
    select(exon_id = X1, exon_length = X2) %>%
    mutate(label = "True exons")

exons <- rbind(raw_exons, length_exons, extensibility_exons, true_exons)

library(ggplot2)

q <- ggplot(data = exons, aes(x = exon_length, fill= label, color = label)) +
    geom_histogram(position="dodge", bins= 30) +
    scale_x_log10()
ggsave(q, filename = "results/dist/exon_histogram.pdf", device = "pdf")

q <- ggplot(data = exons, aes(x = exon_length, color = label)) +
    geom_density() +
    scale_x_log10()
ggsave(q, filename = "results/dist/exon_density.pdf", device = "pdf")
