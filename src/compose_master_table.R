#!/usr/bin/env Rscript

library(tidyverse)
library(lubridate)

species_df <- data.frame(
    species = c(
        "amex", "char", "char_ref", "drer", "hsap", "haxy", "ttin", "ssal",
        "plam"
    ),
    species_binom = c(
        "Ambystoma mexicanus", "Clupea harengus", "Clupea harengus", 
        "Danio rerio", "Homo sapiens", "Harmonia axyridis", "Tinca tinca",
        "Salmo salar", "Pinus lambertiana"
    ),
    species_plain = c(
        "Axolotl", "Atl. herring de novo", "Atl. herring ref", "Zebrafish",
        "Human", "Ladybug", "Tinca", "Salmon", "Sugar pine"
    ),
    reference_type = c(
        "de novo", "de novo", "ncbi", "ensembl", "ensembl", "de novo", 
        "de novo", "ncbi", "de novo"
    )
)

# bwa stats
bwa_initial <- read_tsv(
        file = 'bwa_initial.stats.tsv',
        col_names = c(
            'file', 'exons', 'unmapped', 'mapped', 'perfect_match',
            'uniquely_mapped', 'multimapped', 'transcripts',
            'transcripts_with_exons'
        )
    ) %>%
    mutate(
        mapped_rate = mapped / exons,
        perfect_rate = perfect_match / exons
    )

bwa_review <- read_tsv(
        file = 'bwa_review.stats.tsv',
        col_names = c(
            'file', 'exons', 'unmapped', 'mapped', 'perfect_match',
            'uniquely_mapped', 'multimapped', 'transcripts',
            'transcripts_with_exons'
        )
    ) %>%
    mutate(
        mapped_rate = mapped / exons,
        perfect_rate = perfect_match / exons
    )

# pr
pr_initial <- read_tsv(
        file = "pr_initial.tsv",
        col_names = c(
            'file', 'type', 'min_similarity', 'true', 'predicted', 'true_positives',
            'false_positives', 'false_negatives', 'precision', 'recall', 'f_1'
        )
    ) %>%
    gather(key = 'metric', value = 'value', precision, recall, f_1) %>%
    mutate(
        metric = paste(metric, type, sep = '_'),
    ) %>%
    select(-type, -min_similarity, -true:-false_negatives) %>%
    spread(key = metric, value = value)

pr_review <- read_tsv(
        file = "pr_review.tsv",
        col_names = c(
            'file', 'type', 'min_similarity', 'true', 'predicted', 'true_positives',
            'false_positives', 'false_negatives', 'precision', 'recall', 'f_1'
        )
    ) %>%
    gather(key = 'metric', value = 'value', precision, recall, f_1) %>%
    mutate(
        metric = paste(metric, type, sep = '_'),
    ) %>%
    select(-type, -min_similarity, -true:-false_negatives) %>%
    spread(key = metric, value = value)






# exfi
exfi_build_times <- read_tsv(
        file = 'exfi.build.bmk',
        col_names = c('file', 'build_time', 'build_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'l', 'm', 'sampling', 'baiting'),
        regex = 'results/exfi/([a-z_]+).k(\\d+).l(\\d+).m(\\d+).(\\d+).(\\w+).bf',
        convert = TRUE
    ) %>%
    mutate(
        build_memory = build_memory / 1024
    )

exfi_build_ir <- read_tsv(file = 'exfi.ir.tsv', col_names = c('file', 'insertion_rate')) %>%
    extract(
        col = file,
        into = c('species', 'k', 'l', 'm', 'sampling', 'baiting'),
        regex = 'results/exfi/([a-z_]+).k(\\d+).l(\\d+).m(\\d+).(\\d+).(\\w+).bf.log',
        convert = TRUE
    )

exfi_build_fpr <- read_tsv('exfi.fpr.tsv', col_names = c('file', 'fpr_percent')) %>%
    extract(
        col = file,
        into = c('species', 'k', 'l', 'm', 'sampling', 'baiting'),
        regex = 'results/exfi/([a-z_]+).k(\\d+).l(\\d+).m(\\d+).(\\d+).(\\w+).bf',
        convert = TRUE
    ) %>%
    extract(
        col = fpr_percent,
        into = 'fpr',
        regex = '([0-9.]+)%',
        convert = TRUE
    ) %>%
    mutate(fpr = fpr / 100.0)
exfi_build_fpr['fpr_type'] = c('fpr1', 'fpr2')
exfi_build_fpr <- exfi_build_fpr %>%
    spread(key = fpr_type, value = fpr)

exfi_build_master <- left_join(exfi_build_fpr, exfi_build_ir) %>%
    mutate(insertion_rate = replace_na(insertion_rate, 1)) %>%
    left_join(exfi_build_times) %>%
    mutate(
        hash_fn1 = 4,
        hash_fn2 = 4
    )

exfi_predict_master <- read_tsv(
        file = 'exfi.predict.bmk',
        col_names = c('file', 'predict_time', 'predict_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'l', 'm', 'sampling', 'baiting'),
        regex = 'results/exfi/([a-z_]+).k(\\d+).l(\\d+).m(\\d+).(\\d+).(\\w+).gfa',
        convert = TRUE
    ) %>%
    mutate(predict_memory = predict_memory / 1024)

exfi_bwa <- bwa_initial %>%
    filter(grepl(pattern = 'exfi', x = file)) %>%
    extract(
        col = file,
        into = c('species', 'k', 'l', 'm', 'sampling', 'baiting'),
        regex = 'results/bwa/([a-z_]+).k(\\d+).l(\\d+).m(\\d+).(\\d+).(\\w+).exfi.stats.tsv',
        convert = TRUE
    )

exfi_pr <- pr_initial %>%
    filter(grepl('exfi', file)) %>%
    extract(
        col = file,
        into = c('species', 'k', 'l', 'm', 'sampling', 'baiting', 'min_similarity'),
        regex = 'results/pr/([a-z_]+).k(\\d+).l(\\d+).m(\\d+).(\\d+).(\\w+).exfi.(\\d+.\\d+).tsv',
        convert = TRUE
    )


exfi_master <- exfi_build_master %>%
    left_join(exfi_predict_master) %>%
    left_join(exfi_pr) %>%
    left_join(exfi_bwa) 

rm(exfi_build_fpr, exfi_build_ir, exfi_build_times, exfi_build_master, 
   exfi_predict_master, exfi_pr, exfi_bwa)



# Chopstitch build
chopstitch_build_initial <- read_tsv(
        file = "chopstitch_build_data_initial.tsv"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'fpr', 'sampling', 'baiting'),
        regex = 'results/chopstitch/([a-z_]+).k(\\d+).fpr([0-9]+.[0-9]+).(\\w+).(\\w+).bf',
        convert = TRUE
    ) %>%
    mutate(sampling = 100, target_fpr1 = fpr, target_fpr2 = fpr) %>%
    select(-fpr)

chopstitch_build_review <- read_tsv(
        file = "chopstitch_build_data_review.tsv"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'target_fpr1', 'target_fpr2', 'sampling', 'baiting'),
        regex = 'results/chopstitch_review/([a-z_]+).k_(\\d+).fpr1_([0-9]+.[0-9]+).fpr2_([0-9]+.[0-9]+).(\\w+).(\\w+).bf',
        convert = TRUE
    ) %>%
    mutate(sampling = 100)

chopstitch_build_times_initial <- read_tsv(
        file = 'chopstitch_initial.bf.bmk', 
        col_names = c('file', 'build_time', 'build_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'fpr', 'sampling', 'baiting'),
        regex = 'results/chopstitch/([a-z_]+).k(\\d+).fpr([0-9]+.[0-9]+).(\\w+).(\\w+).bf',
        convert = TRUE
    ) %>%
    mutate(
        sampling = 100,
        build_memory = build_memory / 1024,
        target_fpr1 = fpr,
        target_fpr2 = fpr
    ) %>%
    select(-fpr)

chopstitch_build_times_review <- read_tsv(
        file = 'chopstitch_review.bf.bmk', 
        col_names = c('file', 'build_time', 'build_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'target_fpr1', 'target_fpr2', 'sampling', 'baiting'),
        regex = 'results/chopstitch_review/([a-z_]+).k_(\\d+).fpr1_([0-9]+.[0-9]+).fpr2_([0-9]+.[0-9]+).(\\w+).(\\w+).bf',
        convert = TRUE
    ) %>%
    mutate(
        sampling = 100,
        build_memory = build_memory / 1024
    )

chopstitch_build_initial <- left_join(chopstitch_build_initial, chopstitch_build_times_initial)
chopstitch_build_review <- left_join(chopstitch_build_review, chopstitch_build_times_review)
chopstitch_build <- bind_rows(chopstitch_build_initial, chopstitch_build_review)
rm(chopstitch_build_initial, chopstitch_build_times_initial, chopstitch_build_review, chopstitch_build_times_review)


# Chopstitch predict
chopstitch_predict_initial <- read_tsv(
        file = 'chopstitch_initial.predict.bmk',
        col_names = c('file', 'predict_time', 'predict_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'fpr'),
        regex = 'results/chopstitch/([a-z_]+).k(\\d+).fpr([0-9]+.[0-9]+).find.exons',
        convert = TRUE
    ) %>%
    mutate(
        predict_memory = predict_memory / 1024,
        target_fpr1 = fpr,
        target_fpr2 = fpr
    ) %>%
    select(-fpr)

chopstitch_predict_review <- read_tsv(
        file = 'chopstitch_review.predict.bmk',
        col_names = c('file', 'predict_time', 'predict_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = c('species', 'k', 'target_fpr1', "target_fpr2"),
        regex = 'results/chopstitch_review/([a-z_]+).k_(\\d+).fpr1_([0-9]+.[0-9]+).fpr2_([0-9]+.[0-9]+).find.exons',
        convert = TRUE
    ) %>%
    mutate(
        predict_memory = predict_memory / 1024
    )

chopstitch_predict <- bind_rows(
    chopstitch_predict_initial, chopstitch_predict_review
)
rm(chopstitch_predict_initial, chopstitch_predict_review)


chopstitch_bwa_initial <- bwa_initial %>%
    filter(grepl(pattern = 'chopstitch', x = file)) %>%
    extract(
        col = file,
        into = c('species', 'k', 'fpr'),
        regex = 'results/bwa/([a-z_]+).k(\\d+).fpr(\\d+\\.\\d+).chopstitch.stats.tsv',
        convert = TRUE
    ) %>%
    mutate(target_fpr1 = fpr, target_fpr2 = fpr) %>%
    select(-fpr)

chopstitch_bwa_review <- bwa_review %>%
    filter(grepl(pattern = 'chopstitch', x = file)) %>%
    extract(
        col = file,
        into = c('species', 'k', 'target_fpr1', 'target_fpr2'),
        regex = 'results/bwa_review/([a-z_]+).k_(\\d+).fpr1_(\\d+\\.\\d+).fpr2_(\\d+\\.\\d+).chopstitch.stats.tsv',
        convert = TRUE
    )

chopstitch_bwa <- bind_rows(chopstitch_bwa_initial, chopstitch_bwa_review)
rm(chopstitch_bwa_initial, chopstitch_bwa_review)

chopstitch_pr_initial <- pr_initial %>%
    filter(grepl('chopstitch', file)) %>%
    extract(
        col = file,
        into = c('species', 'k', 'fpr', 'min_similarity'),
        regex = 'results/pr/([a-z_]+).k(\\d+).fpr(\\d+\\.\\d+).chopstitch.(\\d+.\\d+).tsv',
        convert = TRUE
    ) %>%
    mutate(target_fpr1 = fpr, target_fpr2 = fpr) %>%
    select(-fpr)

chopstitch_pr_review <- pr_review %>%
    filter(grepl('chopstitch', file)) %>%
    extract(
        col = file,
        into = c('species', 'k', 'target_fpr1', "target_fpr2", 'min_similarity'),
        regex = 'results/pr_review/([a-z_]+).k_(\\d+).fpr1_(\\d+\\.\\d+).fpr2_(\\d+\\.\\d+).chopstitch.(\\d+.\\d+).tsv',
        convert = TRUE
    )

chopstitch_pr <- bind_rows(chopstitch_pr_initial, chopstitch_pr_review)
rm(chopstitch_pr_initial, chopstitch_pr_review)


chopstitch_master <- chopstitch_build %>%
    left_join(chopstitch_predict) %>%
    left_join(chopstitch_pr) %>%
    left_join(chopstitch_bwa)
rm(chopstitch_build, chopstitch_predict, chopstitch_pr, chopstitch_bwa)


# gmap build
gmap_build_master <- read_tsv(
        file = 'gmap.build.bmk',
        col_names = c('file', 'build_time', 'build_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file,
        into = 'species',
        regex = 'results/gmap/(\\w+).build'
    ) %>%
    mutate(build_memory = build_memory / 1024)

# gmap predict
gmap_predict_master <- read_tsv(
        file = 'gmap.predict.bmk',
        col_names = c('file', 'predict_time', 'predict_memory'),
        col_types = "ccd"
    ) %>%
    extract(
        col = file, into = 'species',
        regex = 'results/gmap/(\\w+).gff3'
    ) %>%
    mutate(predict_memory = predict_memory / 1024)

gmap_bwa <- bwa_initial %>%
    filter(grepl(pattern = 'gmap', x = file)) %>%
    extract(
        col = file,
        into = 'species',
        regex = 'results/bwa/([a-z_]+).gmap.stats.tsv',
        convert = TRUE
    )

gmap_pr <- pr_initial %>%
    filter(grepl('gmap', file)) %>%
    extract(
        col = file,
        into = c('species', 'min_similarity'),
        regex = 'results/pr/(\\w+).gmap.(\\d+\\.\\d+).tsv',
        convert = TRUE
    )



gmap_master <- left_join(gmap_build_master, gmap_predict_master) %>%
    left_join(gmap_pr) %>%
    left_join(gmap_bwa)
rm(gmap_build_master, gmap_predict_master, gmap_pr, gmap_bwa)

rm(bwa_initial, bwa_review, pr_initial, pr_review)

exfi_master['method'] <- 'exfi'
chopstitch_master['method'] <- 'chopstitch'
gmap_master['method'] <- 'gmap'

master_table <- bind_rows(exfi_master, chopstitch_master, gmap_master) %>%
    extract(  # Convert HMS time
        col = build_time,
        into = "build_hms",
        regex = "(\\d+:\\d+:\\d+)",
        convert = TRUE,
        remove = FALSE
    ) %>%
    extract(
        col = build_time,
        into = "build_days",
        regex = "(\\d+) \\w+",
        convert = TRUE,
        remove = TRUE
    ) %>%
    replace_na(list(build_days = 0)) %>%
    mutate(
        build_hms = lubridate::hms(build_hms),
        aux = "H 0M 0S",
        build_hms_days_to_hours = build_days * 24
    ) %>%
    unite(
        "build_hms_days",
        c("build_hms_days_to_hours", "aux"),
        sep = ""
    ) %>%
    mutate(
        build_time = build_hms + hms(build_hms_days),
        predict_time = hms(predict_time),
        total_time = 
            (period_to_seconds(build_time) + period_to_seconds(predict_time)) %>%
            seconds_to_period()
    ) %>%
    rowwise() %>%
    mutate(
        max_memory = max(build_memory, predict_memory)
    ) %>%
    select(-build_days, -build_hms, -build_hms_days) %>%  # Finish converting times
    left_join(species_df) %>% # Add species names, abbreviatures, etc.
    select(species, method, everything())

rm(exfi_master, gmap_master, chopstitch_master, species_df)

write_tsv(x = master_table, path = 'master.tsv')
    