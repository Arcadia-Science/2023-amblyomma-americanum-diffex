#!/usr/bin/env Rscript

library(dplyr)
library(purrr)
library(readr)
library(tidyr)

# command line args -------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
output_tsv   <- args[1]
input_counts <- args[2:length(args)]


# read in and combine count files -----------------------------------------

raw_counts <- input_counts %>%
  set_names() %>% 
  map_dfr(function(x) read_tsv(x, col_names = c("gene", "count")), .id = "sample") %>%
  mutate(sample = gsub("outputs/counts/htseq_count/", "", sample)) %>%
  mutate(sample = gsub("\\.out", "", sample)) %>%
  pivot_wider(id_cols = gene, names_from = sample, values_from = count)

write_tsv(raw_counts, output_tsv)