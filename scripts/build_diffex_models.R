library(tidyverse)
library(DESeq2)
library(tximport)

# read in and format metadata ---------------------------------------------

metadata <- read_tsv(snakemake@input[['metadata']], show_col_types = F) %>%
  filter(tissue != "cell_line") %>% # filter out samples that will make it so model is not full rank
  filter(!study_title %in% "Amblyomma americanum RNA-seq following E. coli treatment") %>% # rm e. coli study, too batched
  filter(!library_name %in% c("Um", "Uf", "Im", "If")) %>% # filter e chaf exposed samples
  select(library_name, experiment_title, study_title, sex, host_meal, tissue, blood_meal_hour, blood_meal_hour_range, total_spots) %>% # keep columns of interest
  # combine samples that had multiple SRR* accessions by library_name
  group_by(library_name) %>%
  mutate(total_spots = sum(total_spots)) %>%
  # define a variable for range of time sample was taken during blood meal hour
  mutate(blood_meal_hour_range = factor(blood_meal_hour_range, levels = c("0", "12_48", "72_144", "168_264", "none"))) %>%
  # keep only distinct rows
  distinct()

# create a vector of filenames for quant.sf files
# solve by parsing library name so that there's a 1:1 match between filename and the file that's read in for that sample
metadata <- metadata %>%
  mutate(filepath = paste0("outputs/quantification/salmon/", library_name, "_quant/quant.sf"))

# read in tx2gene file ----------------------------------------------------

# read in the tx2gene file created by assigning transcripts gene names they overlapped with when mapped to to the reference genome
# restore missing transcripts that don't have gene assignments and that didn't map to the genome
tx2gene <- read_tsv(snakemake@input[['tx2gene']],
                    skip = 1, col_names = c("tx", "gene"), show_col_types = FALSE) %>%
  filter(grepl(pattern = "evm", x = gene))

# build diffex model from sex_tissue_blood_meal_hour ----------------------

# create a new variable that combines most variables of interest
metadata_sex_tissue_blood_meal_hour <- metadata %>% 
  mutate(sex_tissue_blood_meal_hour = paste0(sex, "_x_", tissue, "_x_", blood_meal_hour_range))

# tally how many replicates there are for each observed combo of sex_tissue_blood_meal_hour
metadata_sex_tissue_blood_meal_hour_tally <- metadata_sex_tissue_blood_meal_hour %>% 
  group_by(sex_tissue_blood_meal_hour) %>% 
  tally() %>% arrange(desc(n)) %>% 
  filter(n > 1)

# filter to samples that have > 2 in the above table
metadata_sex_tissue_blood_meal_hour <- metadata_sex_tissue_blood_meal_hour %>% 
  filter(sex_tissue_blood_meal_hour %in% metadata_sex_tissue_blood_meal_hour_tally$sex_tissue_blood_meal_hour)

# read in salmon results and only keep counts that mapped to transcripts with gene assignments
salmon_counts_sex_tissue_blood_meal_hour <- tximport(files = metadata_sex_tissue_blood_meal_hour$filepath,
                                                     type = "salmon", txOut = FALSE, tx2gene = tx2gene)

dds_sex_tissue_blood_meal_hour <- DESeqDataSetFromTximport(salmon_counts_sex_tissue_blood_meal_hour,
                                                           colData = metadata_sex_tissue_blood_meal_hour,
                                                           design = ~ sex_tissue_blood_meal_hour)

ds_sex_tissue_blood_meal_hour <- DESeq(dds_sex_tissue_blood_meal_hour, test="Wald")

# write out differential expression models for use in shiny app
saveRDS(dds_sex_tissue_blood_meal_hour, snakemake@output[['dds_stb']])
saveRDS(ds_sex_tissue_blood_meal_hour, snakemake@output[['ds_stb']])

# build diffex models from sex_tissue -------------------------------------

# create a new metadata table with sex_tissue as a column
# every combo of those variables that exists has replicates so no samples need to be filtered
metadata_sex_tissue <- metadata %>% 
  mutate(sex_tissue = paste0(sex, "_x_", tissue))

# read in salmon results and only keep counts that mapped to transcripts with gene assignments
salmon_counts_sex_tissue <- tximport(files = metadata_sex_tissue$filepath,
                                     type = "salmon", txOut = FALSE, tx2gene = tx2gene)

dds_sex_tissue <- DESeqDataSetFromTximport(salmon_counts_sex_tissue,
                                           colData = metadata_sex_tissue,
                                           design = ~ sex_tissue)

ds_sex_tissue <- DESeq(dds_sex_tissue, test="Wald")

# write out differential expression models for use in shiny app
saveRDS(dds_sex_tissue, snakemake@output[['dds_st']])
saveRDS(ds_sex_tissue, snakemake@output[['ds_st']])
