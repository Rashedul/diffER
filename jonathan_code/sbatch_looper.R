.libPaths(c("/home/jsteif/R_temp_masked/x86_64-centos7-linux-gnu-library/3.5_V2","/home/jsteif/R_temp_masked/x86_64-pc-linux-gnu-library/3.2_backup", "/projects/epigenomics3/epigenomics3_results/users/esu/python3.7/lib/R/library"))
library(data.table)
library(getopt)
library(magrittr)
library(purrr)
library(stringr)
library(parallel)
library(R.utils)
library(dplyr)
library(tidyr)
library(tibble)
library(rlang)

source("/projects/epigenomics3/epigenomics3_results/users/jsteif/scripts/functions_to_source/helper_functions.R")
p_value_threshold <- 0.05
chromosome <- "chr22"
mark = "H3K27ac"

marks <- c("H3K27ac", "H3K27me3", "H3K9me3", "H3K36me3", "H3K4me3", "H3K4me1")

input_directory <- "/projects/blueprint_CLL_dart2/jsteif/symlinked_final_tables"
out_directory <- "/projects/blueprint_CLL_dart2/jsteif/differential_results"

final_tables_available_int <- list.files(input_directory)[list.files(input_directory) %>% str_detect("chr5.gz")] %>% lapply(function(x){str_split(x,"_")[[1]][3]}) %>% unlist

meta_df <- read.table("/projects/blueprint_CLL_dart2/jsteif/metadata/metadata_table", stringsAsFactors = FALSE, header = T, sep=",") %>% 
  select(-library_name) %>% distinct %>%  mutate(one = 1) %>% spread(mark, one)

meta_df[is.na(meta_df)] <- 0
meta_df <- meta_df %>% filter(sample_id_int %in% final_tables_available_int) %>% filter(peak_caller != "MACS2NoInput")

sample_type_combinations <- meta_df$sample_type %>% unique %>% combn(2)

sbatch_file_generator <- function(the_out_directory, the_sample_type1, the_sample_type2,
                                  the_chromosome, the_p_value_threshold, the_mark){
  system(str_c("bash /projects/blueprint_CLL_dart2/jsteif/scripts/sbatch_file_generator_two_group_comparision.sh -o ", 
               the_out_directory, " -t ", the_sample_type1, " -T ", the_sample_type2, " -c ", the_chromosome, 
               " -p ", the_p_value_threshold, " -m ", the_mark)) 
}

sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = "CLL",
                      the_sample_type2 = "NBC", the_chromosome = "chr1", the_p_value_threshold = 0.05, 
                      the_mark = "H3K27ac")

chromosomes <- str_c("chr", seq(1,22,1))

for (i in chromosomes){
  mclapply(seq(1,ncol(sample_type_combinations), 1), function(x){
    sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = sample_type_combinations[1,x],
                          the_sample_type2 = sample_type_combinations[2,x], the_chromosome = i, the_p_value_threshold = 0.03, 
                          the_mark = "H3K27me3")
  }, mc.cores = 10)
}


mclapply(seq(1,ncol(sample_type_combinations), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = sample_type_combinations[1,x],
                        the_sample_type2 = sample_type_combinations[2,x], the_chromosome = "chr22", the_p_value_threshold = 0.03, 
                        the_mark = "H3K27ac")
}, mc.cores = 10)


sbatch_file_generator_all_chr <- function(the_out_directory, the_sample_type1, the_sample_type2,
                                  the_p_value_threshold, the_mark){
  system(str_c("bash /projects/blueprint_CLL_dart2/jsteif/scripts/sbatch_file_generator_two_group_comparision_all_chromosomes.sh -o ", 
               the_out_directory, " -t ", the_sample_type1, " -T ", the_sample_type2, " -c ", 
               " -p ", the_p_value_threshold, " -m ", the_mark)) 
}

mclapply(seq(1,ncol(sample_type_combinations), 1), function(x){
  sbatch_file_generator_all_chr(the_out_directory = out_directory, the_sample_type1 = sample_type_combinations[1,x],
                        the_sample_type2 = sample_type_combinations[2,x], the_p_value_threshold = 0.03, 
                        the_mark = "H3K27ac")
}, mc.cores = 10)


