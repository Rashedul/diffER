.libPaths(c("/home/jsteif/R_temp_masked/x86_64-centos7-linux-gnu-library/3.5_V2","/home/jsteif/R_temp_masked/x86_64-pc-linux-gnu-library/3.2_backup", "/projects/epigenomics3/epigenomics3_results/users/esu/python3.7/lib/R/library"))
library(parallel)
library(getopt)
library(stringr)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(data.table)
library(magrittr)
library(R.utils)
library(rlang)

spec = matrix(c(
  "sample_type1", "t", 1, "character", "the first sample type",
  "sample_type2", "T", 1, "character", "the first sample type",
  "input_tables_directory", "i", 1, "character", "path to input final tables directory",
  "out_directory",   "o",    1, "character", "path to output directory where tables will be written",
  "chromosome", "c", 1, "character", "the chromosome (e.g chr10)",
  "p_value_threshold", "p", 1, "numeric", "the p-value to filter final table",
  "mark", "m", 1, "character", "the mark to be compared",
  "help",         "h",    0, "logical",   "print usage"
), byrow=TRUE, ncol=5);

opt = getopt(spec)

# if the user specifies --help or -h in the command line, print the details from above and exit (won't run your script)
if (!is.null(opt$help)) {
  cat(getopt(spec, usage=T))
  q(status=1)
}


meta_df <- read.table("/projects/blueprint_CLL_dart2/jsteif/metadata/metadata_table_corrected", stringsAsFactors = FALSE, header = T, sep=",") %>% 
  select(-library_name) %>% distinct %>%  mutate(one = 1) %>% spread(mark, one)

meta_df[is.na(meta_df)] <- 0

final_tables_available_int <- list.files(opt$input_tables_directory)[list.files(opt$input_tables_directory) %>% str_detect("chr5.gz")] %>% lapply(function(x){str_split(x,"_")[[1]][3]}) %>% unlist
meta_df <- meta_df %>% filter(sample_id_int %in% final_tables_available_int) %>% filter(peak_caller != "MACS2NoInput")

df_reader <- function(table_directory, sample_id, the_chromosome) {
  fread(str_c(table_directory, "/final_table_", sample_id, "_", the_chromosome, ".gz"), header=T, 
        col.names = c("start_pos", "h3k27ac_rpkm","h3k27ac_peak","h3k27me3_rpkm","h3k27me3_peak",
                      "h3k9me3_rpkm","h3k9me3_peak","h3k36me3_rpkm","h3k36me3_peak","h3k4me3_rpkm","h3k4me3_peak",
                      "h3k4me1_rpkm","h3k4me1_peak","ID")) %>% mutate(chromosome = the_chromosome) %>% 
    mutate(ID = as.character(ID))
}

group_peak_count_table_maker_chr <- function(the_meta_df, the_chromosome, the_mark){
  group_size <- the_meta_df$sample_id_int %>% unique %>% length
  group_df <- the_meta_df$sample_id_int %>% mclapply(function(x){
  df_reader(opt$input_tables_directory, x, the_chromosome)}, mc.cores=as.integer(group_size)) %>%
  do.call("rbind", .)
  group_df <- data.table(group_df) %>% distinct
  a_single_start_position <- group_df[1, "start_pos"] %>% unlist %>% unname
  group_size_df <- group_df[start_pos == a_single_start_position, list(h3k27ac_n = sum(!is.na(h3k27ac_peak)),h3k27me3_n = sum(!is.na(h3k27me3_peak)),
                                                                       h3k9me3_n = sum(!is.na(h3k9me3_peak)),h3k36me3_n = sum(!is.na(h3k36me3_peak)),
                                                                       h3k4me1_n = sum(!is.na(h3k4me1_peak)),h3k4me3_n = sum(!is.na(h3k4me3_peak))),
                            by=list(chromosome, start_pos)]
  
  
  group_df <- group_df[,list(h3k27ac_x=sum(h3k27ac_peak, na.rm = T), h3k27ac_n = group_size_df$h3k27ac_n, h3k27ac_rpkm=mean(h3k27ac_rpkm, na.rm = T), 
                             h3k27me3_x=sum(h3k27me3_peak, na.rm = T), h3k27me3_n = group_size_df$h3k27me3_n, h3k27me3_rpkm=mean(h3k27me3_rpkm, na.rm = T),
                             h3k9me3_x=sum(h3k9me3_peak, na.rm = T ), h3k9me3_n = group_size_df$h3k9me3_n, h3k9me3_rpkm=mean(h3k9me3_rpkm, na.rm = T),
                             h3k36me3_x=sum(h3k36me3_peak, na.rm = T), h3k36me3_n = group_size_df$h3k36me3_n, h3k36me3_rpkm=mean(h3k36me3_rpkm, na.rm = T),  
                             h3k4me1_x=sum(h3k4me1_peak, na.rm = T), h3k4me1_n = group_size_df$h3k4me1_n, h3k4me3_rpkm=mean(h3k4me3_rpkm, na.rm = T),
                             h3k4me3_x=sum(h3k4me3_peak, na.rm = T), h3k4me3_n = group_size_df$h3k4me3_n, h3k4me1_rpkm=mean(h3k4me1_rpkm, na.rm = T)),
                       by=list(chromosome, start_pos)]
  
  the_mark_x = quo(str_c(str_to_lower(the_mark), "_x"))
  the_mark_n = quo(str_c(str_to_lower(the_mark), "_n"))
  the_mark_rpkm = quo(str_c(str_to_lower(the_mark), "_rpkm"))
  
  return(group_df %>% select(chromosome, start_pos, UQ(the_mark_x), UQ(the_mark_n), UQ(the_mark_rpkm)))
}


p_value_extractor <- function(x){
  input_matrix <- matrix(x, nrow=2)
  input_matrix[2,] <- input_matrix[2,] - input_matrix[1,]
  return(fisher.test(input_matrix)$p.value)
}


fisher_tester <- function(prob_table_1, prob_table_2, the_mark, group1, group2){
  the_mark_xA = quo(str_c(str_to_lower(the_mark), "_x.", group1))
  the_mark_nA = quo(str_c(str_to_lower(the_mark), "_n.", group1))
  the_mark_rpkmA = quo(str_c(str_to_lower(the_mark), "_rpkm.", group1))
  the_mark_xB = quo(str_c(str_to_lower(the_mark), "_x.", group2))
  the_mark_nB = quo(str_c(str_to_lower(the_mark), "_n.", group2))
  the_mark_rpkmB = quo(str_c(str_to_lower(the_mark), "_rpkm.", group2))
  chromosome_A = quo(str_c("chromosome.", group1))
  merged_df <- merge(prob_table_1,prob_table_2, by=c("start_pos"), suffixes = c(str_c(".", group1),str_c(".", group2))) %>% 
	       select(chromosome = UQ(chromosome_A), start_pos,
                      UQ(the_mark_xA), UQ(the_mark_nA), UQ(the_mark_rpkmA) ,UQ(the_mark_xB), UQ(the_mark_nB), UQ(the_mark_rpkmB)) 
  
  input_combinations_table <- merged_df %>% select(UQ(the_mark_xA), UQ(the_mark_nA), UQ(the_mark_xB), UQ(the_mark_nB)) %>% na.omit %>% as.data.frame() %>% 
    distinct()
  output_table <- cbind(input_combinations_table, p_value = apply(input_combinations_table, 1,p_value_extractor)) 
  output_table[,str_c("prop_", group1)] <- output_table[,1] / output_table[,2]
  output_table[,str_c("prop_", group2)] <- output_table[,3] / output_table[,4]
  output_table[,str_c("prop", "A")] <- output_table[,1] / output_table[,2]
  output_table[,str_c("prop", "B")] <- output_table[,3] / output_table[,4]
  
  
  output_table <- output_table %>% mutate(enriched_group = ifelse( propA > propB, group1, 
                                                           ifelse( propB > propA, group2, "Equal" )),
                                          enriched_over_group = ifelse( propA < propB, group1, 
                                                                   ifelse( propB < propA, group2, "Equal" ))) %>% 
    select(-propA, -propB)
  
  final_df <- merged_df %>% full_join(output_table)
  return(final_df)
}

a <- group_peak_count_table_maker_chr(meta_df %>% filter(ig_status == opt$sample_type1), opt$chromosome, opt$mark)
b <- group_peak_count_table_maker_chr(meta_df %>% filter(ig_status == opt$sample_type2), opt$chromosome, opt$mark)
fisher_tester(a,b,opt$mark,opt$sample_type1, opt$sample_type2) %>%
filter(p_value < opt$p_value_threshold) %>% 
  write.table( str_c(opt$out_directory, "/", opt$mark, "_" , opt$sample_type1, "_", opt$sample_type2, "_", opt$chromosome, "_", opt$p_value_threshold ) , quote = F, row.names = F, sep=",")

