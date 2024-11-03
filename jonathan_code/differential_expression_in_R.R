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
p_value_threshold <- 0.005
chromosome <- "chr22"
sample_type1 <- "GCBC"
sample_type2 <- "CLL"
mark = "H3K27ac"

input_directory <- "/projects/blueprint_CLL_dart2/jsteif/symlinked_final_tables"

final_tables_avaialble <- list.files(input_directory)[list.files(input_directory) %>% str_detect("chr5.gz")] %>% lapply(function(x){str_split(x,"_")[[1]][3]}) %>% unlist

meta_df %>% filter(sample_id_int %in% setdiff(meta_df$sample_id_int, final_tables_avaialble))



out_directory <- "/projects/blueprint_CLL_dart2/jsteif/differential_results"

meta_df <- read.table("/projects/blueprint_CLL_dart2/jsteif/metadata/metadata_table", stringsAsFactors = FALSE, header = T, sep=",") %>% 
  select(-library_name) %>% distinct %>%  mutate(one = 1) %>% spread(mark, one)

meta_df[is.na(meta_df)] <- 0

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
    df_reader(input_directory, x, the_chromosome)}, mc.cores=as.integer(group_size)) %>%
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

all_group_1 <- group_peak_count_table_maker_chr(meta_df %>% filter(sample_type == sample_type1), chromosome, mark)
all_group_2 <- group_peak_count_table_maker_chr(meta_df %>% filter(sample_type == sample_type2), chromosome, mark)



all_group_reader_chr <- function(input_tables_directory, the_meta_df, the_chromosome ){
  group_size <- the_meta_df$sample_id_int %>% unique %>% length
  print(group_size)
  group_df <- the_meta_df$sample_id_int %>% mclapply(function(x){
    df_reader(input_tables_directory, x, the_chromosome)}, mc.cores=as.integer(group_size)) %>%
    do.call("rbind", .)
  group_df <- data.table(group_df %>% distinct)
  return(group_df)
}

group_peak_count_table_maker_chr <- function(all_group_df, the_meta_df, the_chromosome, the_mark){
#  sub_meta_df <- the_meta_df %>% filter(get(the_mark) == 1) %>% sample_n(subsample_size)
#  group_df <- all_group_df[ID %in% sub_meta_df$sample_id_int,]
  # Needs one row to calculate group size
  a_single_start_position <- all_group_df[1, "start_pos"] %>% unlist %>% unname
  group_size_df <- all_group_df[start_pos == a_single_start_position, list(h3k27ac_n = sum(!is.na(h3k27ac_peak)),h3k27me3_n = sum(!is.na(h3k27me3_peak)),
                                                                       h3k9me3_n = sum(!is.na(h3k9me3_peak)),h3k36me3_n = sum(!is.na(h3k36me3_peak)),
                                                                       h3k4me1_n = sum(!is.na(h3k4me1_peak)),h3k4me3_n = sum(!is.na(h3k4me3_peak))),
                            by=list(chromosome, start_pos)]
  
  
  group_df <- all_group_df[,list(h3k27ac_x=sum(h3k27ac_peak, na.rm = T), h3k27ac_n = group_size_df$h3k27ac_n, h3k27ac_rpkm=mean(h3k27ac_rpkm, na.rm = T),  
                             h3k27me3_x=sum(h3k27me3_peak, na.rm = T), h3k27me3_n = group_size_df$h3k27me3_n, h3k27me3_rpkm=mean(h3k27me3_rpkm, na.rm = T),
                             h3k9me3_x=sum(h3k9me3_peak, na.rm = T ), h3k9me3_n = group_size_df$h3k9me3_n, h3k9me3_rpkm=mean(h3k9me3_rpkm, na.rm = T),
                             h3k36me3_x=sum(h3k36me3_peak, na.rm = T), h3k36me3_n = group_size_df$h3k36me3_n, h3k36me3_rpkm=mean(h3k36me3_rpkm, na.rm = T), 
                             h3k4me1_x=sum(h3k4me1_peak, na.rm = T), h3k4me1_n = group_size_df$h3k4me1_n, h3k4me3_rpkm=mean(h3k4me3_rpkm, na.rm = T),
                             h3k4me3_x=sum(h3k4me3_peak, na.rm = T), h3k4me3_n = group_size_df$h3k4me3_n, h3k4me1_rpkm=mean(h3k4me1_rpkm, na.rm = T)),
                       by=list(chromosome, start_pos)]
  
  the_mark_x = quo(str_c(str_to_lower(the_mark), "_x"))
  the_mark_rpkm = quo(str_c(str_to_lower(the_mark), "_rpkm"))
  the_mark_n = quo(str_c(str_to_lower(the_mark), "_n"))
  
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
  the_mark_xB = quo(str_c(str_to_lower(the_mark), "_x.", group2))
  the_mark_nB = quo(str_c(str_to_lower(the_mark), "_n.", group2))
  chromosome_A = quo(str_c("chromosome.", group1))
  merged_df <- merge(prob_table_1,prob_table_2, by=c("start_pos"), suffixes = c(str_c(".", group1),str_c(".", group2))) %>% select(chromosome = UQ(chromosome_A), start_pos,
                                                                                                                                   UQ(the_mark_xA), UQ(the_mark_nA), UQ(the_mark_xB), UQ(the_mark_nB)) 
  
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

all_group_1 <- all_group_reader_chr(input_directory, meta_df %>% filter(sample_type == sample_type1), chromosome)
all_group_2 <- all_group_reader_chr(input_directory, meta_df %>% filter(sample_type == sample_type2), chromosome)

fisher_tester(all_group_1, all_group_2)

group_1_sub_list <- mclapply(seq(1,number_of_comparisions, 1), function(x){
  group_peak_count_table_maker_chr(all_group_1, meta_df %>% filter(sample_type == sample_type1), chromosome, mark, 3)
}, mc.cores = 5)

group_2_sub_list <- mclapply(seq(1,number_of_comparisions, 1), function(x){
  group_peak_count_table_maker_chr(all_group_2, meta_df %>% filter(sample_type == sample_type2), chromosome, mark, 3)
}, mc.cores = 5)

fisher_list <- map2(group_1_sub_list, group_2_sub_list, function(a,b){
  fisher_tester(a,b, mark, sample_type1, sample_type2)
})

fisher_list %>% do.call("rbind", .) %>% filter(p_value < p_value_threshold) %>% group_by(chromosome, start_pos, enriched_group, enriched_over_group) %>% 
  dplyr::summarize(Comparisions_enriched = n()) %>% mutate(number_of_comparisions = number_of_comparisions) %>% head


out_directory <- "/projects/epigenomics_assembly/jsteif/results/fisher_simulations"


sbatch_file_generator <- function(the_out_directory, the_sample_type1, the_sample_type2, the_subsample_size, min_subample_size,
                                  the_chromosome, the_p_value_threshold, the_number_of_comparisions, the_mark){
  system(str_c("bash /projects/epigenomics3/epigenomics3_results/users/jsteif/scripts/TERM3/Enhancer_detection_sample_size/tmp_sbatch_file_generator_two_group_comparision_simulation_v3.sh -o ", 
               the_out_directory, " -t ", the_sample_type1, " -T ", the_sample_type2, " -s " , the_subsample_size, " -c ", the_chromosome, 
               " -n ", the_number_of_comparisions, " -p ", the_p_value_threshold, " -m ", the_mark, " -l ", min_subample_size)) 
}

sbatch_file_generator(the_out_directory = "/projects/epigenomics_assembly/jsteif/results/fisher_simulations", the_sample_type1 = "Dendritic",
                      the_sample_type2 = "Plasma", the_subsample_size = 5, the_chromosome = "chr18", the_p_value_threshold = 0.05, 
                      the_number_of_comparisions = 5, the_mark = "H3K27ac")

h3k27ac_16_or_more <- meta_df %>% group_by(sample_type) %>% dplyr::summarize(n = n(), H3K27ac = sum(H3K27ac)) %>% filter(H3K27ac > 15)
h3k27ac_16_or_more_comb <- h3k27ac_16_or_more$sample_type %>% combn(2)

mclapply(seq(1,ncol(h3k27ac_16_or_more_comb), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = h3k27ac_16_or_more_comb[1,x],
                        the_sample_type2 = h3k27ac_16_or_more_comb[2,x], the_subsample_size = 25, the_chromosome = "chr18", the_p_value_threshold = 0.05, 
                        the_number_of_comparisions = 5, the_mark = "H3K27ac", min_subample_size = 5)
}, mc.cores = 40)



h3k27me3_6_or_more <- meta_df %>% group_by(sample_type) %>% dplyr::summarize(n = n(), H3K27me3 = sum(H3K27me3)) %>% filter(H3K27me3 > 5)
h3k27me3_6_or_more_comb <- h3k27me3_6_or_more$sample_type %>% combn(2)

mclapply(seq(1,ncol(h3k27me3_6_or_more_comb), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = h3k27me3_6_or_more_comb[1,x],
                        the_sample_type2 = h3k27me3_6_or_more_comb[2,x], the_subsample_size = 25, the_chromosome = "chr17", the_p_value_threshold = 0.05, 
                        the_number_of_comparisions = 5, the_mark = "H3K27me3", min_subample_size = 5)
}, mc.cores = 40)


h3k4me1_6_or_more <- meta_df %>% group_by(sample_type) %>% dplyr::summarize(n = n(), H3K4me1 = sum(H3K4me1)) %>% filter(H3K4me1 > 5)
h3k4me1_6_or_more_comb <- h3k4me1_6_or_more$sample_type %>% combn(2)

mclapply(seq(1,ncol(h3k4me1_6_or_more_comb), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = h3k4me1_6_or_more_comb[1,x],
                        the_sample_type2 = h3k4me1_6_or_more_comb[2,x], the_subsample_size = 25, the_chromosome = "chr18", the_p_value_threshold = 0.05, 
                        the_number_of_comparisions = 5, the_mark = "H3K4me1", min_subample_size = 5)
}, mc.cores = 40)



h3k4me3_6_or_more <- meta_df %>% group_by(sample_type) %>% dplyr::summarize(n = n(), H3K4me3 = sum(H3K4me3)) %>% filter(H3K4me3 > 5)
h3k4me3_6_or_more_comb <- h3k4me3_6_or_more$sample_type %>% combn(2)

mclapply(seq(1,ncol(h3k4me3_6_or_more_comb), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = h3k4me3_6_or_more_comb[1,x],
                        the_sample_type2 = h3k4me3_6_or_more_comb[2,x], the_subsample_size = 5, the_chromosome = "chr17", the_p_value_threshold = 0.05, 
                        the_number_of_comparisions = 5, the_mark = "H3K4me3")
}, mc.cores = 40)

h3k36me3_6_or_more <- meta_df %>% group_by(sample_type) %>% dplyr::summarize(n = n(), H3K36me3 = sum(H3K36me3)) %>% filter(H3K36me3 > 5)
h3k36me3_6_or_more_comb <- h3k36me3_6_or_more$sample_type %>% combn(2)


mclapply(seq(1,ncol(h3k36me3_6_or_more_comb), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = h3k36me3_6_or_more_comb[1,x],
                        the_sample_type2 = h3k36me3_6_or_more_comb[2,x], the_subsample_size = 5, the_chromosome = "chr17", the_p_value_threshold = 0.05, 
                        the_number_of_comparisions = 5, the_mark = "H3K36me3")
}, mc.cores = 40)

h3k9me3_6_or_more <- meta_df %>% group_by(sample_type) %>% dplyr::summarize(n = n(), H3K9me3 = sum(H3K9me3)) %>% filter(H3K9me3 > 5)
h3k9me3_6_or_more_comb <- h3k9me3_6_or_more$sample_type %>% combn(2)


mclapply(seq(1,ncol(h3k9me3_6_or_more_comb), 1), function(x){
  sbatch_file_generator(the_out_directory = out_directory, the_sample_type1 = h3k9me3_6_or_more_comb[1,x],
                        the_sample_type2 = h3k9me3_6_or_more_comb[2,x], the_subsample_size = 5, the_chromosome = "chr17", the_p_value_threshold = 0.05, 
                        the_number_of_comparisions = 5, the_mark = "H3K9me3")
}, mc.cores = 40)



fisher_comparisions <-  list.files(out_directory)[list.files(out_directory) %>% str_detect("fisher")]

fisher_simulation_reader <- function(simulation_dir, file){
  fread(str_c(prob_table_dir, "/prob_table_", sample_type, "_", chromosome), header=T)
}

fisher_simulation_list <- lapply(fisher_comparisions, function(x){
  fread(str_c(out_directory, "/", x))
})

#saveRDS(fisher_simulation_list, str_c(out_directory, "/H3K27ac_fisher_simulation_list_chr18.RDS"))

marks <- c("H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3", "H3K4me1")

fisher_simulation_df_minus_h3k27ac <- mclapply(marks, function(m){
  mark_list <- lapply(fisher_comparisions[fisher_comparisions %>% str_detect(m)], function(x){
    fread(str_c(out_directory, "/", x)) 
  })
  return(mark_list %>% do.call("rbind", .))
}, mc.cores = 5) %>% do.call("rbind", .)

#fisher_simulation_df_h3k27ac <- readRDS(str_c(out_directory, "/H3K27ac_fisher_simulation_list_chr18.RDS")) %>% do.call("rbind", .) %>% mutate(mark = "H3K27ac")
fisher_simulation_df_h3k27ac_chr18 <- readRDS(str_c(out_directory, "/H3K27ac_fisher_simulation_list_chr18.RDS")) %>% do.call("rbind", .) %>% 
  mutate(mark = "H3K27ac")
chr17_subsample_5 <- fisher_comparisions[fisher_comparisions %>% str_detect("chr17$")]
fisher_simulation_df_h3k27ac_chr17 <- chr17_subsample_5[chr17_subsample_5 %>% str_detect("H3K27ac")] %>% mclapply(function(x){
  fread(str_c(out_directory, "/", x)) 
}, mc.cores = 50) %>% do.call("rbind", .)

plot_df <- fisher_simulation_list %>% do.call("rbind", .) %>% group_by(enriched_group, enriched_over_group, Comparisions_enriched, mark ) %>% 
  dplyr::summarize(number_of_bins = n())

plot_df <- fisher_simulation_all_marks_df %>% group_by(enriched_group, enriched_over_group, Comparisions_enriched, mark ) %>% 
  dplyr::summarize(number_of_bins = n())

blood_sample_types <- c("AML", "B_cell", "CLL", "Dendritic", "Macrophage", "Monocyte", "Neutrophil", "Plasma", "T_cell")
non_blood_sample_types <- c("Basal", "Brain", "Brain_Glioma", "Colon", "Epithelial", "Hepatocyte", "Luminal", "Luminal_progenitor", 
                            "Placenta", "Renal_Cancer","Muscle", "Stromal", "Thyroid")

plot_df <- fisher_simulation_list_h3k4me1 %>% do.call("rbind", .) %>% group_by(enriched_group, enriched_over_group, Comparisions_enriched ) %>% 
  dplyr::summarize(number_of_bins = n()) %>% ungroup %>% mutate(nc = Comparisions_enriched * number_of_bins) %>%  
  group_by(enriched_group, enriched_over_group) %>% dplyr::summarise(avg_bins= sum(nc)/5) %>% 
  mutate(enriched_group = ifelse(enriched_group == "Acute_Myeloid_Leukemia", "AML",
                                 ifelse(enriched_group == "Chronic_Lymphocytic_Leukemia", "CLL", enriched_group)),
         enriched_over_group = ifelse(enriched_over_group == "Acute_Myeloid_Leukemia", "AML",
                                      ifelse(enriched_over_group == "Chronic_Lymphocytic_Leukemia", "CLL", enriched_over_group))) %>% 
  mutate(enriched_group = factor(enriched_group, levels = c(blood_sample_types, non_blood_sample_types)) ,
         enriched_over_group = factor(enriched_over_group, levels = c(rev(non_blood_sample_types), rev(blood_sample_types))))  %>% 
  filter(enriched_group != "Epithelial") %>% filter(enriched_over_group != "Epithelial") %>% 
  mutate(facet_enriched = ifelse(enriched_group %in% blood_sample_types, "Blood", "Non-blood")) %>% 
  mutate(facet_enriched_over = ifelse(enriched_over_group %in% blood_sample_types, "Blood", "Non-blood")) %>% 
  mutate(enriched_color = ifelse(facet_enriched == "Blood", "red", "black"),
         enriched_over_color = ifelse(facet_enriched_over == "Blood", "red", "black")) %>% ungroup

# H3K27ac color 51,102,255 (blue) #3366FF
# H3K4me1 color 255,102,51 (orange)  #FF6633
# H3K4me3 color 255,0,0 (red) #FF0000
# H3K9me3 color 0,0,102 (dark blue)#000066
# H3K27me3 color 102,51,0 (brown) #663300
# H3K36me3 color 153,0,153 (purple) #990099
# H3K9ac color 153,0,153 (purple)
# Input color 0,0,0 (black)

mark_list <- c("H3K27ac", "H3K27me3", "H3K36me3", "H3K9me3", "H3K4me3", "H3K4me1" )

mark_colors <- c("#3366FF", "#663300", "#990099" ,"#000066", "#FF0000", "#FF6633")



fisher_simulation_all_marks_df <- rbind(fisher_simulation_df_h3k27ac_chr17, fisher_simulation_df_h3k27ac_chr18, fisher_simulation_df_minus_h3k27ac)

sample_type_name_replacer <- function(a_character_vector){
  str_replace_all(a_character_vector, c("Luminal_progenitor" = "Luminal progenitor", 
                                        "Brain_Glioma" = "Brain Glioma",
                                        "B_cell" = "B-Cell", "T_cell" = "T-Cell",
                                        "Muscle" = "Skeletal Muscle",
                                        "Renal_Cancer" = "Renal Cancer",
                                        "Chronic_Lymphocytic_Leukemia" = "CLL"))
}


plot_df <- fisher_simulation_all_marks_df %>% group_by(enriched_group, enriched_over_group, Comparisions_enriched, mark ) %>% 
  dplyr::summarize(number_of_bins = n()) %>% filter(mark == mark_list[3]) %>%  ungroup %>% mutate(nc = Comparisions_enriched * number_of_bins) %>%  
  group_by(enriched_group, enriched_over_group) %>% dplyr::summarise(avg_bins= sum(nc)/5) %>% 
  mutate(enriched_group = ifelse(enriched_group == "Acute_Myeloid_Leukemia", "AML",
                                 ifelse(enriched_group == "Chronic_Lymphocytic_Leukemia", "CLL", enriched_group)),
         enriched_over_group = ifelse(enriched_over_group == "Acute_Myeloid_Leukemia", "AML",
                                      ifelse(enriched_over_group == "Chronic_Lymphocytic_Leukemia", "CLL", enriched_over_group))) %>% 
  mutate(enriched_group = factor(enriched_group, levels = c(blood_sample_types, non_blood_sample_types)) ,
         enriched_over_group = factor(enriched_over_group, levels = c(rev(non_blood_sample_types), rev(blood_sample_types))))  %>% 
  filter(enriched_group != "Epithelial") %>% filter(enriched_over_group != "Epithelial") %>% 
  mutate(facet_enriched = ifelse(enriched_group %in% blood_sample_types, "Blood", "Non-blood")) %>% 
  mutate(facet_enriched_over = ifelse(enriched_over_group %in% blood_sample_types, "Blood", "Non-blood")) %>% 
  mutate(enriched_color = ifelse(facet_enriched == "Blood", "red", "black"),
         enriched_over_color = ifelse(facet_enriched_over == "Blood", "red", "black")) %>% ungroup %>% 
  mutate(enriched_group = sample_type_name_replacer(enriched_group), enriched_over_group = sample_type_name_replacer(enriched_over_group))

plot_df_2 <- plot_df %>% 
  filter(enriched_group != "Renal Cancer") %>% filter(enriched_over_group != "Renal Cancer") %>% 
  # filter(enriched_group != "Hepatocyte") %>% filter(enriched_over_group != "Hepatocyte") %>% 
  filter(enriched_group != "Stromal") %>% filter(enriched_over_group != "Stromal") %>% 
  mutate(enriched_group = factor(enriched_group, levels = c(sample_type_name_replacer(blood_sample_types), sample_type_name_replacer(non_blood_sample_types[-c(10,12)]))),
         enriched_over_group = factor(enriched_over_group, levels = c(rev(sample_type_name_replacer(non_blood_sample_types[-c(10,12)])), 
                                                                      rev(sample_type_name_replacer(blood_sample_types)))))


y_axis_col <- plot_df_2 %>% select(enriched_over_group , enriched_over_color) %>% distinct %>% arrange(enriched_over_group) %>% ungroup %>% select(enriched_over_color) %>% unlist %>% unname

x_axis_col <- plot_df_2 %>% select(enriched_group , enriched_color) %>% distinct %>% arrange(enriched_group) %>% ungroup %>% select(enriched_color) %>% unlist %>% unname

a <- ggplot(data = plot_df_2, aes(x=enriched_group, y=enriched_over_group, fill=avg_bins)) + 
  geom_tile() + scale_fill_gradient(low = "white",
                                    high = "darkgreen")  + theme_bw() + 
  # high = mark_colors[3])  + theme_bw() + 
  geom_point(data = plot_df_2, aes(shape=NA ,col = facet_enriched)) +
  labs(x = "Enriched Group", y = "Comparision Group", fill = "Average Enriched Bins", color = "",
       title = str_c( mark_list[3], " Enriched Regions in 10 Simulated Random Pairwise Comparisions (Chr 17+18)")) +
  theme(axis.text.y = element_text(color = y_axis_col), axis.text.x = element_text(color = x_axis_col))  +
  scale_color_manual(values = c("red", "black")) +  geom_point(data = plot_df_2, aes(shape=NA, col = facet_enriched)) +
  guides(colour = guide_legend(override.aes = list(size=10))) 



tile_theme <- function(a_ggplot,x_angle=0, y_angle=0){
  a_ggplot <- a_ggplot + theme(axis.text.x = element_text(size=14, angle=x_angle),
                               axis.text.y = element_text(size=14, angle=y_angle),
                               axis.title = element_text(size=14, face = "bold"),
                               legend.text = element_text(size=14),
                               legend.title = element_text(size=14, face="bold"),
                               strip.text = element_text(face="bold", size=14),
                               plot.title = element_text(face="bold", size=16, hjust=0.5),
                               plot.caption = element_text(size=14, face="bold", hjust=0))
  return(a_ggplot)
}


a %>% tile_theme(x_angle = 90)



fisher_simulation_df_h3k27ac_chr18 <- readRDS(str_c(out_directory, "/H3K27ac_fisher_simulation_list_chr18.RDS")) %>% do.call("rbind", .) %>% 
  mutate(mark = "H3K27ac", subsample_size = 5)

chr17_subsample_5 <- fisher_comparisions[fisher_comparisions %>% str_detect("chr17$")]

fisher_simulation_df_h3k27ac_chr17 <- chr17_subsample_5[chr17_subsample_5 %>% str_detect("H3K27ac")] %>% mclapply(function(x){
  fread(str_c(out_directory, "/", x)) 
}, mc.cores = 50) %>% do.call("rbind", .) %>% mutate(subsample_size = 5)

fisher_simulation_df_h3k27ac_chr17_n <- fisher_comparisions[fisher_comparisions %>% str_detect("chr17_n_")] %>% mclapply(function(x){
  fread(str_c(out_directory, "/", x)) 
}, mc.cores = 50) %>% do.call("rbind", .)

fisher_simulation_df_h3k27ac_chr17 <- rbind(fisher_simulation_df_h3k27ac_chr17, fisher_simulation_df_h3k27ac_chr17_n)
fisher_simulation_df_h3k27ac_chr17 %>% filter(enriched_group == "Thyroid" & enriched_over_group == "T_cell")

enhancer_atlas_df_reader <- function(enhanceratlas_directory = enhanceratlas_directory, chromosome = chromosome, sample_type = sample_type, fantom = F){
  if (fantom){
    df <- fread(str_c(enhanceratlas_directory, "/FANTOM5_hg38_", sample_type ,
                      "_differentially_expressed_enhancers_50bp_", chromosome, ".bed"), header=F) %>% select(-V4)
    colnames(df) <- c("chromosome","start_pos","end_pos" ,"enahncer_start", "enhancer_end")
    return(df %>% select(-end_pos)) 
  }
  else{
    df <- fread(str_c(enhanceratlas_directory, "/hglft_enhancers_", sample_type ,"_50bp","_", chromosome, ".bed"), header=F) %>% select(-V4)
    colnames(df) <- c("chromosome","start_pos","end_pos" ,"enahncer_start", "enhancer_end", "enhancer_score")
    return(df %>% select(-end_pos))
  }
}

joined_df_maker <- function(enhanceratlas_directory = enhanceratlas_directory, enhancer_atlas_sample_type = enhancer_atlas_sample_type,
                            chromosome = chromosome, fisher_simulation_df = fisher_simulation_df, fantom = fantom){
  ea_df <- enhancer_atlas_df_reader(enhanceratlas_directory, chromosome, enhancer_atlas_sample_type, fantom)
  inner_join(ea_df, fisher_simulation_df, by=c("chromosome", "start_pos")) %>% return()
}

enhanceratlas_directory <- "/projects/epigenomics_assembly/jsteif/results/sample_size_enhancer_detection"
prob_table_directory <- "/projects/epigenomics_assembly/jsteif/IHEC/prob_tables"

enhancer_atlas_sample_types <- c("B_cell", "Monocyte", "Skeletal_muscle", "GCB", "Trophoblast", "Macrophages", "Dendritic_cell")
fantom_sample_types <- c("monocyte", "skeletal_muscle", "macrophage", "dendritic", "brain", "T_cell")
prob_table_sampl_types <- c( "Monocyte", "Muscle", "Macrophage", "Dendritic", "Brain", "T_cell")

fisher_monocyte_h3k27ac_chr17_df <- fisher_simulation_df_h3k27ac_chr17 %>% filter(enriched_group == "Monocyte")
fisher_macrophage_h3k27ac_chr17_df <- fisher_simulation_df_h3k27ac_chr17 %>% filter(enriched_group == "Macrophage")
fisher_T_cell_h3k27ac_chr17_df <- fisher_simulation_df_h3k27ac_chr17 %>% filter(enriched_group == "T_cell")
fisher_B_cell_h3k27ac_chr17_df <- fisher_simulation_df_h3k27ac_chr17 %>% filter(enriched_group == "B_cell")

joined_fisher_monocyte_h3k27ac_chr17_df <-  joined_df_maker(enhanceratlas_directory = enhanceratlas_directory, enhancer_atlas_sample_type = "monocyte", chromosome = "chr17",
                                                            fisher_monocyte_h3k27ac_chr17_df , fantom=T)
joined_fisher_macrophage_h3k27ac_chr17_df <-  joined_df_maker(enhanceratlas_directory = enhanceratlas_directory, enhancer_atlas_sample_type = "macrophage", chromosome = "chr17",
                                                              fisher_macrophage_h3k27ac_chr17_df , fantom=T)
joined_fisher_T_cell_h3k27ac_chr17_df <- joined_df_maker(enhanceratlas_directory = enhanceratlas_directory, enhancer_atlas_sample_type = "T_ell", chromosome = "chr17",
                                                         fisher_macrophage_h3k27ac_chr17_df , fantom=T)


joined_fisher_B_cell_h3k27ac_chr17_df <-  joined_df_maker(enhanceratlas_directory = enhanceratlas_directory, enhancer_atlas_sample_type = "B_cell", chromosome = "chr17",
                                                          fisher_B_cell_h3k27ac_chr17_df , fantom=F)


joined_fisher_monocyte_h3k27ac_chr17_df %>% filter(Comparisions_enriched > 4) %>% group_by(subsample_size, enriched_over_group ) %>% 
  dplyr::summarise(n =n()) %>% filter(enriched_over_group == "T_cell") %>% 
  ggplot(aes(x = subsample_size, y = n)) + geom_point()

joined_fisher_monocyte_h3k27ac_chr17_df %>% filter(Comparisions_enriched > 2) %>% group_by(chromosome, subsample_size,enriched_group, enriched_over_group, enahncer_start, enhancer_end ) %>% 
  dplyr::summarise(n =n()) %>% ungroup %>% group_by(chromosome, subsample_size, enriched_group, enriched_over_group ) %>%   dplyr::summarise(n2 =n()) %>% 
  filter(enriched_over_group == "Thyroid") %>% 
  ggplot(aes(x = subsample_size, y = n2)) + geom_point()

p <- rbind(joined_fisher_monocyte_h3k27ac_chr17_df, joined_fisher_macrophage_h3k27ac_chr17_df) %>% filter(Comparisions_enriched > 1) %>% group_by(chromosome, subsample_size,enriched_group, enriched_over_group, enahncer_start, enhancer_end ) %>% 
  dplyr::summarise(n =n()) %>% ungroup %>% group_by(chromosome, subsample_size, enriched_group, enriched_over_group ) %>%   dplyr::summarise(n2 =n()) %>% 
  filter(enriched_over_group %in% c("Chronic_Lymphocytic_Leukemia", "T_cell", "B_cell")) %>% mutate(enriched_over_group = sample_type_name_replacer(enriched_over_group)) %>% 
  ggplot(aes(x = subsample_size, y = n2, color = enriched_over_group)) + geom_point() + facet_wrap(~ enriched_group) + geom_smooth() + theme_bw() +
  labs(x = "Subsample size", y = "Differentially marked FANTOM Enhancers", color = "Comparision Group") 

p %>% plot_theme()v


joined_fisher_T_cell_h3k27ac_chr17_df %>% filter(Comparisions_enriched > 1) %>% group_by(chromosome, subsample_size,enriched_group, enriched_over_group, enahncer_start, enhancer_end ) %>% 
  dplyr::summarise(n =n()) %>% ungroup %>% group_by(chromosome, subsample_size, enriched_group, enriched_over_group ) %>%   dplyr::summarise(n2 =n()) %>% 
  filter(enriched_over_group == "Monocyte") %>% 
  ggplot(aes(x = subsample_size, y = n2)) + geom_point() + theme_bw()




joined_fisher_B_cell_h3k27ac_chr17_df %>% filter(Comparisions_enriched > 3) %>% group_by(chromosome, subsample_size,enriched_group, enriched_over_group, enahncer_start, enhancer_end ) %>% 
  dplyr::summarise(n =n()) %>% ungroup %>% group_by(chromosome, subsample_size, enriched_group, enriched_over_group ) %>%   dplyr::summarise(n2 =n()) %>% 
  filter(enriched_over_group == "Neutrophil") %>% 
  ggplot(aes(x = subsample_size, y = n2)) + geom_point()


p_df <- fisher_simulation_df_h3k27ac_chr17 %>% filter(Comparisions_enriched > 3) %>%  
  group_by(enriched_group, enriched_over_group, mark, subsample_size ) %>% 
  dplyr::summarize(number_of_bins = n()) %>% filter(mark == "H3K27ac") %>%  ungroup 

mutate(nc = Comparisions_enriched * number_of_bins) 

plot_df <- plot_df %>% group_by(enriched_group, enriched_over_group, mark, subsample_size) %>% dplyr::summarise(avg_bins= sum(nc)/5) %>% 
  mutate(enriched_group = ifelse(enriched_group == "Acute_Myeloid_Leukemia", "AML",
                                 ifelse(enriched_group == "Chronic_Lymphocytic_Leukemia", "CLL", enriched_group)),
         enriched_over_group = ifelse(enriched_over_group == "Acute_Myeloid_Leukemia", "AML",
                                      ifelse(enriched_over_group == "Chronic_Lymphocytic_Leukemia", "CLL", enriched_over_group))) %>% 
  mutate(enriched_group = factor(enriched_group, levels = c(blood_sample_types, non_blood_sample_types)) ,
         enriched_over_group = factor(enriched_over_group, levels = c(rev(non_blood_sample_types), rev(blood_sample_types))))  %>% 
  filter(enriched_group != "Epithelial") %>% filter(enriched_over_group != "Epithelial") %>% 
  mutate(facet_enriched = ifelse(enriched_group %in% blood_sample_types, "Blood", "Non-blood")) %>% 
  mutate(facet_enriched_over = ifelse(enriched_over_group %in% blood_sample_types, "Blood", "Non-blood")) %>% 
  mutate(enriched_color = ifelse(facet_enriched == "Blood", "red", "black"),
         enriched_over_color = ifelse(facet_enriched_over == "Blood", "red", "black")) %>% ungroup




plot_df %>% filter(enriched_group == "Thyroid" & enriched_over_group == "T_cell")









plot_df %>% ggplot(aes(x = n))


group_peak_count_table_maker_chr_2 <- function(input_tables_directory, the_chromosome, the_meta_df, the_mark, subsample_size ){
  sub_meta_df <- the_meta_df %>% filter(get(the_mark) == 1) %>% sample_n(subsample_size)
  group_size <- sub_meta_df$sample_id_int %>% unique %>% length
  print(group_size)
  group_df <- sub_meta_df$sample_id_no_dec %>% mclapply(function(x){
    df_reader(input_tables_directory, x, the_chromosome)}, mc.cores=as.integer(group_size)) %>%
    do.call("rbind", .)
  group_df <- data.table(group_df) %>% distinct
  # Needs one row to calculate group size
  a_single_start_position <- group_df[1, "start_pos"] %>% unlist %>% unname
  group_size_df <- group_df[start_pos == a_single_start_position, list(h3k27ac_n = sum(!is.na(h3k27ac_peak)),h3k27me3_n = sum(!is.na(h3k27me3_peak)),
                                                                       h3k9me3_n = sum(!is.na(h3k9me3_peak)),h3k36me3_n = sum(!is.na(h3k36me3_peak)),
                                                                       h3k4me1_n = sum(!is.na(h3k4me1_peak)),h3k4me3_n = sum(!is.na(h3k4me3_peak))),
                            by=list(chromosome, start_pos)]
  
  
  group_df <- group_df[,list(h3k27ac_x=sum(h3k27ac_peak, na.rm = T), h3k27ac_n = group_size_df$h3k27ac_n, 
                             h3k27me3_x=sum(h3k27me3_peak, na.rm = T), h3k27me3_n = group_size_df$h3k27me3_n, 
                             h3k9me3_x=sum(h3k9me3_peak, na.rm = T ), h3k9me3_n = group_size_df$h3k9me3_n,
                             h3k36me3_x=sum(h3k36me3_peak, na.rm = T), h3k36me3_n = group_size_df$h3k36me3_n, 
                             h3k4me1_x=sum(h3k4me1_peak, na.rm = T), h3k4me1_n = group_size_df$h3k4me1_n,
                             h3k4me3_x=sum(h3k4me3_peak, na.rm = T), h3k4me3_n = group_size_df$h3k4me3_n),
                       by=list(chromosome, start_pos)]
  
  the_mark_x = quo(str_c(str_to_lower(the_mark), "_x"))
  the_mark_n = quo(str_c(str_to_lower(the_mark), "_n"))
  
  return(group_df %>% select(chromosome, start_pos, UQ(the_mark_x), UQ(the_mark_n)))
}



a <- group_peak_count_table_maker_chr_2(input_directory, "chr18", meta_df %>% filter(sample_type == sample_type1) , "H3K27ac", 5)
b <- group_peak_count_table_maker_chr_2(input_directory, "chr18", meta_df %>% filter(sample_type == sample_type2) , "H3K27ac", 5)
d <- fisher_tester(a,b, "H3K27ac", sample_type1, sample_type2)

final_tables <- list()
for (i in 1:number_of_comparisions){
  print(i)
  final_tables[[i]] <- fisher_tester(a,b, "H3K27ac", sample_type1, sample_type2)
}



