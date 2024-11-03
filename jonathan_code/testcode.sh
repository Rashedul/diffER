# compare with Rashedul's version

# cll run codes
cd /projects/blueprint_CLL_dart2/jsteif/
less /projects/blueprint_CLL_dart2/jsteif/differential_results/tmp_bash_sbatch_H3K27me3_cll_mutated_cll_unmutated_0.03.sh
bash /projects/blueprint_CLL_dart2/jsteif/scripts/run_opt_two_group_comparision.sh -t cll_mutated -T cll_unmutated -p 0.03 -c chr11 -m H3K27me3 -i /projects/blueprint_CLL_dart2/jsteif/symlinked_final_tables -g ig_status -o /projects/blueprint_CLL_dart2/jsteif/differential_results

# compare number of bins (okay)
cd /projects/epigenomics_assembly/jsteif/IHEC/postgresql_tables/ALL/final_tables
g=$cll/genome_feature/CHROMSIZES_hg38_formatted.txt
chr22=$cll/genome_feature/CHROMSIZES_hg38_chr22.txt
chr10=$cll/genome_feature/CHROMSIZES_hg38_chr10.txt

# cll jonathan 
less /projects/blueprint_CLL_dart2/analysis/chipSeq/bamCoverage_and_peaks_database_prep/final_tables/final_table_1358551_chr22.gz | wc -l
# 1016371
# 2675950 (chr 10)
bedtools makewindows -g $chr22 -w 50 | wc -l
#1016370
bedtools makewindows -g $chr10 -w 50 | wc -l
# 2675949

bedtools makewindows -g $g -w 50 | wc -l
#64185939 

# compare number of DE bins before merging (does not match)
cd /projects/blueprint_CLL_dart2/jsteif/differential_results/tmp_bash
less H3K27me3_cll_mutated_cll_unmutated_chr10_0.03 | wc -l 
# 16036

less H3K27me3_cll_mutated_cll_unmutated_chr10_0.03 | grep -v "chromosome" | tr ',' "\t" | awk '{print $1 "\t" $2 "\t" $2+50}' | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 | wc -l
# 5097


################################
# fisher test
p_value_extractor <- function(x) {
  # Check if the input is a vector of length 4
  if (length(x) != 4) {
    stop("Input vector must contain exactly 4 elements.")
  }
  
  # Create a 2x2 matrix from the input vector
  input_matrix <- matrix(x, nrow = 2)
  
  # Adjust the second row
  input_matrix[2, ] <- input_matrix[2, ] - input_matrix[1, ]
  
  # Perform Fisher's exact test and return the p-value
  return(fisher.test(input_matrix)$p.value)
}


# Example vector with 4 elements
x <- c(10, 15, 5, 20)

# Extract p-value
p_value <- p_value_extractor(x)
print(p_value)
