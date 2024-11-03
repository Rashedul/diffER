# cll run codes
cd /projects/blueprint_CLL_dart2/jsteif/
less /projects/blueprint_CLL_dart2/jsteif/differential_results/tmp_bash_sbatch_H3K27me3_cll_mutated_cll_unmutated_0.03.sh

# compare number of bins 
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
