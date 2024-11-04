# compare with Rashedul's version

#################
# cll run codes
cd /projects/blueprint_CLL_dart2/jsteif/
less /projects/blueprint_CLL_dart2/jsteif/differential_results/tmp_bash_sbatch_H3K27me3_cll_mutated_cll_unmutated_0.03.sh
bash /projects/blueprint_CLL_dart2/jsteif/scripts/run_opt_two_group_comparision.sh -t cll_mutated -T cll_unmutated -p 0.03 -c chr11 -m H3K27me3 -i /projects/blueprint_CLL_dart2/jsteif/symlinked_final_tables -g ig_status -o /projects/blueprint_CLL_dart2/jsteif/differential_results

##################
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

###################
# compare number of DE bins before merging (does not match)
cd /projects/blueprint_CLL_dart2/jsteif/differential_results/tmp_bash
less H3K27me3_cll_mutated_cll_unmutated_chr10_0.03 | wc -l 
# 16,036

# js
less H3K27me3_cll_mutated_cll_unmutated_chr10_0.03 | grep -v "chromosome" | tr ',' "\t" | awk '{print $1 "\t" $2 "\t" $2+50}' | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 | wc -l
# 5097

# js merge 100, min len 300bp
less H3K27me3_cll_mutated_cll_unmutated_chr10_0.03 | grep -v "chromosome" | tr ',' "\t" | awk '{print $1 "\t" $2 "\t" $2+50}' | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 | awk '$3-$2 >=300{print $0}' >js.bed 
#785

# ri
cat <(cat ../*bed) | awk '$3-$2 >= 300 {print $0}' >ri.bed 
# 459

# 60% overlap
bedtools intersect -a ri.bed -b js.bed| wc
# 272


# 
python diffER.py \
    --genome_file ./data/CHROMSIZES_hg38_chr10.txt \
    --group_A_beds ./data/uCLL/*.bed \
    --group_B_beds ./data/mCLL/*.bed 

ri: chr10 68650 68700 2/5 0/13 0.1105263157894737
js: chr10,68650,      3/7 13/13 0.00722394220846233

js:
input_matrix
     [,1] [,2]
[1,]   13    3 # intersect; in my data i find 5 intersect and not significant 
[2,]    0    4 # non-intersect
    #mCLL  #uCLL

> fisher.test(input_matrix)$p.value
[1] 0.001354489

Performing Fisher's exact test...
[[5, 2], [13, 0]]
# 0.1105263157894737

[[5, 2], [0, 13]]
# 0.001354489

#
cd /Users/rashedulislam/Documents/diffER/temp/
bedtools intersect -a chr10_68650_68700 -b ../data/mCLL/*bed -wo | wc -l
bedtools intersect -a chr10_68650_68700 -b ../data/uCLL/*bed -wo | wc -l

# only 50% overlap at chr10
less H3K27me3_cll_mutated_cll_unmutated_chr10_0.03 | grep -v "chromosome" | tr ',' "\t" | awk '{print $1 "\t" $2 "\t" $2+50}' | sort -k1,1 -k2,2n | bedtools merge -i stdin -d 100 | bedtools intersect -b stdin -a <(cat ../*bed) -u | wc -l


# js dir
/projects/epigenomics_assembly/jsteif/IHEC/postgresql_tables
/projects/epigenomics3/temp/jsteif

less /projects/epigenomics3/epigenomics3_results/users/jsteif/scripts/TERM3/reference_epigenome/binned_reference_epigenome_w_bamCoverage_V3/IHEC/sample_database_prep/run_primary_table_maker.sh

# table maker
less /projects/epigenomics3/epigenomics3_results/users/jsteif/scripts/TERM3/reference_epigenome/binned_reference_epigenome_w_bamCoverage_V3/IHEC/sample_database_prep/opt_CEMT_primary_table_maker.R

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
x <- c(13, 13, 3, 3)

# Extract p-value
p_value <- p_value_extractor(x)
print(p_value)
