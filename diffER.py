import argparse
from util import (
    make_windows_build, 
    make_windows_file, 
    intersect_bedfiles, 
    perform_fisher_test, 
    merge_groups, 
    enriched_regions,
    remove_directory
) 

def diffER(genome_build, genome_file, group_A_beds, group_B_beds, window_size, p_value, distance, outfile):
    # Create windows of specified size
    if genome_build: 
        make_windows_build(genome_build, window_size)
    elif genome_file:
        make_windows_file(genome_file, window_size)
    else:
        raise ValueError("Either genome_build or genome_file must be provided.")

    # Intersect primary BED with group A and B BED files
    intersect_bedfiles('./temp/windows.bed', group_A_beds, "group_A_intersect.bed")
    intersect_bedfiles('./temp/windows.bed', group_B_beds, "group_B_intersect.bed")

    # Merge two files containing the number of intersections (and non-intersection) of samples in group_A and group_B
    merge_groups('./temp/group_A_intersect.bed', './temp/group_B_intersect.bed', './temp/merged_group_A_B.bed')
    
    # Fisher's exact test
    perform_fisher_test('./temp/merged_group_A_B.bed', './temp/merged_group_A_B_fisher.bed')

    # Generate ERs per group
    enriched_regions(p_value, distance, outfile)

    # Remove intermediate files
    remove_directory('temp')

## test
# cd /home/rashedul/project/diffER
# less temp/merged_group_A_B_fisher.bed | awk '$8<.5{print $0, $4/($4+$5) - $6/($6+$7)}' | awk '$10 > 0{print}' | bedtools merge -i stdin -d 100 >group_A_enriched.bed
# less temp/merged_group_A_B_fisher.bed | awk '$8<.5{print $0, $4/($4+$5) - $6/($6+$7)}' | awk '$10 < 0{print}' | bedtools merge -i stdin -d 100 >group_B_enriched.bed
# ##

def main():
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Differentially Enriched Regions (diffER)")
    parser.add_argument("--genome_build", required=False, help="Input genome build name")
    parser.add_argument("--genome_file", required=False, help="Input genome file")
    parser.add_argument("--group_A_beds", required=True, nargs='+', help="Input group A BED files")
    parser.add_argument("--group_B_beds", required=True, nargs='+', help="Input group B BED files")
    parser.add_argument("--window_size", type=int, default=50, help="Size of the windows; default [50]")
    parser.add_argument("--p_value", type=float, default=0.05, help="p-value threshold for Fisher's exact test; default [0.05]")
    parser.add_argument("--distance", type=int, default=100, help=" Maximum distance between intervals allowed to be merged; default [100]")
    parser.add_argument("--outfile", required=False, default="diffER", help="Output file name")
    args = parser.parse_args()

    # Call the diffER function with parsed arguments
    diffER(args.genome_build, 
        args.genome_file, 
        args.group_A_beds, 
        args.group_B_beds, 
        args.window_size, 
        args.p_value, 
        args.distance, 
        args.outfile)

if __name__ == "__main__":
    main()