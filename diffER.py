from util import make_windows_build, make_windows_file, intersect_bedfiles, perform_fisher_test, merge_groups, enriched_regions
import argparse

# Create windows from a genome and intersect BED files
def diffER(genome_build, genome_file, group_A_beds, group_B_beds, window_size, p_value, distance):
    # Create windows of specified size
    if genome_build: 
        make_windows_build(genome_build, window_size)
    elif genome_file:
        make_windows_file(genome_file, window_size)
    else:
        raise ValueError("Either genome_build or genome_file must be provided.")

    # Intersect primary BED with group A and B BED files
    primary_bed = './temp/windows.bed'
    intersect_bedfiles(primary_bed, group_A_beds, "group_A_intersect.bed")
    intersect_bedfiles(primary_bed, group_B_beds, "group_B_intersect.bed")

    # Number of intersections (and non-intersection) of samples in both group_A and group_B
    group_A = './temp/group_A_intersect.bed'
    group_B = './temp/group_B_intersect.bed'
    output_file_groups = './temp/merged_group_A_B.bed'
    merge_groups(group_A, group_B, output_file_groups)
    
    # Fisher's exact test
    merge = './temp/merged_group_A_B.bed'
    output_file_fisher = './temp/merged_group_A_B_fisher.bed'
    perform_fisher_test(merge, output_file_fisher)
    # 4/(4+5) - 6/(6+7); if +ve group_A, if -ve group_B
    # split the file for up and dn regions
    # merge neighboring bins by 100bp

    # Generates ERs per group
    enriched_regions(p_value, distance)

## test
#cd /home/rashedul/project/diffER
#less temp/merged_group_A_B_fisher.bed | awk '$8<.5{print $0, $4/($4+$5) - $6/($6+$7)}' | awk '$10 > 0{print}' | bedtools merge -i stdin -d 100 >group_A_enriched.bed
#less temp/merged_group_A_B_fisher.bed | awk '$8<.5{print $0, $4/($4+$5) - $6/($6+$7)}' | awk '$10 < 0{print}' | bedtools merge -i stdin -d 100 >group_B_enriched.bed
##

# def main():
#     # Check for correct number of arguments
#     if len(sys.argv) < 4:
#         print(f"Usage: {sys.argv[0]} <primary_bed_file> <multiple_bed_file1> [<multiple_bed_file2> ...] <output_filename>")
#         sys.exit(1)

if __name__ == "__main__":
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Differentially Enriched Regions (diffER)")
    parser.add_argument("--genome_build", required=False, help="Input genome build name")
    parser.add_argument("--genome_file", required=False, help="Input genome file")
    parser.add_argument("--group_A_beds", required=True, nargs='+', help="Input group A BED files")
    parser.add_argument("--group_B_beds", required=True, nargs='+', help="Input group B BED files")
    parser.add_argument("--window_size", type=int, default=50, help="Size of the windows; default [50]")
    parser.add_argument("--p_value", type=float, default=0.05, help="p-value threshold for Fisher's exact test; default [0.05]")
    parser.add_argument("--distance", type=int, default=200, help=" Maximum distance between intervals allowed to be merged; default [200]")
    args = parser.parse_args()

    # Call the diffER function with parsed arguments
    diffER(args.genome_build, args.genome_file, args.group_A_beds, args.group_B_beds, args.window_size, args.p_value, args.distance)
