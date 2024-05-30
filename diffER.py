from util import make_windows_build, make_windows_file, intersect_bedfiles, perform_fisher_test 
import argparse
import pandas as pd

# Create windows from a genome and intersect BED files
def diffER(genome_build, genome_file, group_A_beds, group_B_beds, window_size, out_file):
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

    # Merge bed files to have the number of intersections (and non-intersection) of samples for group_A and group_B
    file_1 = pd.read_csv('./temp/group_A_intersect.bed', sep='\t', header=None)
    file_2 = pd.read_csv('./temp/group_B_intersect.bed', sep='\t', header=None)
    columns_to_add = file_2.iloc[:, [3, 4]]
    result = pd.concat([file_1, columns_to_add], axis=1)
    result.to_csv('./temp/merged_group_A_B.bed', sep='\t', index=False, header=False)

    # Fisher's exact test
    merge = './temp/merged_group_A_B.bed'
    perform_fisher_test(merge, out_file)

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
    parser.add_argument("--window_size", type=int, default=50, help="Size of the windows")
    parser.add_argument("--out_file", required=False, help="Output file name")
    args = parser.parse_args()

    # Call the diffER function with parsed arguments
    diffER(args.genome_build, args.genome_file, args.group_A_beds, args.group_B_beds, args.window_size, args.out_file)
