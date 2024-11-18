import argparse
import tempfile
from util import (
    make_windows_build, 
    make_windows_file, 
    intersect_bedfiles, 
    perform_fisher_test, 
    merge_groups, 
    enriched_regions,
)

def diffER(genome_build, genome_file, group_A_beds, group_B_beds, window_size, p_value, distance, outfile, outdir):
    # Use a temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        print(f"To store temporary files, using temporary directory: {temp_dir}")
        
        # Create windows of specified size
        if genome_build: 
            make_windows_build(genome_build, window_size, temp_dir)
        elif genome_file:
            make_windows_file(genome_file, window_size, temp_dir)
        else:
            raise ValueError("Either genome_build or genome_file must be provided.")

        # Intersect primary BED with group A and B BED files
        intersect_bedfiles(f"{temp_dir}/windows.bed", group_A_beds, "group_A_intersect.bed", temp_dir)
        intersect_bedfiles(f"{temp_dir}/windows.bed", group_B_beds, "group_B_intersect.bed", temp_dir)

        # Merge two files containing the number of intersections (and non-intersection) of samples in group_A and group_B
        merge_groups(f"{temp_dir}/group_A_intersect.bed", f"{temp_dir}/group_B_intersect.bed", 'merged_group_A_B.bed', temp_dir)
        
        # Perform Fisher's exact test
        perform_fisher_test(f"{temp_dir}/merged_group_A_B.bed", f"{temp_dir}/merged_group_A_B_fisher.bed")

        # Generate enriched regions per group
        enriched_regions(p_value, distance, outfile, outdir, temp_dir)

        # Temporary directory and its contents are automatically cleaned up here
        print("Processing complete. Temporary files removed.")

def main():
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Differentially Enriched Regions (diffER)")
    parser.add_argument("--genome_build", required=False, help="Input genome build name such as hg38 [required, if genome_file is not provided]")
    parser.add_argument("--genome_file", required=False, help="Input genome file [required, if genome_build is not provided]")
    parser.add_argument("--group_A_beds", required=True, nargs='+', help="Input group A BED files (e.g., 'group_A/*.bed') [required]")
    parser.add_argument("--group_B_beds", required=True, nargs='+', help="Input group B BED files (e.g., 'group_B/*.bed') [required]")
    parser.add_argument("--window_size", type=int, default=50, help="Size of the windows to bin the genome; default [50]")
    parser.add_argument("--p_value", type=float, default=0.03, help="P-value threshold for Fisher's exact test; default [0.03]")
    parser.add_argument("--distance", type=int, default=100, help="Maximum distance (in base pairs) between differentially enriched regions to be merged; default [100]")
    parser.add_argument("--outfile", required=False, default="diffER", help="Output file prefix; default [diffER]")
    parser.add_argument("--outdir", required=False, default=".", help="Output direcory name; default  [current direcory]")
    args = parser.parse_args()

    # Call the diffER function with parsed arguments
    diffER(args.genome_build, 
        args.genome_file, 
        args.group_A_beds, 
        args.group_B_beds, 
        args.window_size, 
        args.p_value, 
        args.distance, 
        args.outfile,
        args.outdir)

if __name__ == "__main__":
    main()