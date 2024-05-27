import pybedtools

def intersect_multiple_bedfiles(primary_bed, group_A, group_B=None, output_file="intersected_counts.bed"):
    """
    Intersect multiple sets of BED files with a primary BED file and save the count of intersecting BED files.

    Parameters:
        primary_bed (str): Path to the primary BED file.
        group_A (str): Paths to the BED files in group A to intersect with the primary BED file.
        group_B (str): Paths to the BED files in group B.
        output_file (str): Path to save the intersected counts BED file.

    Returns:
        None
    """
    # Load the primary BED file
    primary = pybedtools.BedTool(primary_bed)

    # Initialize dictionaries to store intersect counts for each genomic interval
    intersect_counts = {}

    # Function to initialize intersect counts for a set of BED files
    def initialize_intersect_counts(bed_files, column):
        for bed_file in bed_files:
            for feature in pybedtools.BedTool(bed_file):
                interval = tuple(feature.fields[:3])
                if interval not in intersect_counts:
                    intersect_counts[interval] = [0, 0]
                intersect_counts[interval][column] = 0

    # Initialize intersect counts for group_A
    initialize_intersect_counts(group_A, 0)

    # Initialize intersect counts for group_B if provided
    if group_B:
        initialize_intersect_counts(group_B, 1)

    # Function to update intersect counts for a set of BED files
    def update_intersect_counts(bed_files, column):
        for bed_file in bed_files:
            intersected = primary.intersect(bed_file, c=True)
            for feature in intersected:
                interval = tuple(feature.fields[:3])
                count = int(feature.fields[-1])
                intersect_counts[interval][column] += count

    # Update intersect counts for group_A
    update_intersect_counts(group_A, 0)

    # Update intersect counts for group_B if provided
    if group_B:
        update_intersect_counts(group_B, 1)

    # Write the combined counts to the output file
    with open(output_file, "w") as outfile:
        for interval, counts in intersect_counts.items():
            # Append count of intersecting BED files for group_A and group_B as the fourth and fifth columns
            outfile.write("\t".join(interval + (str(counts[0]), str(counts[1]))) + "\n")

    print(f"Intersected counts BED file saved as {output_file}")

if __name__ == "__main__":
    import argparse

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Intersect multiple sets of BED files with a primary BED file and save the count of intersecting BED files")
    parser.add_argument("primary_bed", help="Path to the primary BED file")
    parser.add_argument("--group_A", nargs="+", required=True, help="Paths to the BED files in group A to intersect with the primary BED file")
    parser.add_argument("--group_B", nargs="+", required=True, help="Paths to the BED files in group B")
    parser.add_argument("--output_file", default="intersected_counts.bed", help="Path to save the intersected counts BED file")
    args = parser.parse_args()

    intersect_multiple_bedfiles(args.primary_bed, args.group_A, group_B=args.group_B, output_file=args.output_file)

