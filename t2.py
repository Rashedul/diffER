import pybedtools

def intersect_multiple_bedfiles(primary_bed, *bed_files, set_b_bed_files=None, output_file="intersected_counts.bed"):
    """
    Intersect multiple sets of BED files with a primary BED file and save the count of intersecting BED files.

    Parameters:
        primary_bed (str): Path to the primary BED file.
        bed_files (str): Paths to the BED files in set A to intersect with the primary BED file.
        set_b_bed_files (str): Paths to the BED files in set B.
        output_file (str): Path to save the intersected counts BED file.

    Returns:
        None
    """
    # Load the primary BED file
    primary = pybedtools.BedTool(primary_bed)

    # Initialize dictionaries to store intersect counts for each genomic interval
    intersect_counts_set_a = {}
    intersect_counts_set_b = {}

    # Function to update intersect counts for a set of BED files
    def update_intersect_counts(intersect_counts, bed_files):
        for bed_file in bed_files:
            intersected = primary.intersect(bed_file, c=True)
            for feature in intersected:
                interval = tuple(feature.fields[:3])
                count = int(feature.fields[-1])
                intersect_counts[interval] = intersect_counts.get(interval, 0) + count

    # Update intersect counts for set A
    update_intersect_counts(intersect_counts_set_a, bed_files)

    # Update intersect counts for set B if provided
    if set_b_bed_files:
        update_intersect_counts(intersect_counts_set_b, set_b_bed_files)

    # Write the combined counts to the output file
    with open(output_file, "w") as outfile:
        for interval, count_a in intersect_counts_set_a.items():
            # Get the count for set B if available, otherwise set it to 0
            count_b = intersect_counts_set_b.get(interval, 0)
            # Append counts of intersecting BED files for set A and set B as the fifth and sixth columns
            outfile.write("\t".join(interval + (str(count_a), str(count_b))) + "\n")

    print(f"Intersected counts BED file saved as {output_file}")

if __name__ == "__main__":
    import argparse

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Intersect multiple sets of BED files with a primary BED file and save the count of intersecting BED files")
    parser.add_argument("primary_bed", help="Path to the primary BED file")
    parser.add_argument("bed_files", nargs="+", help="Paths to the BED files in set A")
    parser.add_argument("--set_b_bed_files", nargs="+", help="Paths to the BED files in set B")
    parser.add_argument("--output_file", default="intersected_counts.bed", help="Path to save the intersected counts BED file")
    args = parser.parse_args()

    intersect_multiple_bedfiles(args.primary_bed, *args.bed_files, set_b_bed_files=args.set_b_bed_files, output_file=args.output_file)

