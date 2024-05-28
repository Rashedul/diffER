import pybedtools

def intersect_multiple_bedfiles(primary_bed, *bed_files, output_file="intersected_counts.bed"):
    """
    Intersect multiple BED files with a primary BED file and save the count of intersecting BED files.

    Parameters:
        primary_bed (str): Path to the primary BED file.
        bed_files (str): Paths to the additional BED files to intersect with the primary BED file.
        output_file (str): Path to save the intersected counts BED file.

    Returns:
        None
    """
    # Load the primary BED file
    primary = pybedtools.BedTool(primary_bed)

    # Initialize a dictionary to store intersect counts for each genomic interval
    intersect_counts = {}

    # Iterate through additional BED files
    for bed_file in bed_files:
        # Intersect the primary BED file with each additional BED file and count the number of overlaps
        intersected = primary.intersect(bed_file, c=True)
        # Store the counts for each genomic interval
        for feature in intersected:
            interval = tuple(feature.fields[:3])
            count = int(feature.fields[-1])
            intersect_counts[interval] = intersect_counts.get(interval, 0) + count

    # Write the combined counts to the output file
    with open(output_file, "w") as outfile:
        for interval, count in intersect_counts.items():
            # Append the count of intersecting BED files as the fourth column
            outfile.write("\t".join(interval + (str(count),)) + "\n")

    print(f"Intersected counts BED file saved as {output_file}")

if __name__ == "__main__":
    import argparse

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Intersect multiple BED files with a primary BED file and save the count of intersecting BED files")
    parser.add_argument("primary_bed", help="Path to the primary BED file")
    parser.add_argument("bed_files", nargs="+", help="Paths to additional BED files")
    parser.add_argument("--output_file", default="intersected_counts.bed", help="Path to save the intersected counts BED file")
    args = parser.parse_args()

    intersect_multiple_bedfiles(args.primary_bed, *args.bed_files, output_file=args.output_file)

