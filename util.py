import pybedtools
import os
import sys
import pandas as pd

"""
Create windows of specific size from a genome build 
"""

def make_windows_build(genome_build, window_size):

# Load the input (temp) BED file
    a = pybedtools.BedTool()

    # Create windows of specified size
    windows = a.window_maker(genome=genome_build, w=window_size)
    output_bed = f"windows.bed"
    
    # Save the windows to a new BED file in temp directory
    os.makedirs("temp", exist_ok=True)
    filepath = os.path.join("temp", output_bed)
    windows.saveas(filepath)
    print(f"Windows created and saved to {output_bed}")

"""
Create windows of specific size from a genome file
"""

def make_windows_file(genome_file, window_size):

# Load the input (temp) BED file
    a = pybedtools.BedTool()

    # Create windows of specified size
    windows = a.window_maker(g=genome_file, w=window_size)
    output_bed = f"windows.bed"
    
    # Save the windows to a new BED file in temp directory
    os.makedirs("temp", exist_ok=True)
    filepath = os.path.join("temp", output_bed)
    windows.saveas(filepath)
    print(f"Windows created and saved to {output_bed}")

"""
Count number of bed files intersected (not-intersected) to the windows 
"""

def intersect_bedfiles(primary_bed, multiple_beds, output_filename):
    # Read primary bed file
    primary = pybedtools.BedTool(primary_bed)

    # Initialize DataFrame for the primary bed file
    primary_df = primary.to_dataframe()

    # Add columns for intersected and non-intersected counts
    primary_df['intersect_count'] = 0
    primary_df['non_intersect_count'] = 0

    # Process each multiple_bed file
    for bedfile in multiple_beds:
        # Read each multiple bed file
        bed = pybedtools.BedTool(bedfile)

        # Intersect primary with the multiple bed file
        intersected = primary.intersect(bed, c=True)
        intersected_df = intersected.to_dataframe()

        # Update the intersect count
        primary_df['intersect_count'] += intersected_df.iloc[:, -1]

    # Calculate the non-intersect count
    primary_df['non_intersect_count'] = len(multiple_beds) - primary_df['intersect_count']

    # Convert DataFrame back to BedTool object
    result_bed = pybedtools.BedTool.from_dataframe(primary_df)

    # Save the result to the specified file
    filepath = os.path.join("temp", output_filename)
    result_bed.saveas(filepath)

    return result_bed

    # Extract arguments
    primary_bed_file = sys.argv[1]
    output_filename = sys.argv[-1]
    multiple_bed_files = sys.argv[2:-1]

    # Perform intersection
    result_bed = intersect_bedfiles(primary_bed_file, multiple_bed_files, output_filename)

    print(f"Intersection result saved to {output_filename}")

"""
Fisher's Exact Test 
"""

