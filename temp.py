import pybedtools
import pandas as pd
import argparse

a = pybedtools.BedTool()
genome_file='./data/genome_chr1'
windows = a.window_maker(g=genome_file, w=1000)
output_bed = f"windows.bed"
windows.saveas(output_bed)

#######################
# test of intersection and count 

#!/usr/bin/env python3

import argparse
import pybedtools
import pandas as pd
import glob
import os

#python temp.py windows.bed './data/uCLL/*.bed' intersected_output.bed --output_directory mytemp2

def intersect_bedfiles(primary_bed, bed_files_pattern, output_filename, output_directory):
    print("Intersecting bed files with genomic windows...")
    
    # Read primary bed file and sort
    window = pd.read_csv(primary_bed, sep='\t', header=None)
    primary = pybedtools.BedTool.from_dataframe(window).sort()

    # Initialize DataFrame for the primary bed file
    primary_df = primary.to_dataframe()
    primary_df['intersect_count'] = 0
    primary_df['non_intersect_count'] = 0

    # Use glob to get all BED files matching the pattern
    multiple_beds = glob.glob(bed_files_pattern)

    # Process multiple bed files
    for bedfile in multiple_beds:
        bed = pybedtools.BedTool(bedfile).sort()
        intersected = primary.intersect(bed, c=True, sorted=True)
        intersected_df = intersected.to_dataframe()

        # Update the intersect count 
        primary_df['intersect_count'] += (intersected_df.iloc[:, -1] > 0).astype(int)

    # Calculate the non-intersect count
    primary_df['non_intersect_count'] = len(multiple_beds) - primary_df['intersect_count']

    # Convert DataFrame back to BedTool object
    result_bed = pybedtools.BedTool.from_dataframe(primary_df)

    # Determine the output file path
    output_filepath = os.path.join(output_directory, output_filename)

    # Save the result to the specified file
    result_bed.saveas(output_filepath)

    print(f"Intersection results saved to {output_filepath}")

def main():
    parser = argparse.ArgumentParser(description='Intersect BED files with genomic windows.')
    parser.add_argument('primary_bed', help='Path to the primary BED file.')
    parser.add_argument('bed_files_pattern', help='Pattern to match multiple BED files.')
    parser.add_argument('output_filename', help='Filename for the intersected output.')
    parser.add_argument('--output_directory', '-o', default='.', help='Directory to save the output file (default: current directory).')

    args = parser.parse_args()

    intersect_bedfiles(args.primary_bed, args.bed_files_pattern, args.output_filename, args.output_directory)

if __name__ == '__main__':
    main()

# mkdir mytemp
# python temp.py windows.bed './data/uCLL/*.bed' intersected_output.bed --output_directory mytemp
# tested and worked

#######################
# test of sort function 
# # Load BED file into a BedTool object
# bed = pybedtools.BedTool('input.bed')

# # Sort the BedTool object
# sorted_bed = bed.sort()

# # Save the sorted intervals to a new file
# sorted_bed.saveas('sorted_output.bed')
