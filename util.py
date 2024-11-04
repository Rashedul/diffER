import pybedtools
import os
import sys
import pandas as pd
from scipy.stats import fisher_exact
import shutil
import glob

"""
Create windows of specific size from a genome build 
"""

def make_windows_build(genome_build, window_size):
    print("Generating windows from genome build...")
    # Create a BED object
    a = pybedtools.BedTool()

    # Create windows of specified size
    windows = a.window_maker(genome=genome_build, w=window_size)
    output_bed = f"windows.bed"
    
    # Ensure the temp directory exists
    temp_dir = "temp"
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)
        print(f"Directory created: {temp_dir}")
    else:
        print(f"Directory already exists: {temp_dir}")

    # Save the windows to a new BED file in temp directory
    filepath = os.path.join("temp", output_bed)
    windows.saveas(filepath)

"""
Create windows of specific size from a genome file
"""

def make_windows_file(genome_file, window_size):
    print("Generating windows from genome file...")
    # Create a BED object
    a = pybedtools.BedTool()

    # Create windows of specified size
    windows = a.window_maker(g=genome_file, w=window_size)
    output_bed = f"windows.bed"
    
    # Save the windows to a new BED file in temp directory
    os.makedirs("temp", exist_ok=True)
    filepath = os.path.join("temp", output_bed)
    windows.saveas(filepath)

"""
Count number of bed files intersected (not-intersected) to the windows 
"""

def intersect_bedfiles(primary_bed, bed_files_pattern, output_filename, output_directory):
    print("Intersecting bed files with genomic windows...")
    
    # Read primary bed file and sort
    window = pd.read_csv(primary_bed, sep='\t', header=None)
    primary = pybedtools.BedTool.from_dataframe(window).sort()

    # Initialize DataFrame for the primary bed file
    primary_df = primary.to_dataframe()
    primary_df['intersect_count'] = 0
    primary_df['non_intersect_count'] = 0

    # Process each bed file individually
    for bedfile in bed_files_pattern:
        bed = pybedtools.BedTool(bedfile).sort()
        intersected = primary.intersect(bed, c=True, sorted=True)
        intersected_df = intersected.to_dataframe()

        # Update the intersect count 
        primary_df['intersect_count'] += (intersected_df.iloc[:, -1] > 0).astype(int)

    # Calculate the non-intersect count
    primary_df['non_intersect_count'] = len(bed_files_pattern) - primary_df['intersect_count']

    # Convert DataFrame back to BedTool object
    result_bed = pybedtools.BedTool.from_dataframe(primary_df)

    # Determine the output file path
    output_filepath = os.path.join(output_directory, output_filename)

    # Save the result to the specified file
    result_bed.saveas(output_filepath)

    print(f"Intersection counts saved to {output_filepath}")

"""
Concate BED files of two groups  
"""

def merge_groups(file1, file2, merge_output, output_directory):
    print("Counting number of intersects...")
    # Read the BED files without considering the first row as column names
    file_1 = pd.read_csv(file1, sep='\t', header=None)
    file_2 = pd.read_csv(file2, sep='\t', header=None)
    
    # Select the 4th and 5th columns from file_2.tsv
    columns_to_add = file_2.iloc[:, [3, 4]]
    
    # Concatenate the columns to file_1
    result = pd.concat([file_1, columns_to_add], axis=1)

    # Remove rows where there are no intersections in both groups  
    result = result[(result.iloc[:, 3] != 0) & (result.iloc[:, 5] != 0)]
    output_filepath = os.path.join(output_directory, merge_output)
    result.to_csv(output_filepath, sep='\t', index=False, header=False)
    
    print(f"Merged output saved to {output_filepath}")

"""
Fisher's Exact Test 
"""

def perform_fisher_test(input_file, fisher_output):
    print("Performing Fisher's exact test...")
    # Read file without considering the first row as column names
    data = pd.read_csv(input_file, sep='\t', header=None)
    
    # Initialize lists to store p-values and odds ratios
    p_values = []
    odds_ratios = []
    
    # Iterate over each row in the DataFrame
    for index, row in data.iterrows():
        # Create a 2x2 contingency table from columns 4 to 7 (zero-indexed as 3 to 6)
        contingency_table = [[row[3], row[4]], [row[5], row[6]]]
        # print(contingency_table)
        
        # Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(contingency_table)
        
        # Append the results to the lists
        p_values.append(p_value)
        odds_ratios.append(odds_ratio)
    
    # Add the p-values and odds ratios as new columns to the DataFrame
    data[7] = p_values
    data[8] = odds_ratios
    
    # Save the resulting DataFrame to a new file
    data.to_csv(fisher_output, sep='\t', index=False, header=False)

"""
Merge neighboring enriched bins 
"""    

def enriched_regions(fisher_p_value, merge_intervals, differ_output_filename, output_directory):

    # Ensure the output directory exists
    os.makedirs(output_directory, exist_ok=True)

    # Read the file (assuming you are using pandas)
    df = pd.read_csv('./temp/merged_group_A_B_fisher.bed', sep='\t', header=None)

    # Filter rows based on the p-value threshold 
    df = df[df[7].astype(float) <= fisher_p_value]

    # Compute the new column using the formula
    df[9] = df[3] / (df[3] + df[4]) - df[5] / (df[5] + df[6])

    # Split the DataFrame into two groups
    df_pos = df[df[9] > 0]
    df_neg = df[df[9] < 0]

    # Merge neighboring bins
    df_pos_bed = pybedtools.BedTool.from_dataframe(df_pos)
    merged_intervals_pos = df_pos_bed.merge(d=merge_intervals)
    df_neg_bed = pybedtools.BedTool.from_dataframe(df_neg)
    merged_intervals_neg = df_neg_bed.merge(d=merge_intervals)

    # Construct output file paths
    output_file_pos = os.path.join(output_directory, f'{differ_output_filename}_group_A_enriched_regions.bed')
    output_file_neg = os.path.join(output_directory, f'{differ_output_filename}_group_B_enriched_regions.bed')

    # Write the merged intervals to BED files
    merged_intervals_pos.saveas(output_file_pos)
    merged_intervals_neg.saveas(output_file_neg)

    print(f"Output created and saved as: \n - {output_directory}/{differ_output_filename}_group_A_enriched_regions.bed \n - {output_directory}/{differ_output_filename}_group_B_enriched_regions.bed")

def remove_directory(dir_path):
    if os.path.exists(dir_path):
        try:
            shutil.rmtree(dir_path)
        except OSError as e:
            print(f"Error: {e.strerror}")
    else:
        print(f"Directory '{dir_path}' does not exist.")