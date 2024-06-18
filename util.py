import pybedtools
import os
import sys
import pandas as pd
from scipy.stats import fisher_exact
import shutil

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

def ensure_consistent_columns(filepath):
    with open(filepath, 'r') as file:
        lines = file.readlines()

    num_columns = len(lines[0].strip().split())
    
    parsed_lines = []
    for line in lines:
        columns = line.strip().split()[:3]  # Keep only the first three columns

        # Check for missing values
        if len(columns) < 3:
            raise ValueError(f"Missing values in line: {line.strip()}")
        
        # Ensure the second and third columns are integers (start and end positions)
        try:
            start = int(columns[1])
            end = int(columns[2])
        except ValueError:
            raise ValueError(f"Second or third column contains non-integer value in line: {line.strip()}")

        # Check for negative values
        if start < 0 or end < 0:
            raise ValueError(f"Negative values found in line: {line.strip()}")
        
        parsed_lines.append((columns[0], start, end))
    
    # Sort lines first by chromosome, then by start, then by end
    sorted_lines = sorted(parsed_lines, key=lambda x: (x[0], x[1], x[2]))

    with open(filepath, 'w') as file:
        for chromosome, start, end in sorted_lines:
            file.write(f"{chromosome}\t{start}\t{end}\n")


def intersect_bedfiles(primary_bed, multiple_beds, output_filename):
    print("Intersecting bed files with genomic windows...")

    # Ensure the temporary directory exists
    os.makedirs("temp", exist_ok=True)

    # Ensure primary_bed has consistent columns
    ensure_consistent_columns(primary_bed)
    
    # Read primary bed file
    window = pd.read_csv(primary_bed, sep='\t', header=None)
    primary = pybedtools.BedTool.from_dataframe(window)

    # Initialize DataFrame for the primary bed file
    primary_df = window.copy()
    primary_df['intersect_count'] = 0
    primary_df['non_intersect_count'] = 0

    # Process each multiple_bed file
    for bedfile in multiple_beds:
        # Ensure bedfile has consistent columns
        ensure_consistent_columns(bedfile)
        
        # Read each multiple bed file
        bed = pybedtools.BedTool(bedfile)

        # Intersect primary with the multiple bed file
        intersected = primary.intersect(bed, c=True, sorted=True)
        
        # Check if intersected has any lines
        # if len(intersected) == 0:
        #     continue
        
        intersected_df = intersected.to_dataframe()

        # Check if intersected_df has the expected column
        if intersected_df.shape[1] <= 0:
            raise ValueError(f"Intersected DataFrame from {bedfile} has no columns.")

        # Update the intersect count
        primary_df['intersect_count'] += (intersected_df.iloc[:, -1] > 0).astype(int)

    # Calculate the non-intersect count
    primary_df['non_intersect_count'] = len(multiple_beds) - primary_df['intersect_count']

    # Convert DataFrame back to BedTool object
    result_bed = pybedtools.BedTool.from_dataframe(primary_df)

    # Save the result to the specified file
    filepath = os.path.join("temp", output_filename)
    result_bed.saveas(filepath)

    print(f"Intersection results saved to {filepath}")

"""
Concate BED files of two groups  
"""

def merge_groups(file1, file2, merge_output):
    print("Counting number of intersects...")
    # Read the BED files without considering the first row as column names
    file_1 = pd.read_csv(file1, sep='\t', header=None)
    file_2 = pd.read_csv(file2, sep='\t', header=None)
    
    # Select the 4th and 5th columns from file_2.tsv
    columns_to_add = file_2.iloc[:, [3, 4]]
    
    # Concatenate the columns to file_1
    result = pd.concat([file_1, columns_to_add], axis=1)
    
    # Save the resulting DataFrame to a new BED file
    result.to_csv(merge_output, sep='\t', index=False, header=False)

"""
Fisher's Exact Test 
"""

def perform_fisher_test(input_file, fisher_output):
    print("Performing Fisher's exact test...")
    # Read the TSV file without considering the first row as column names
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
    
    # Save the resulting DataFrame to a new TSV file
    data.to_csv(fisher_output, sep='\t', index=False, header=False)

"""
Merge neighboring enriched bins 
"""    
def enriched_regions(fisher_p_value, merge_intervals, differ_output_filename):
    # Read the TSV file (assuming you are using pandas)
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

    # Write the merged intervals to BED files
    merged_intervals_pos.saveas(f'{differ_output_filename}_group_A_enriched_regions.bed')
    merged_intervals_neg.saveas(f'{differ_output_filename}_group_B_enriched_regions.bed')

    print(f"Output created and saved as: \n - {differ_output_filename}_group_A_enriched_regions.bed \n - {differ_output_filename}_group_B_enriched_regions.bed")

def remove_directory(dir_path):
    if os.path.exists(dir_path):
        try:
            shutil.rmtree(dir_path)
        except OSError as e:
            print(f"Error: {e.strerror}")
    else:
        print(f"Directory '{dir_path}' does not exist.")