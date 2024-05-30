#!/usr/bin/env python3

import pandas as pd
from scipy.stats import fisher_exact
import sys

def perform_fisher_test(input_file, output_file):
    # Read the TSV file without considering the first row as column names
    data = pd.read_csv(input_file, sep='\t', header=None)
    
    # Initialize lists to store p-values and odds ratios
    p_values = []
    odds_ratios = []
    
    # Iterate over each row in the DataFrame
    for index, row in data.iterrows():
        # Create a 2x2 contingency table from columns 4 to 7 (zero-indexed as 3 to 6)
        contingency_table = [[row[3], row[4]], [row[5], row[6]]]
        print(contingency_table)
        
        # Perform Fisher's exact test
        odds_ratio, p_value = fisher_exact(contingency_table)
        
        # Append the results to the lists
        p_values.append(p_value)
        odds_ratios.append(odds_ratio)
    
    # Add the p-values and odds ratios as new columns to the DataFrame
    data[7] = p_values
    data[8] = odds_ratios
    
    # Save the resulting DataFrame to a new TSV file
    data.to_csv(output_file, sep='\t', index=False, header=False)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input_file.tsv> <output_file.tsv>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    perform_fisher_test(input_file, output_file)

