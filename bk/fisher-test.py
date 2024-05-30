import pandas as pd

# Read the TSV files
file_1 = pd.read_csv('./temp/group_A_intersect.bed', sep='\t', header=None)
file_2 = pd.read_csv('./temp/group_B_intersect.bed', sep='\t', header=None)

# Select the 4th and 5th columns from file_2.tsv
# Note: Columns in pandas are zero-indexed, so 4th column is index 3 and 5th column is index 4
columns_to_add = file_2.iloc[:, [3, 4]]

# Concatenate the columns to file_1
result = pd.concat([file_1, columns_to_add], axis=1)

# Save the resulting DataFrame to a new TSV file
result.to_csv('./temp/merged_group_A_B.bed', sep='\t', index=False, header=False)
