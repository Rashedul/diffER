import pandas as pd

# Read the TSV files
file_1 = pd.read_csv('./temp/group_A_intersect.bed', sep='\t')
file_2 = pd.read_csv('./temp/group_A_intersect.bed', sep='\t')

# Print the first few rows and data types for debugging
print("File 1 - Initial:")
print(file_1.head())
print(file_1.dtypes)

print("File 2 - Initial:")
print(file_2.head())
print(file_2.dtypes)

# Select the 4th and 5th columns from file_2.tsv
columns_to_add = file_2.iloc[:, [3, 4]]

# Print the selected columns for debugging
print("Columns to Add:")
print(columns_to_add.head())
print(columns_to_add.dtypes)

# Concatenate the columns to file_1
result = pd.concat([file_1, columns_to_add], axis=1)

# Print the resulting DataFrame for debugging
print("Resulting DataFrame:")
print(result.head())
print(result.dtypes)

# Save the resulting DataFrame to a new TSV file
result.to_csv('./temp/merged_group_A_B.bed', sep='\t', index=False)

