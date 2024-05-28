import argparse
import pandas as pd
import pybedtools
from scipy.stats import fisher_exact

def read_bed(file_path):
    """Reads a BED file and returns a pybedtools BedTool object."""
    return pybedtools.BedTool(file_path)

def count_overlaps(primary_bed, group_bed):
    """Counts overlaps between primary_bed and a group_bed."""
    return primary_bed.intersect(group_bed, c=True)

def fishers_exact_test(primary_bed, group_A_beds, group_B_beds):
    """Perform Fisher's exact test for each region in primary_bed based on overlap counts."""
    primary = read_bed(primary_bed)
    group_A = [read_bed(bed) for bed in group_A_beds]
    group_B = [read_bed(bed) for bed in group_B_beds]

    results = []

    for interval in primary:
        row = {
            'chrom': interval.chrom,
            'start': interval.start,
            'end': interval.end
        }

        # Count overlaps with group A and B
        overlap_counts_A = [count_overlaps(interval, group_bed).to_dataframe().shape[0] for group_bed in group_A]
        overlap_counts_B = [count_overlaps(interval, group_bed).to_dataframe().shape[0] for group_bed in group_B]

        # Calculate Fisher's exact test
        count_A = sum(overlap_counts_A)
        count_B = sum(overlap_counts_B)

        contingency_table = [
            [count_A, sum(overlap_counts_A) - count_A],
            [count_B, sum(overlap_counts_B) - count_B]
        ]

        oddsratio, pvalue = fisher_exact(contingency_table)

        # Add results to the row
        row.update({
            'overlap_A': sum(overlap_counts_A),
            'overlap_B': sum(overlap_counts_B),
            'count_A': count_A,
            'count_B': count_B,
            'odds_ratio': oddsratio,
            'p_value': pvalue
        })

        results.append(row)

    return pd.DataFrame(results)

def save_results(results, output_file):
    """Saves results DataFrame to a TSV file."""
    results.to_csv(output_file, sep='\t', index=False)

def main(primary_bed, group_A, group_B, output_file=None):
    results = fishers_exact_test(primary_bed, group_A, group_B)
    print("Fisher's exact test results:")
    print(results)
    
    if output_file:
        save_results(results, output_file)
        print(f"\nResults saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform Fisher\'s exact test for primary_bed based on overlap counts with group_A and group_B.')
    parser.add_argument('primary_bed', type=str, help='Primary BED file')
    parser.add_argument('--group_A', nargs='+', type=str, required=True, help='List of group A BED files')
    parser.add_argument('--group_B', nargs='+', type=str, required=True, help='List of group B BED files')
    parser.add_argument('--output_file', type=str, help='Output TSV file path')
    args = parser.parse_args()

    main(args.primary_bed, args.group_A, args.group_B, args.output_file)
