import pybedtools
import pandas as pd
import argparse
import glob
import os
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches  # For creating legend handles

def main(intervals_path, group_A_pattern, group_B_pattern, output_file, image_file):
    # Define the main intervals file and load it as a BedTool object
    main_intervals = pybedtools.BedTool(intervals_path)

    # Expand the sample file patterns for group A and group B
    group_A_files = glob.glob(group_A_pattern)
    group_B_files = glob.glob(group_B_pattern)

    # Check if files are found for each group
    if not group_A_files:
        print(f"No files found for Group A pattern: {group_A_pattern}")
        return
    if not group_B_files:
        print(f"No files found for Group B pattern: {group_B_pattern}")
        return

    # Combine the files and label them for the heatmap annotation
    sample_files = group_A_files + group_B_files
    sample_labels = ['Group A'] * len(group_A_files) + ['Group B'] * len(group_B_files)

    # Initialize a dictionary to store results (fraction of interval overlapped for each sample)
    results = {interval: [0] * len(sample_files) for interval in main_intervals}

    # Loop through each sample to calculate the fraction of interval overlapped
    for idx, sample_file in enumerate(sample_files):
        sample_bed = pybedtools.BedTool(sample_file)
        
        # Intersect sample with main intervals using the -wo option to get overlap length
        intersected = main_intervals.intersect(sample_bed, wo=True)
        
        # Calculate the fraction of enriched regions occupied by peaks
        for feature in intersected:
            # Identify the main interval from the intersection output
            interval = pybedtools.create_interval_from_list([feature[0], feature[1], feature[2]])
            
            # Get the overlap length (last field in the feature with -wo option)
            overlap_length = int(feature[-1])
            
            # Calculate the length of the interval
            interval_length = int(interval.end) - int(interval.start)
            
            # Add the fraction of overlap to the appropriate interval and sample
            results[interval][idx] += overlap_length / interval_length

    # Convert the results to a DataFrame
    interval_names = [f"{iv.chrom}:{iv.start}-{iv.end}" for iv in main_intervals]
    df = pd.DataFrame(results.values(), index=interval_names, columns=[os.path.basename(f) for f in sample_files])

    # Save the DataFrame to a CSV file
    df.to_csv(output_file, index=True)
    print(f"Data has been written to '{output_file}'")

    # Create an annotation bar for the heatmap
    annotation = pd.DataFrame(sample_labels, index=df.columns, columns=["Group"])
    lut = {"Group A": "skyblue", "Group B": "salmon"}  # Define colors for each group
    col_colors = annotation["Group"].map(lut)

    # Plot heatmap with annotation
    g = sns.clustermap(df, cmap="viridis", col_cluster=False, row_cluster=False, figsize=(10, 8), 
                       xticklabels=False, yticklabels=False, col_colors=col_colors)

    # Add y-axis label
    g.ax_heatmap.set_ylabel("Differentially Enriched Regions", labelpad=15, fontsize=20)
    g.ax_heatmap.yaxis.set_label_position("left")

    # Create custom legend
    legend_labels = [mpatches.Patch(color='skyblue', label='Group A'),
                     mpatches.Patch(color='salmon', label='Group B')]
    # plt.legend(handles=legend_labels, loc='lower left', fontsize=12, title="Groups", bbox_to_anchor=(0, -0.8))
    plt.legend(handles=legend_labels, loc='lower left', fontsize=12, title="Groups", bbox_to_anchor=(2, 0))

    # Save the heatmap as an image file
    plt.savefig(image_file, dpi=300, bbox_inches='tight')  # Save as specified image file with 300 DPI and tight bounding box
    print(f"Heatmap has been written to '{image_file}'")

    # Display the plot
    # plt.show()

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Calculate fraction of enriched regions covered by peaks using multiple BED files.")
    parser.add_argument("-i", "--intervals", required=True, help="Path to the main intervals BED file (e.g., diffER_group_A_enriched_regions.bed)")
    parser.add_argument("--group_A_beds", required=True, help="Glob pattern for Group A sample BED files (e.g., 'group_A/*.bed')")
    parser.add_argument("--group_B_beds", required=True, help="Glob pattern for Group B sample BED files (e.g., 'group_B/*.bed')")
    parser.add_argument("-o", "--output", default="interval_fraction_coverage.csv", help="Output CSV file name")
    parser.add_argument("-img", "--image", default="heatmap.png", help="Output image file name for the heatmap")

    # Parse arguments
    args = parser.parse_args()

    main(args.intervals, args.group_A_beds, args.group_B_beds, args.output, args.image)
