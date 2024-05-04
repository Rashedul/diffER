import pybedtools
import genomepy

# def make_windows(genome_build, window_size):
#     # Load the input (temp) BED file
a = pybedtools.example_bedtool('a.bed')
g = 'genome.txt'
genome_build = BedTool(g)
window_size = 1000


    # Create windows of specified size
windows = a.window_maker(genome = genome_build, w = window_size)
for i in windows:
    print(i)

#     # Save the windows to a new BED file
#     output_bed = f"{genome_build}_windows.bed"
#     windows.saveas(output_bed)

#     print(f"Windows created and saved to {output_bed}")

# if __name__ == "__main__":
#     import argparse

#     # Parsing command line arguments
#     parser = argparse.ArgumentParser(description="Create windows from a BED file")
#     parser.add_argument("--genome_build", required=True, help="Input genome build name")
#     parser.add_argument("--window_size", type=int, default=1000, help="Size of the windows")
#     args = parser.parse_args()

#     make_windows(args.genome_build, args.window_size)
