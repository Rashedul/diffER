import pybedtools

def make_windows(genome_build, genome_file, window_size):
    """
    Create windows from a BED file.

    Parameters:
        genome_build (str): Input genome build name.
        genome_file (str): Input genome file.
        window_size (int): Size of the windows.

    Raises:
        ValueError: If neither genome_build nor genome_file is provided.

    Returns:
        None
    """
    # Load the input (temp) BED file
    a = pybedtools.example_bedtool('a.bed')

    # Create windows of specified size
    if genome_build: 
        windows = a.window_maker(genome=genome_build, w=window_size)
        output_bed = f"{genome_build}_windows.bed"
    elif genome_file:
        windows = a.window_maker(g=genome_file, w=window_size)
        output_bed = f"{genome_file}_windows.bed"
    else:
        raise ValueError("Either genome_build or genome_file must be provided.")

    # Save the windows to a new BED file
    windows.saveas(output_bed)
    print(f"Windows created and saved to {output_bed}")

if __name__ == "__main__":
    import argparse

    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Create windows from a BED file")
    parser.add_argument("--genome_build", required=False, help="Input genome build name")
    parser.add_argument("--genome_file", required=False, help="Input genome file")
    parser.add_argument("--window_size", type=int, default=1000, help="Size of the windows")
    args = parser.parse_args()

    make_windows(args.genome_build, args.genome_file, args.window_size)

