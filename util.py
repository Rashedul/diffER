import pybedtools

def make_windows_build(genome_build, window_size):

# Load the input (temp) BED file
    a = pybedtools.example_bedtool('a.bed')

    # Create windows of specified size
    windows = a.window_maker(genome=genome_build, w=window_size)
    output_bed = f"windows.bed"
    
    # Save the windows to a new BED file
    windows.saveas(output_bed)
    print(f"Windows created and saved to {output_bed}")

def make_windows_file(genome_file, window_size):

# Load the input (temp) BED file
    a = pybedtools.example_bedtool('a.bed')

    # Create windows of specified size
    windows = a.window_maker(g=genome_file, w=window_size)
    output_bed = f"windows.bed"
    
    # Save the windows to a new BED file
    windows.saveas(output_bed)
    print(f"Windows created and saved to {output_bed}")
