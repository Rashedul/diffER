import argparse
import multiprocessing
import tempfile
import os
from util import (
    make_windows_build,
    make_windows_file,
    intersect_bedfiles, 
    perform_fisher_test, 
    merge_groups, 
    enriched_regions,
    remove_directory
)


def intersect_bedfiles_worker(args):
    primary_bed, group_beds, output_file = args
    intersect_bedfiles(primary_bed, group_beds, output_file)

def perform_fisher_test_worker(args):
    input_file, output_file_chunk = args
    perform_fisher_test(input_file, output_file_chunk)

def split_file(input_file, num_chunks):
    with open(input_file, 'r') as f:
        lines = f.readlines()
    chunk_size = len(lines) // num_chunks
    chunks = [lines[i:i + chunk_size] for i in range(0, len(lines), chunk_size)]
    return chunks

def write_chunks(chunks, base_filename):
    chunk_files = []
    for i, chunk in enumerate(chunks):
        chunk_filename = f"{base_filename}_chunk_{i}.bed"
        with open(chunk_filename, 'w') as f:
            f.writelines(chunk)
        chunk_files.append(chunk_filename)
    return chunk_files

def combine_files(output_files, combined_filename):
    with open(combined_filename, 'w') as outfile:
        for fname in output_files:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
            os.remove(fname)  # Remove chunk file after combining

def diffER(genome_build, genome_file, group_A_beds, group_B_beds, window_size, p_value, distance, outfile, threads):
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create windows of specified size
        if genome_build: 
            make_windows_build(genome_build, window_size, temp_dir)
        elif genome_file:
            make_windows_file(genome_file, window_size, temp_dir)
        else:
            raise ValueError("Either genome_build or genome_file must be provided.")

        primary_bed = os.path.join(temp_dir, 'windows.bed')

        # Use multiprocessing to intersect primary BED with group A and B BED files
        pool = multiprocessing.Pool(processes=threads)

        intersect_tasks = [
            (primary_bed, group_A_beds, os.path.join(temp_dir, 'group_A_intersect.bed')),
            (primary_bed, group_B_beds, os.path.join(temp_dir, 'group_B_intersect.bed'))
        ]
        
        pool.map(intersect_bedfiles_worker, intersect_tasks)
        pool.close()
        pool.join()

        # Merge two files containing the number of intersections (and non-intersection) of samples in group_A and group_B
        merged_file = os.path.join(temp_dir, 'merged_group_A_B.bed')
        merge_groups(intersect_tasks[0][2], intersect_tasks[1][2], merged_file)

        # Split merged file into chunks for parallel processing
        chunks = split_file(merged_file, threads)
        chunk_files = write_chunks(chunks, os.path.join(temp_dir, 'merged_group_A_B'))

        # Perform Fisher's exact test in parallel
        pool = multiprocessing.Pool(processes=threads)
        fisher_tasks = [(chunk_file, f'{chunk_file}_fisher.bed') for chunk_file in chunk_files]
        pool.map(perform_fisher_test_worker, fisher_tasks)
        pool.close()
        pool.join()

        # Combine Fisher's test results
        fisher_output_files = [f'{chunk_file}_fisher.bed' for chunk_file in chunk_files]
        combined_fisher_output = os.path.join(temp_dir, 'merged_group_A_B_fisher.bed')
        combine_files(fisher_output_files, combined_fisher_output)

        # Generate ERs per group
        enriched_regions(p_value, distance, outfile)

def main():
    # Parsing command line arguments
    parser = argparse.ArgumentParser(description="Differentially Enriched Regions (diffER)")
    parser.add_argument("--genome_build", required=False, help="Input genome build name")
    parser.add_argument("--genome_file", required=False, help="Input genome file")
    parser.add_argument("--group_A_beds", required=True, nargs='+', help="Input group A BED files")
    parser.add_argument("--group_B_beds", required=True, nargs='+', help="Input group B BED files")
    parser.add_argument("--window_size", type=int, default=50, help="Size of the windows; default [50]")
    parser.add_argument("--p_value", type=float, default=0.05, help="p-value threshold for Fisher's exact test; default [0.05]")
    parser.add_argument("--distance", type=int, default=100, help=" Maximum distance between intervals allowed to be merged; default [100]")
    parser.add_argument("--outfile", required=False, default="diffER", help="Output file prefix")
    parser.add_argument("--threads", type=int, default=2, help="Number of threads (processors) to use for parallel processing; default [2]")
    args = parser.parse_args()

    # Call the diffER function with parsed arguments
    diffER(args.genome_build, 
           args.genome_file, 
           args.group_A_beds, 
           args.group_B_beds, 
           args.window_size, 
           args.p_value, 
           args.distance, 
           args.outfile, 
           args.threads)

if __name__ == "__main__":
    main()
