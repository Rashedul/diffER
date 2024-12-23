## diffER: Analysis of differentially enriched regions from ChIP-seq peaks

<table>
  <tr>
    <td>
      <h3> Differentially Enriched Regions (diffER)</h3>
      <p>Workflow:</p>
      <ol>
        <li> Split the genome into non-overlapping bins.</li>
        <li> Calculate the number of samples per group that have (or have not) peaks in each bin.</li>
        <li> Fisher’s exact test to identify the differentially enriched bins between two groups.</li>
		<li> Merge the neighboring bins that are differentially enriched.</li>
      </ol>
      <br>
    </td>
    <td>
      <img src="./data/pipeline.gif" alt="Project Image" width="400"/>
    </td>
  </tr>
</table>

### 1. Requirements 
- Python >= 3.10
- Linux
- macOS >= 13.5

### 2. Installation

 - Clone the repository

```
git clone https://github.com/Rashedul/diffER.git
```

 - Create a virtual environment

```
cd diffER
python -m venv environment
source ./environment/bin/activate
```

 - Install python packages

```
pip install -r requirements.txt
```

### 3. Run diffER

#### Note:
- Provide a `--genome_file` of your interest. Example of a genome file is provided [here](./data/CHROMSIZES_hg38.txt). 
- Alternatively, `--genome_build` such as `hg38`, `hg19`, `mm9`, `mm10` etc. can be used. Find the available genomes in `pybedtools` by searching the genome build such as `hg38`.
```
import pybedtools
pybedtools.genome_registry.hg38
or
pybedtools.chromosomesizes("hg38")
```
- `--group_A_beds` and `--group_B_beds` bed files can be provided as list and/or wildcard (`*`) character. 
- The required number of samples per group is at least 4.

#### Usage
```
python diffER.py -h

Usage:
    python diffER.py [--genome_build GENOME_BUILD | --genome_file GENOME_FILE] --group_A_beds GROUP_A_BEDS [--group_B_beds GROUP_B_BEDS] [--window_size WINDOW_SIZE] [--p_value P_VALUE] [--distance DISTANCE] [--outfile OUTFILE] [--outdir OUTDIR]  

Arguments:
  -h, --help            Show this help message and exit
  --genome_build GENOME_BUILD
                        Input genome build name such as hg38 [required, if genome_file is not provided]
  --genome_file GENOME_FILE
                        Input genome file [required, if genome_build is not provided]
  --group_A_beds GROUP_A_BEDS [GROUP_A_BEDS ...]
                        Input group A BED files (e.g., 'group_A/*.bed') [required]
  --group_B_beds GROUP_B_BEDS [GROUP_B_BEDS ...]
                        Input group B BED files (e.g., 'group_B/*.bed') [required]
  --window_size WINDOW_SIZE
                        Size of the windows to bin the genome; default [50]
  --p_value P_VALUE     P-value threshold for Fisher's exact test; default [0.03]
  --distance DISTANCE   Maximum distance (in base pairs) between differentially enriched regions to be merged; default [100]
  --outfile OUTFILE     Output file prefix; default [diffER]
  --outdir OUTDIR       Output directory name; default [current directory]
```

#### Example command with genome file (hg38:chr10)
```
python diffER.py \
    --genome_file ./data/CHROMSIZES_hg38_chr10.txt \
    --group_A_beds ./data/uCLL/*.bed \
    --group_B_beds ./data/mCLL/*.bed 
```

#### Example command with genome build 
```
python diffER.py \
    --genome_build hg38 \
    --group_A_beds ./data/uCLL/*.bed \
    --group_B_beds ./data/mCLL/*.bed
```

#### Output files
- For each of the two groups, there are two output bed files containing enriched regions. 
- The output files contain 3 columns (`chr` `start` `end`) 

```
  - diffER_group_A_enriched_regions.bed 
  - diffER_group_B_enriched_regions.bed
```

#### Note:

- Depending on the assay, you may want to exclude regions shorter than a specified length. For instance, we recommend filtering out regions under 300bp for broad marks.

### 4. Plot heatmap

- You may experiment with different p-values and check quality by visualizing the results in a heatmap or profile plot. You can use `plot_heatmap.py` to generate a heatmap. 

```
# Show help message
python plot_heatmap.py -h

# Example heatmaps
python plot_heatmap.py \
  -i diffER_group_A_enriched_regions.bed \
  --group_A_beds "./data/uCLL/*.bed" \
  --group_B_beds "./data/mCLL/*.bed" \
  -o output1.csv \
  -img heatmap1.png

python plot_heatmap.py \
  -i diffER_group_B_enriched_regions.bed \
  --group_A_beds "./data/uCLL/*.bed" \
  --group_B_beds "./data/mCLL/*.bed" \
  -o output2.csv \
  -img heatmap2.png
```

<div style="display: flex; justify-content: space-between;">
  <img src="./data/heatmap1.png" alt="Figure 1" style="width: 45%;"/>
  <img src="./data/heatmap2.png" alt="Figure 2" style="width: 45%;"/>
</div>

- Heatmaps for group_A (n = 7 samples) and group_B (n = 13 samples) enriched regions at chr10.
- The fraction of enriched regions occupied by peaks is used to generate heatmaps.
- `output1.csv` and `output2.csv` contain the data used in heatmap.

### 5. Contact  
Rashedul Islam, PhD (rashedul.gen@gmail.com)

### 6. Citation  
Islam R. et al., Integrative analysis of aberrant epigenomic landscape in chronic lymphocytic leukemia, (2024).

### 7. License
This project is licensed under the [MIT license](LICENSE). 