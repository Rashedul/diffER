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
		<li> Merged the neighboring bins that are differentially enriched.</li>
      </ol>
      <br>
    </td>
    <td>
      <img src="./data/pipeline.gif" alt="Project Image" width="400"/>
    </td>
  </tr>
</table>

### 1. Requirements 
- python >= 3.10
- Linux

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

#### Notes
- Examples of `genome_build` are `hg38`, `hg19`, `mm9`, `mm10` etc. Check available genomes in `genomepy`. Alternatively provide a `genome file` of your interest. Example of `genome file` is provided [here](./data/genome_file). 
- `group_A_beds` and `group_B_beds` bed files can be provided as list and/or wildcard (`*`) character. 
- Required number of samples per group is at least 4.

#### Usage
```
python diffER.py -h

Usage:
    python diffER.py [--genome_build GENOME_BUILD | --genome_file GENOME_FILE] --group_A_beds GROUP_A_BEDS [--group_B_beds GROUP_B_BEDS] [--window_size WINDOW_SIZE] [--p_value P_VALUE] [--distance DISTANCE]

Arguments:
    --genome_build GENOME_BUILD   Input genome build name
    --genome_file GENOME_FILE     Input genome file
    --group_A_beds GROUP_A_BEDS   Input group A BED files (multiple allowed)
    --group_B_beds GROUP_B_BEDS   Input group B BED files (multiple allowed)
    --window_size WINDOW_SIZE     Size of the windows; default [50]
    --p_value P_VALUE             p-value threshold for Fisher's exact test; default [0.05]
    --distance DISTANCE           Maximum distance between intervals allowed to be merged; default [100]
    --outfile OUTFILE             Output file prefix; default [diffER]
```

#### Example command with genome build
```
python diffER.py \
    --genome_build hg38 \
    --group_A_beds ./data/sample_A*.bed \
    --group_B_beds ./data/sample_B*.bed
```

#### Example command with genome file 
```
python diffER.py \
    --genome_file ./data/genome_file \
    --group_A_beds ./data/sample_A*.bed \
    --group_B_beds ./data/sample_B*.bed 
```

#### Output files
For each of the two groups, there are two output bed files containing enriched regions.   
```
  - diffER_group_A_enriched_regions.bed 
  - diffER_group_B_enriched_regions.bed
```

### 4. Contact  
Rashedul Islam, PhD (rashedul.gen@gmail.com)