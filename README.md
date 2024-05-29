# diffER (under development)
Identify differentially Enriched Regions (diffER) from ChIP-seq peaks. 

### 1. Requirements 
- python >= 3.8
- Linux

### 2. Installation

 - Cloning the repository

```
git clone https://github.com/Rashedul/diffER.git
```

 - Creation of Python Virtual Environment

```
cd diffER
python -m venv environment
source ./environment/bin/activate
```

 - Installation of Python Packages

```
pip install -r requirements.txt
```

### 3. Run diffER

- Provide either a genome build or genome file 
- Example of `--genome_build` is `hg38`, `hg19`, `mm9`, `mm10` etc.  You can check available genomes in `genomepy`. 
- Provide a genome file of your interest. Example of genome file is provided [here](./data/genome_file). 


#### Find option details
```
python diffER.py -h
```

#### Example command with genome build
```
python diffER.py \
    --genome_build hg38 \
    --group_A_beds ./data/sample_A*.bed \
    --group_B_beds ./data/sample_B*.bed \
    --window_size 50
```

#### Example command with genome file 
```
python diffER.py \
    --genome_file ./data/genome_file \
    --group_A_beds ./data/sample_A*.bed \
    --group_B_beds ./data/sample_B*.bed \
    --window_size 50
```

#### Evaluate output in a heatmap
```
python diffER_heatmap.py \
	--outfile diffER_out.bed \
	--group_A_beds ./data/sample_A*.bed \
    --group_B_beds ./data/sample_B*.bed \
    --filename diffER_heatmap 
``` 

### 4. Notes
- It's recommended to evaluate the output carefylly, particularly for samples with higher sample-to-sample variations. 
- Check your results at different p-values and/or odd ratios. 
- Required number of samples per group is at least 3.

### 5. Citation 