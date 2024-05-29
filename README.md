# diffER (under development)
Identify differentially enriched regions from ChIP-seq peaks. 

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
Example of `--genome_build` is `hg38`, `hg19`, `mm9`, `mm10` etc.  You can check available genomes in `genomepy`. Otherwise, provide a genome file of your interest. Example is added [here](genome_file). 

```
python diffER.py --genome_file genome_file --group_A_beds ./data/sample_A*.bed  --group_B_beds ./data/sample_B*.bed --window_size 200
```
