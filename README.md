# diffER (under development)
Identify differentially enriched regions from ChIP-seq peaks. 

### 1. Running diffER 

#### Requirements 
- python >= 3.8

#### Installation

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

### 2. build database
Example of `--genome_build` is `hg38`, `hg19`, `mm9`, `mm10` etc.  You can check available genomes in `genomepy`. Otherwise, provide a genome file of your interest. Example is added [here](genome_file). 

```
python build.py --genome_build <genome build> | --genome_file <genome file> --window_size <window size>
```

### 3. run diffER

```
python diffER.py --group_A <bed files> --group_B <bed files> --genome <genome_build> --window_size <window size> --p_value <p-value cutoff> 

python t3.py --primary_bed genome_file_windows.bed --group_A ./data/*A*.bed --group_B ./data/*B*.bed
```
