# diffER
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
python3 -m venv environment
source ./environment/bin/activate
```

 - Installation of Python Packages

```
pip install -r requirements.txt
```

### 2. build database
Example of `--genome` is `hg38`, `hg19`, `mm9`, `mm10` etc.  You can check available genomes in `genomepy`. Otherwise, provide a genome file of your interest. Example is added [here](genome_file). 

```
build.py --genome_build <genome build> | --genome_file <genome file> --window_size <window size>
```

### 3. run diffER

```
diffER.py --group_A <bed files> --group_B <bed files> --genome <genome_build> --window_size <window size> --p_value <p-value cutoff> 
```
