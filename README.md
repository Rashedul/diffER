# diffER
Identify differentially enriched regions from ChIP-seq data

#### build database
Example of `--genome` is `hg38`, `hg19`.  You can check available genomes `genomepy genomes`

```
build.py --genome <genome_build> --window_size <window size>
```

#### build database

```
build.py --genome --window_size 
```

#### run diffER

```
diffER.py --group_A <bed files> --group_B <bed files> --p_value <p-value cutoff> 
```
