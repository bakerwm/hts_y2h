# hts_y2h

for high throughput sequencing of yeast two-hybrid experiment

## Purpose
extract read1,read2 pairs from high through sequencing data
the following graph explain the library structure
`---(AD/BD)----[vec-1][attL][vec-2]----(BD/AD)---`

The sequences (5'->3'):
1. vec-1: `TAGAACCCAGCTTTCTTGTACAAAGTGGTGAGCTTGGGCCCGTTTAAAC`
2. vec-2: `GATTATAAGGATGACGACGATAAAGGGCACTCGAGATATCTAGACCCAGCTTTCTTGTACAAAGTGGTGAGCTC` (rev-comp)
3. attL: `GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG` (sens)
4. attL_rc: `CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC` (rev-comp)

## How-To

1. fetch attL: extract all reads contain full-length 33-bp-attL or its rev-comp sequence
2. remove attL: remove attL and downstream sequence from each 3' end of read (minlen=69)
3. remove vector: remove AD/BD vector sequence from 3' end of each read (minlen=20)
4. annotate the AD/BD sequence by gene_name
5. pairing read12: output format, ID,gene_1,gene_2

## Getting Started 

Run the program: 
```
$ bash hts_y2h.sh results demo_1.fq.gz demo_2.fq.gz
```

## Prerequisites 

+ `python` >=3.9.13
+ `cutadapt` >=v4.1 
+ `hisat2` >=2.2.1 
+ `samtools` >=v1.16.1 
+ `pytabix` 

## Installation

No need to install the program actually. You need just update requirements as follow:

1. Install Prerequisites
    + Using `conda` for example
    ```
    $ conda install -c bioconda cutadapt hisat2 samtools pytabix
    ```
    
2. Clone the repo to your local machine

> update the two files from `hts_y2h.sh` by realpath on your computer. `HG38_IDX` and `GENE_BED`;

```
$ git clone https://github.com/bakerwm/hts_y2h
$ cd hts_y2h
$ nano hts_y2h.sh 
```

Update the following two files
```
# file: hts_y2h.sh
# line 32-33
...
HG38_IDX="/data/yulab/hiseq001/data/genome/hg38/hisat2_index/hg38"
GENE_BED="/data/yulab/hiseq001/user/wangming/hts_y2h/data/db/Homo_sapiens.GRCh38.106.gene.bed.gz"
...
```

3. Check the program

```
$ bash hts_y2h.sh a b c
------------------------------
Required tools:
 yes : cutadapt     : /data/yulab/hiseq001/miniconda3/envs/hiseq/bin/cutadapt
 yes : hisat2       : /data/yulab/hiseq001/miniconda3/envs/hiseq/bin/hisat2
 yes : samtools     : /data/yulab/hiseq001/miniconda3/envs/hiseq/bin/samtools
 yes : tabix        : python module 'pytabix'
 yes : anno_py      : /data/yulab/hiseq001/user/wangming/hts_y2h/scripts/anno_bed.py
 yes : gene_bed     : /data/yulab/hiseq001/user/wangming/hts_y2h/data/db/Homo_sapiens.GRCh38.106.gene.bed.gz
 yes : hisat2_index : /data/yulab/hiseq001/data/genome/hg38/hisat2_index/hg38
------------------------------
fq1, not .fq.gz: b
```


## License
Distributed under the MIT License. See LICENSE.txt for more information.


## Contact
Ming Wang - @github.com/bakerwm - wangm08(AT)hotmail.com

## Acknowledgments

+ ...
