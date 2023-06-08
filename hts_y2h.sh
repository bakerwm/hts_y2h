#!/usr/bin/bash

# Name: hts_y2h, high throughput sequencing for Yeast 2 Hybrid
# Date: 2023-06-07
# Author: Wang Ming
# version: 1.0

# Purpose
# extract AD/BD pairs from high through sequencing data
# the following graph explain the library structure
# ---(AD/BD)----[vec-1][attL][vec-2]----(BD/AD)---
#
# The sequences (5'->3'):
# 1. vec-1: TAGAACCCAGCTTTCTTGTACAAAGTGGTGAGCTTGGGCCCGTTTAAAC
# 2. vec-2: GATTATAAGGATGACGACGATAAAGGGCACTCGAGATATCTAGACCCAGCTTTCTTGTACAAAGTGGTGAGCTC (rev-comp)
# 3. attL: GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG (sens)
# 4. attL_rc: CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC (rev-comp)

# How-To
# 1. fetch attL: extract all reads contain full-length 33-bp-attL or its rev-comp sequence
# 2. remove attL: remove attL and downstream sequence from each 3' end of read (minlen=69)
# 3. remove vector: remove AD/BD vector sequence from 3' end of each read (minlen=20)
# 4. pairing read12: output format, ID,read1_seq,read2_seq

## required tools
## cutadapt, hisat2, samtools 
## python modules: pytabix

################################################################################
## Global variables
CPU=16
HG38_IDX="/data/yulab/hiseq001/data/genome/hg38/hisat2_index/hg38"
GENE_BED="/data/yulab/hiseq001/user/wangming/hts_y2h/data/db/Homo_sapiens.GRCh38.106.gene.bed.gz"
SRC_DIR=$(dirname $(realpath -s $0)) # path to current scirpt
ANNO_PY="${SRC_DIR}/anno_bed.py"
# [[ ! -f ${ANNO_PY} ]] && echo "script not found: ${ANNO_PY}" && exit 1
# [[ ! -f ${GENE_BED} ]] && echo "hg38 gene bed not found: ${GENE_BED}" && exit 1
# [[ ! -f ${HG38_IDX}.1.ht2 ]] && echo "hisat2 index not found: ${HG38_IDX}" && exit 1

################################################################################
## modules
## SE mode
function has_command() {
    if command -v $1 >/dev/null 2>&1
    then 
        echo "yes"
    else
        echo "no"
    fi
}
export -f has_command


function has_pymodule() {
    if python -c "import $1" &>/dev/null
    then 
        echo "yes"
    else
        echo "no"
    fi
}
export -f has_pymodule


function check_commands() {
    >&2 echo "------------------------------"
    >&2 echo "Required tools:"
    for cmd in cutadapt hisat2 samtools 
    do 
        ss=$(has_command ${cmd})
        >&2 printf "%4s : %-12s : %s\n" ${ss} ${cmd} $(which ${cmd})
        echo ${ss}
    done 
    # python module
    pp=$(has_pymodule "tabix")
    >&2 printf "%4s : %-12s : %s\n" ${pp} "tabix" "python module 'pytabix'"
    echo ${pp}
    # genome files
    [[ -f ${ANNO_PY} ]] && f1="yes" || f1="no"
    [[ -f ${GENE_BED} ]] && f2="yes" || f2="no"
    [[ -f ${HG38_IDX}.1.ht2 ]] && f3="yes" || f3="no"    
    >&2 printf "%4s : %-12s : %s\n" ${f1} "anno_py" ${ANNO_PY}
    >&2 printf "%4s : %-12s : %s\n" ${f2} "gene_bed" "${GENE_BED}"
    >&2 printf "%4s : %-12s : %s\n" ${f3} "hisat2_index" "${HG38_IDX}"
    >&2 echo "------------------------------"
    echo ${f1} ${f2} ${f3}
}
export -f check_commands


function trim_attL() {
    local out_dir=$1
    local fq_in=$2
    local attL=$3 # sens, anti
    # check output 
    [[ ${fq_in} = *_2.fq.gz ]] && rd=2 || rd=1
    if [[ ${attL} = anti ]]
    then
        as="a"
        ad="GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG"
    else
        as="s"
        ad="CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC"
    fi # anti/sens
    local suffix="_${rd}${as}"
    local fname=$(basename ${fq_in/.fq.gz})
    local fq_out="${out_dir}/${fname}${suffix}.fq.gz"
    local log="${out_dir}/${fname}${suffix}.cutadapt.log"
    [[ -f ${fq_out} ]] && echo "... file exists: ${fq_out}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    cutadapt -j ${CPU} --action trim --discard-untrimmed -O 33 -e 0.1 -m 69 -n 1 -a ${ad} -o ${fq_out} ${fq_in} > ${log}
    # echo ${fq_out} # file name
}
export -f trim_attL


function trim_vec() {
    local out_dir=$1
    local fq_in=$2
    local fq_out="${out_dir}/$(basename ${fq_in})"
    local log=${fq_out/.fq.gz/.cutadapt.log}
    local vec_1="TAGAACCCAGCTTTCTTGTACAAAGTGGTGAGCTTGGGCCCGTTTAAAC" # AD
    local vec_2="GATTATAAGGATGACGACGATAAAGGGCACTCGAGATATCTAGACCCAGCTTTCTTGTACAAAGTGGTGAGCTC" # BD
    [[ -f ${fq_out} ]] && echo "... file exists: ${fq_out}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    cutadapt -j ${CPU} --discard-untrimmed -O 12 -e 0.1 -m 20 -n 2 -a ${vec_1} -a ${vec_2} -o ${fq_out} ${fq_in} > ${log}
    # echo ${fq_out} # file name
}
export -f trim_vec


# Annotate fq, by gene_name/gene_id/NULL ...
# human library
function align_hg38() {
    local out_dir=$1
    local fq=$2 # fasta
    ## Global variable: HG38_IDX
    local fname=$(basename ${fq/.fq.gz})
    local bam="${out_dir}/${fname}.bam"
    local bed="${bam/.bam/.bed}"
    local log="${out_dir}/${fname}.hisat2.log"
    # [[ -f ${bam} ]] && echo "file exists: ${bam}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    [[ ! -f ${bam} ]] && \
        hisat2 -p ${CPU} --very-sensitive --add-chrname -x ${HG38_IDX} -U ${fq} 2> ${log} | \
        samtools view -Sub -F 0x4 -F 2048 -q 10 - | \
        samtools sort -@ ${CPU} -o ${bam} - && \
        samtools index ${bam}
    # convert bam to bed
    [[ ! -f ${bed} ]] && \
        bedtools bamtobed -i ${bam} > ${bed} #
    # echo ${bed}
}
export -f align_hg38


function anno_bed() {
    local out_dir=$1
    local in_bed=$2 # query
    local bname=$(basename ${in_bed/.bed/.anno.bed})
    local out_bed="${out_dir}/${bname}" # output 
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    ## Global variables: ANNO_PY, GENE_BED
    # GRCh38, all genes
    # local GENE_BED="/data/yulab/hiseq001/user/wangming/hts_y2h/data/db/Homo_sapiens.GRCh38.106.gene.bed.gz"
    # local ANNO_PY="/data/yulab/hiseq001/user/wangming/hts_y2h/scripts/anno_bed.py"
    [[ ! -f ${out_bed} ]] && python ${ANNO_PY} ${GENE_BED} ${in_bed} ${out_bed}
}
export -f anno_bed


# annotate fastq, from 
function anno_fq() {
    local out_dir=$1
    local fq_in=$2
    local fname="$(basename ${fq_in/.fq.gz})"
    local bed="${out_dir}/${fname}.bed"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # 1. map to hg38 genome
    [[ ! -f ${bed} ]] && align_hg38 ${out_dir} ${fq_in}
    # 2. get annotation (gene_name)
    anno_bed ${out_dir} ${bed} 
    # 3. log
    echo "... save to file: ${bed}"
}
export -f anno_fq


# output common reads (by_id) from two fastq files
# save to new directory
function merge_fq_by_id() {
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # file-1:
    out1="${out_dir}/$(basename ${fq1/.fq/.fa})"
    [[ -f ${out1} ]] && echo "... file exists: ${out1}" || \
        bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print ">"$name"\n"a[$name]"\n"}}' ${fq1} ${fq2} | gzip > ${out1}
    # file-2:
    out2="${out_dir}/$(basename ${fq2/.fq/.fa})"
    [[ -f ${out2} ]] && echo "... file exists: ${out2}" || \
        bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print ">"$name"\n"a[$name]"\n"}}' ${fq2} ${fq1} | gzip > ${out2}
}
export -f merge_fq_by_id


# input: bed6+pos+gene_name
# ouput: read_id,gene_name,gene_name
function merge_bed_by_id() {
    local out_dir=$1
    local bed1=$2
    local bed2=$3
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # file-1:
    out1="${out_dir}/$(basename ${bed1/.bed/.read12.bed})"
    [[ ! -f ${out1} ]] && \
        awk 'FNR==NR {a[$4]=$8; next} {if($4 in a){print $4,a[$4],$8}}' ${bed1} ${bed2} > ${out1}
    # # file-2:
    # out1="${out_dir}/$(basename ${bed1/.bed/.read12.bed})"
    # [[ ! -f ${out2} ]] && \
    #     awk 'FNR==NR {a[$4]=$8; next} {if($4 in a){print $4,$8,a[$4]}}' ${bed2} ${bed1} > ${out2}
    ## log
    echo "... $(wc -l ${out1})"
}
export -f merge_bed_by_id


# wrap all steps
function main() {
    local out_dir=$1
    local fq1=$2
    local fq2=$3

    # 0. check input files
    [[ ! -f ${fq1} || ! ${fq1} = *fq.gz ]] && echo "fq1, not .fq.gz: ${fq1}" && return 1
    [[ ! -f ${fq2} || ! ${fq2} = *fq.gz ]] && echo "fq2, not .fq.gz: ${fq2}" && return 1

    # 1. extract/trim attL
    echo "[1/4] trimming attL"
    out_dir1="${out_dir}/1.trim_attL"
    trim_attL ${out_dir1} ${fq1} sens
    trim_attL ${out_dir1} ${fq1} anti
    trim_attL ${out_dir1} ${fq2} sens
    trim_attL ${out_dir1} ${fq2} anti

    # 2. trim vector
    echo "[2/4] trimming vector"
    out_dir2="${out_dir}/2.trim_vector"
    for fq in ${out_dir1}/*gz 
    do 
        trim_vec ${out_dir2} ${fq}
    done

    # 3. annotate reads
    echo "[3/4] annotating reads"
    out_dir3="${out_dir}/3.annotate_reads"
    for c in ${out_dir2}/*fq.gz 
    do 
        anno_fq ${out_dir3} ${c}
    done

    # 4. paring read12
    echo "[4/4] find read1,2 paires"
    out_dir4="${out_dir}/4.pairing_reads"
    r1_name=$(basename ${fq1/.fq.gz})
    r2_name=$(basename ${fq2/.fq.gz})
    bed1s="${out_dir3}/${r1_name}_1s.anno.bed"
    bed1a="${out_dir3}/${r1_name}_1a.anno.bed"
    bed2s="${out_dir3}/${r2_name}_2s.anno.bed"
    bed2a="${out_dir3}/${r2_name}_2a.anno.bed"
    ## bed format:
    merge_bed_by_id ${out_dir4} ${bed1s} ${bed2a}
    merge_bed_by_id ${out_dir4} ${bed1a} ${bed2s}

    # 5. finish
    echo "Finish!"
}
export -f main

[[ $# -lt 3 ]] && echo "Usage: hts_y2h.sh <out_dir> <fq1> <fq2>" && exit 1
status=$(check_commands)
[[ ${status} = *no* ]] && echo "!!! error, check above missing files" && exit 1
main $1 $2 $3
