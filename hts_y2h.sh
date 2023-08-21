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

# Change log
## vector-1: TAGAACCCAGCTTTC 
## vector-2: GATTATAAGGATGAC
## 
# level-1: attL yes|no
# level-2: read1, vector-1 yes|no
# level-3: read2, vector-2 yes|no
#
# directory structure:
# tn5
# ├── attL+
# │   └── v1 (r1v1+)
# │       └── v2 (r2v2-)
# └── attL-
#     ├── v1 (r1v1+)
#     │   └── v2 (r2v2-)
#     └── v1 (r1v1-)
#         └── v2 (r2v2-)
#
# parameters:
# out_dir, fq1, fq2, 3-adapter, <1|2>
#
# out_dir
# ├── suffix
# │   ├── *r1v1+_1.fq.gz
# │   ├── *r1v1+_2.fq.gz
# │   ├── *r1v1-_1.fq.gz
# │   ├── *r1v1-_2.fq.gz
# │   └── *r1v1.cutadapt.log
# └── suffix-
#     ├── fq.gz
#     └── log
#
# 3. annotate read1, read2 (gene)
#

################################################################################
## Global variables
SRC_DIR=$(dirname $(realpath -s $0))
[[ -z ${N_CPU} ]] && N_CPU=8
SRC_DIR=$(dirname $(realpath -s $0)) # path to current scirpt
PARSE_GTF="${SRC_DIR}/parse_gtf.py"
ANNO_PY="${SRC_DIR}/anno_bed.py"
ANNO_FQ="${SRC_DIR}/anno_fq.sh"
ANNO_FQ_PAIR="${SRC_DIR}/anno_fq_pair.sh"
## Change the following variables accroding to your experiment ##
HG38_IDX="/data/yulab/hiseq001/data/genome/hg38/hisat2_index/hg38"
GENE_GTF="/data/biodata/ensembl/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"
## GLOBAL VECTOR
TN5="CTGTCTCTTATACA"    # tn5 adapters
VEC1="TAGAACCCAGCTTTC"  # vector-1
VEC2="GATTATAAGGATGAC"  # vector-2
ATTLs="GTGCCAGGGCGTGCCCTTGAGTTCTCTCAGTTG" # sens
ATTLa="CAACTGAGAGAACTCAAGGGCACGCCCTGGCAC" # anti
################################################################################

## modules
function has_command() {
    if command -v $1 >/dev/null 2>&1
    then 
        echo "ok"
    else
        echo "failed"
    fi
}
export -f has_command


function has_pymodule() {
    if python -c "import $1" &>/dev/null
    then 
        echo "ok"
    else
        echo "failed"
    fi
}
export -f has_pymodule


function check_commands() {
    cmd_tag=()
    >&2 echo "--------------------------------------------------------------------------------"
    >&2 echo "Required tools:"
    for cmd in cutadapt hisat2 samtools tabix
    do 
        ss=$(has_command ${cmd})
        >&2 printf "%6s : %-12s : %s\n" ${ss} ${cmd} $(which ${cmd})
        cmd_tag+=("${ss}")
        echo ${ss}
    done 
    # python module
    pp=$(has_pymodule "tabix")
    >&2 printf "%6s : %-12s : %s\n" ${pp} "tabix" "python module 'pytabix'"
    echo ${pp}
    # genome files
    [[ -f ${HG38_IDX}.1.ht2 ]] && f1="ok" || f1="failed"
    [[ -f ${GENE_GTF} ]] && f2="ok" || f2="failed"
    [[ -f ${PARSE_GTF} ]] && f3="ok" || f3="failed"
    [[ -f ${ANNO_PY} ]] && f4="ok" || f4="failed"
    [[ -f ${ANNO_FQ} ]] && f5="ok" || f5="failed"
    [[ -f ${ANNO_FQ_PAIR} ]] && f6="ok" || f6="failed"
    >&2 printf "%6s : %-12s : %s\n" ${f1} "hisat2_index" "${HG38_IDX}"
    >&2 printf "%6s : %-12s : %s\n" ${f2} "gene_gtf" "${GENE_GTF}"
    >&2 printf "%6s : %-12s : %s\n" ${f3} "parse_gtf" ${PARSE_GTF}
    >&2 printf "%6s : %-12s : %s\n" ${f4} "anno_py" ${ANNO_PY}
    >&2 printf "%6s : %-12s : %s\n" ${f5} "anno_fq" ${ANNO_FQ}
    >&2 printf "%6s : %-12s : %s\n" ${f6} "anno_fq_pair" ${ANNO_FQ_PAIR}
    >&2 echo "--------------------------------------------------------------------------------"
    # chk_tag+=("${f1} ${f2} ${f3} ${f4} ${f5} ${f6}")
    # echo ${f1} ${f2} ${f3} ${f4} ${f5} ${f6}
    echo ${cmd_tag[@]} ${f1} ${f2} ${f3} ${f4} ${f5} ${f6}
}
export -f check_commands
################################################################################


# convert gtf to bed, and bgzip compress for tabix
# para: <out_dir> <GENE_GTF>
function gtf2bed() {
    [[ $# -lt 2 ]] && echo "Usage: gtf2bed <out_dir> <gtf>" && return 1
    local out_dir=$1
    local gtf=$2 # Ensembl gtf
    local gname=$(basename ${gtf%.gz})
    gname=${gname%.gtf} # trim suffix
    local bed="${out_dir}/${gname}.bed.gz" # bgzipped compressed
    if [[ ! -f ${bed} ]]
    then
        [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
        local bed_dir="${out_dir}/gene_files"
        local bed_raw="${bed_dir}/${gname}.gene.f-1.pc.g100.bed"
        local log_raw="${out_dir}/${gname}.gene.f-1.pc.g100.log"
        [[ ! -d ${bed_dir} ]] && mkdir -p ${bed_dir}
        python ${PARSE_GTF} -o ${bed_dir} -f gene -g protein_coding -gs 100 -fs -1 -i ${gtf} 1> ${log_raw}
        # gzipped, tabix
        bgzip -c ${bed_raw} > ${bed}
        tabix -p bed ${bed}
    fi
    [[ ! -f ${bed} ]] && echo "bed file not exists" && return 1
    echo ${bed}
}
export -f gtf2bed


# Prepare gene_bed file, from GENE_GTF
# output: <out_dir> <GENE_GTF>
function get_gene_bed() {
    local bed=$(gtf2bed $1 $2)
    [[ -f ${bed} ]] && echo ${bed}
}
export -f get_gene_bed


# # wrap all steps
# function main() {
#     local out_dir=$1
#     local fq1=$2
#     local fq2=$3

#     # 0. check input files
#     [[ ! -f ${fq1} || ! ${fq1} = *fq.gz ]] && echo "fq1, not .fq.gz: ${fq1}" && return 1
#     [[ ! -f ${fq2} || ! ${fq2} = *fq.gz ]] && echo "fq2, not .fq.gz: ${fq2}" && return 1

#     # 1. extract/trim attL
#     echo "[1/4] trimming attL"
#     out_dir1="${out_dir}/1.trim_attL"
#     trim_attL ${out_dir1} ${fq1} sens
#     trim_attL ${out_dir1} ${fq1} anti
#     trim_attL ${out_dir1} ${fq2} sens
#     trim_attL ${out_dir1} ${fq2} anti

#     # 2. trim vector
#     echo "[2/4] trimming vector"
#     out_dir2="${out_dir}/2.trim_vector"
#     for fq in ${out_dir1}/*gz 
#     do 
#         trim_vec ${out_dir2} ${fq}
#     done

#     # 3. annotate reads
#     echo "[3/4] annotating reads"
#     out_dir3="${out_dir}/3.annotate_reads"
#     for c in ${out_dir2}/*fq.gz 
#     do 
#         anno_fq ${out_dir3} ${c}
#     done

#     # 4. paring read12
#     echo "[4/4] find read1,2 paires"
#     out_dir4="${out_dir}/4.pairing_reads"
#     r1_name=$(basename ${fq1/.fq.gz})
#     r2_name=$(basename ${fq2/.fq.gz})
#     bed1s="${out_dir3}/${r1_name}_1s.anno.bed"
#     bed1a="${out_dir3}/${r1_name}_1a.anno.bed"
#     bed2s="${out_dir3}/${r2_name}_2s.anno.bed"
#     bed2a="${out_dir3}/${r2_name}_2a.anno.bed"
#     ## bed format:
#     merge_bed_by_id ${out_dir4} ${bed1s} ${bed2a}
#     merge_bed_by_id ${out_dir4} ${bed1a} ${bed2s}

#     # 5. finish
#     echo "Finish!"
# }
# export -f main

################################################################################
## cutadapt
## vector-1: TAGAACCCAGCTTTC 
## vector-2: GATTATAAGGATGAC
# prefix
function fx_prefix() {
    local fq=$1
    basename ${fq} | sed -Ee 's/(_[0-9]+)?.f(ast)?[aq](.gz)?$//i'
}
export -f fx_prefix


function revcomp() {
    local sequence="$1"
    local complement=$(echo "${sequence}" | tr 'ATCG' 'TAGC')
    local reverse_complement=$(echo "${complement}" | rev)
    echo "${reverse_complement}"
}
export -f revcomp


# check vector name: v1,v2
function get_ad_name() {
    if [[ $1 == ${TN5} ]] 
    then 
        echo "tn5"
    elif [[ $1 == ${VEC1} ]] 
    then
        echo "v1"
    elif [[ $1 == ${VEC2} ]]
    then 
        echo "v2"
    elif [[ $1 == ${ATTLs} || $1 == ${ATTLa} ]]
    then
        echo "attL"
    else
        echo "null"
        return 1
    fi
}
export -f get_ad_name


# get the suffix for read1/2, vector1/2
# para: <adapter> <1|2>
function get_suffix() {
    [[ $# -lt 2 ]] && echo "rxvx" && return 1
    vname=$(get_ad_name $1) # v1/2
    echo "r$2${vname}" # r1v1, r1v2, ...
}
export -f get_suffix


# para: <log> <keyword:str>
function parse_by_keyword() {
    local log=$1
    local key="$2" # by "", if spaces exists
    local num=$(grep "${key}" ${log} | sed -Ee 's/\([0-9.%]+\)//' -e 's/,//g' | awk '{print $NF}')
    local n_rec=$(grep -c "${key}" ${log})
    if [[ ${n_rec} -eq 0 ]] 
    then
        # return 0 if not found
        echo 0 # default
    elif [[ ${n_rec} -eq 1 ]]
    then        
        echo ${num}
    elif [[ ${n_rec} -gt 1 ]]
    then
        # the biggest value
        echo ${num} | xargs -n 1 | sort -k1nr | head -n 1
    fi
    # return max if multi matched
}
export -f parse_by_keyword


# para: <log>
function parse_cutadapt_log() {
    local log=$1 # cutadapt
    [[ ! ${log} = *cutadapt.log ]] && echo "invalid log name" && return 1
    [[ ! -f ${log} ]] && echo "log not exists" && return 1
    local fname=$(basename ${log%.cutadapt.log})
    local ad=$(basename $(dirname ${log})) # ad
    # header
    echo name,ad,total,ad+,ad-,output,discard
    local total=$(parse_by_keyword ${log} "Total read pairs processed")
    local ad_plus=$(parse_by_keyword ${log} "with adapter")
    # local ad_minus=$(parse_by_keyword ${log} "as untrimmed")
    # local ad_minus=$(echo ${total} ${ad_plus} | awk '{print $1-$2}')
    local ad_minus=$(echo ${total} ${ad_plus} | awk '{printf "%d",$1-$2}') #
    local output=$(parse_by_keyword ${log} "passing filters")
    local discard=$(parse_by_keyword ${log} "too short")
    ## update ad_minus
    # if [[ ${ad_minus} -eq 0 ]] 
    # then 
    #     ad_minus=$(echo ${total} ${ad_plus} | awk '{printf "%d",t-p}') #
    # fi
    # report
    echo ${fname},${ad},${total},${ad_plus},${ad_minus},${output},${discard}
}
export -f parse_cutadapt_log


# trim 3' adapter
# para: <out_dir> <fq1> <fq2> <ad> <1|2>
function trim_3ad() {
    [[ $# -lt 5 ]] && echo "Usage: trim_3ad <out_dir> <fq1> <fq2> <ad> [extra]" && return 1
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    local ad=$4         # adapter
    local ad_in_read=$5 # adapter in read1/2, [1|2], 3=read1 and read2
    # local suffix=$6
    # local para_extra=$7
    # check parameters
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not exists: ${fq2}" && return 1
    # [[ ! ${fq1} = *1.fq.gz ]] && "fq1 invalid name [_1.fq.gz]: ${fq1}" && return 1
    # [[ ! ${fq2} = *2.fq.gz ]] && "fq2 invalid name [_2.fq.gz]: ${fq2}" && return 1
    [[ ! ${ad_in_read} = [123] ]] && "ad_in_read invalid [1|2|3]: ${ad_in_read}" && return 1 
    # check suffix
    local fname=$(fx_prefix ${fq1}) # prefix, trim _[12].fq.gz
    local vname=$(get_ad_name ${ad}) # v1,v2,attL
    local suffix=$(get_suffix ${ad} ${ad_in_read}) # r1v1, r2v2, ...
    # prepare files
    local out_sub="${out_dir}/${vname}"
    local prefix="${out_dir}/${vname}/${fname}.${suffix}"
    [[ ! -d ${out_sub} ]] && mkdir -p ${out_sub}
    # ad+
    local fq1_plus="${out_sub}/${fname}.${suffix}+_1.fq.gz"
    local fq2_plus="${out_sub}/${fname}.${suffix}+_2.fq.gz"
    # ad-
    local fq1_minus="${out_sub}/${fname}.${suffix}-_1.fq.gz"
    local fq2_minus="${out_sub}/${fname}.${suffix}-_2.fq.gz"
    # log,cmd
    local log="${out_sub}/${fname}.${suffix}.cutadapt.log"
    local cmd="${out_sub}/${fname}.${suffix}.cmd.sh"
    # command, para
    ## adapter:BEGIN
    if [[ ${ad_in_read} == 1 ]]
    then
        para_ad="-a ${ad}" # for vector at 3' end of read1
    elif [[ ${ad_in_read} == 2 ]]
    then 
        para_ad="-A ${ad}" # for vector at 3' end of read2
    elif [[ ${ad_in_read} == 3 ]]
    then 
        # para_ad="-a ${ad} -A ${ad}" # for attL at 3' of read1 or read2
        local ad_rc=$(revcomp ${ad})
        para_ad="-a ${ad} -a ${ad_rc} --times 4" # for attL at 3' of read1 or read2
    else
        para_ad=""
    fi
    ## adapter:END
    [[ -f ${fq1_plus} && -f ${fq2_plus} ]] && run_tag="# " || run_tag="" #
    echo "${run_tag} cutadapt -j ${N_CPU} -m 15 --action trim ${para_ad} \
            --untrimmed-output ${fq1_minus} \
            --untrimmed-paired-output ${fq2_minus} \
            -o ${fq1_plus} -p ${fq2_plus} ${fq1} ${fq2} > ${log}" > ${cmd}
    # run
    bash ${cmd}
    # wrap log file
    local stat="${out_sub}/${fname}.${suffix}.stat"
    parse_cutadapt_log ${log} > ${stat}
    # return prefix
    echo ${prefix}
}
export -f trim_3ad


# trim Tn5 at 3' adapter
# para: <out_dir> <fq1> <fq2> <ad>
function trim_3ad_tn5() {
    [[ $# -lt 4 ]] && echo "Usage: trim_3ad_tn5 <out_dir> <fq1> <fq2> <ad>" && return 1
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    local ad=$4 # adapter
    # check parameters
    [[ ! -f ${fq1} ]] && echo "fq1 not exists: ${fq1}" && return 1
    [[ ! -f ${fq2} ]] && echo "fq2 not exists: ${fq2}" && return 1
    # check suffix
    local fname=$(fx_prefix ${fq1}) # prefix, trim _[12].fq.gz
    local vname=$(get_ad_name ${ad}) # tn5,v1,v2,attL
    # prepare files
    local out_sub="${out_dir}/${vname}"
    local prefix="${out_dir}/${vname}/${fname}"
    [[ ! -d ${out_sub} ]] && mkdir -p ${out_sub}
    local fq1_out="${prefix}_1.fq.gz"
    local fq2_out="${prefix}_2.fq.gz"
    # log,cmd
    local log="${prefix}.cutadapt.log"
    local cmd="${prefix}.cmd.sh"
    # command, para
    local para_ad="-a ${ad} -A ${ad}" # both read1 and read2
    [[ -f ${fq1_out} && -f ${fq2_out} ]] && run_tag="# " || run_tag="" #
    echo "${run_tag} cutadapt -j ${N_CPU} -m 15 --action trim ${para_ad} \
            -o ${fq1_out} -p ${fq2_out} ${fq1} ${fq2} > ${log}" > ${cmd}
    # run
    bash ${cmd}
    # wrap log file
    local stat="${out_sub}/${fname}.stat"
    parse_cutadapt_log ${log} > ${stat}
    # return prefix
    echo ${prefix}
}
export -f trim_3ad_tn5


## format stat
## fmt1: +(1234, 100.0%) -(1234, 100.0%)
## fmt2: output(1234, 100.0%) discard(1234, 100.0%)
#
# input:
# name,ad,total,ad+,ad-,output,discard
# demo,tn5,10000,7285,2715,9784,216
function fmt_stat() {
    local stat=$1 # file
    local fmt=$2  # format: 1,2
    local err="0,0.0% 0,0.0%"
    [[ ! -f ${stat} ]] && echo "${err}" && return 1 # stat not exists
    local s=$(grep -v discard $1 | sed 's/,/ /g') # total=0
    local total=$(echo ${s} | awk '{print $3}')  # total=0
    # echo 3. ${s}
    [[ ${total} -eq 0 ]] && echo "${err}" && return 1 # total=0
    if [[ ${fmt} -eq 1 ]] 
    then 
        # ad+, ad-
        echo ${s} | awk '{printf "%d,%.1f%% %d,%.1f%%\n", $4, $4/$3*100, $5, $5/$3*100}'
    elif [[ ${fmt} -eq 2 ]] 
    then
        # output, discard
        echo ${s} | awk '{printf "%d,%.1f%% %d,%.1f%%\n", $6, $6/$3*100, $7, $7/$3*100}'
    else 
        echo "${err}" 
    fi
}
export -f fmt_stat


################################################################################
# Example of wrapped stat output
# 
# input (1234, 100.0%)
# ├── tn5 (1234, 100.0%)
# │   ├── attL(+) (1234, 100.0%)
# │   │   ├── read1:vector-1(+) (1234, 100.0%)
# │   │   │   ├── read2:vector-2(+) (1234, 100.0%)
# │   │   │   └── read2:vector-2(-) (1234, 100.0%)
# │   │   └── read1:vector-1(-) (1234, 100.0%)
# │   │       ├── read2:vector-2(+) (1234, 100.0%)
# │   │       └── read2:vector-2(-) (1234, 100.0%)
# │   └── attL(-) (1234, 100.0%)
# │       ├── read1:vector-1(+) (1234, 100.0%)
# │       │   ├── read2:vector-2(+) (1234, 100.0%)
# │       │   └── read2:vector-2(-) (1234, 100.0%)
# │       └── read1:vector-1(-) (1234, 100.0%)
# │           ├── read2:vector-2(+) (1234, 100.0%)
# │           └── read2:vector-2(-) (1234, 100.0%)
# └── discard (1234, 100.0%)
################################################################################
# wrap statistics
# para: <out_dir> <fq1>
function wrap_stat() {
    local out_dir=$1
    local fname=$(fx_prefix ${fq1}) # prefix, trim _[12].fq.gz
    # 1. Tn5: output, discard
    local s1="${out_dir}/tn5/${fname}.stat" #
    local f1=($(fmt_stat ${s1} 2)) #
    # 2. attL: ad+, ad-
    local s2="${out_dir}/tn5/attL/${fname}.r3attL.stat" #
    local f2=($(fmt_stat ${s2} 1)) #
    # 3. v1: ad+, ad-
    local s3="${out_dir}/tn5/attL/v1/${fname}.r3attL+.r1v1.stat" #
    local f3=($(fmt_stat ${s3} 1)) #
    local s4="${out_dir}/tn5/attL/v1/${fname}.r3attL-.r1v1.stat" #
    local f4=($(fmt_stat ${s4} 1)) #
    # 4. v2: ad+, ad-
    local s5="${out_dir}/tn5/attL/v1/v2/${fname}.r3attL+.r1v1+.r2v2.stat" #
    local f5=($(fmt_stat ${s5} 1)) #
    local s6="${out_dir}/tn5/attL/v1/v2/${fname}.r3attL+.r1v1-.r2v2.stat" #
    local f6=($(fmt_stat ${s6} 1)) #
    local s7="${out_dir}/tn5/attL/v1/v2/${fname}.r3attL-.r1v1+.r2v2.stat" #
    local f7=($(fmt_stat ${s7} 1)) #
    local s8="${out_dir}/tn5/attL/v1/v2/${fname}.r3attL-.r1v1-.r2v2.stat" #
    local f8=($(fmt_stat ${s8} 1)) #
    ## output, formatted
    echo ${fname} #
    echo "├── Tn5 .........................(${f1[0]})"
    echo "│   ├── attL+ ...................(${f2[0]})*"
    echo "│   │   ├── read1:vector1+ ......(${f3[0]})*"  # !!! bug1
    echo "│   │   │   ├── read2:vector2+ ..(${f5[0]})*"
    echo "│   │   │   └── read2:vector2- ..(${f5[1]})"
    echo "│   │   └── read1:vector1- ......(${f3[1]})"
    echo "│   │       ├── read2:vector2+ ..(${f6[0]})*"
    echo "│   │       └── read2:vector2- ..(${f6[1]})"
    echo "│   └── attL- ...................(${f2[1]})"
    echo "│       ├── read1:vector1+ ......(${f4[0]})*"
    echo "│       │   ├── read2:vector2+ ..(${f7[0]})*"
    echo "│       │   └── read2:vector2- ..(${f7[1]})"
    echo "│       └── read1:vector1- ......(${f4[1]})"
    echo "│           ├── read2:vector2+ ..(${f8[0]})*"
    echo "│           └── read2:vector2- ..(${f8[1]})"
    echo "└── discard .....................(${f1[1]})"
    echo ""
    echo "* indicated the count of output and too-short reads,"
    echo "  it could be larger than sum of (+/-)"

}
export -f wrap_stat


# annotate fq_pair within directory
# para: <gene_bed> <fq_dir>
function anno_dir() {
    # check fq files within <fq_dir>
    local gene_bed=$1
    local fq_dir=$2
    local n_fq=$(ls ${fq_dir}/ 2>/dev/null | grep -c 1.fq.gz) #
    [[ ${n_fq} -eq 0 ]] && echo "No fq files found: ${fq_dir}" && return 1
    for fq1 in ${fq_dir}/*1.fq.gz 
    do 
        local fq2="${fq1%1.fq.gz}2.fq.gz"
        local anno_dir="${fq_dir}/anno"
        bash ${ANNO_FQ_PAIR} ${anno_dir} ${gene_bed} ${fq1} ${fq2} # anno
    done
}
export -f anno_dir


## run analysis
# para: <out_dir> <fq1> <fq2>
function hts_y2h() {
    [[ $# -lt 3 ]] && echo "Usage: hts_y2h <out_dir> <fq1> <fq2>" && return 1
    local out_dir=$1
    local fq1=$2
    local fq2=$3
    # absolute path
    out_dir=$(realpath -s ${out_dir})
    fq1=$(realpath -s ${fq1})
    fq2=$(realpath -s ${fq2})
    # prepare gene bed files
    echo "preparing gene_bed file"
    local gene_bed=$(get_gene_bed ${out_dir}/db ${GENE_GTF}) # global variables
    echo "working in read1-read2 mode"
    echo "[1/x] - remove Tn5 adapters from both read1 and read2"
    local prefix1=$(trim_3ad_tn5 ${out_dir} ${fq1} ${fq2} ${TN5})
    local tn5_fq1="${prefix1}_1.fq.gz"
    local tn5_fq2="${prefix1}_2.fq.gz"
    anno_dir ${gene_bed} "${out_dir}/tn5"
    echo 1. ${prefix1}
    echo "[2/x] - check attL in read1"
    local attL_dir="${out_dir}/tn5"
    local prefix2=$(trim_3ad ${attL_dir} ${tn5_fq1} ${tn5_fq2} ${ATTLs} 3) # 3=read1,read2 
    local attLp_fq1="${prefix2}+_1.fq.gz" # p=plus,+, m=minus,-
    local attLp_fq2="${prefix2}+_2.fq.gz" # p=plus,+, m=minus,-
    local attLm_fq1="${prefix2}-_1.fq.gz" # p=plus,+, m=minus,-
    local attLm_fq2="${prefix2}-_2.fq.gz" # p=plus,+, m=minus,-
    anno_dir ${gene_bed} ${attL_dir}/attL
    echo 2. ${prefix2}
    echo "[3/x] - find vector-1 in read1"
    local vec1_dir="${out_dir}/tn5/attL"
    local prefix31=$(trim_3ad ${vec1_dir} ${attLp_fq1} ${attLp_fq2} ${VEC1} 1)
    local Lpv1p_fq1="${prefix31}+_1.fq.gz" # p=plus,+, m=minus,-
    local Lpv1p_fq2="${prefix31}+_2.fq.gz" # p=plus,+, m=minus,-
    local Lpv1m_fq1="${prefix31}-_1.fq.gz" # p=plus,+, m=minus,-
    local Lpv1m_fq2="${prefix31}-_2.fq.gz" # p=plus,+, m=minus,-
    echo 3. ${prefix31}
    local prefix32=$(trim_3ad ${vec1_dir} ${attLm_fq1} ${attLm_fq2} ${VEC1} 1)
    local Lmv1p_fq1="${prefix32}+_1.fq.gz" # p=plus,+, m=minus,-
    local Lmv1p_fq2="${prefix32}+_2.fq.gz" # p=plus,+, m=minus,-
    local Lmv1m_fq1="${prefix32}-_1.fq.gz" # p=plus,+, m=minus,-
    local Lmv1m_fq2="${prefix32}-_2.fq.gz" # p=plus,+, m=minus,-
    anno_dir ${gene_bed} ${vec1_dir}/v1
    echo 4. ${prefix32}
    echo "[4/x] - find vector-2 in read2"
    local vec2_dir="${out_dir}/tn5/attL/v1"
    # Lpv1pv2(p|m)
    local prefix41=$(trim_3ad ${vec2_dir} ${Lpv1p_fq1} ${Lpv1p_fq2} ${VEC2} 2)
    # Lpv1mv2(p|m)
    local prefix42=$(trim_3ad ${vec2_dir} ${Lpv1m_fq1} ${Lpv1m_fq2} ${VEC2} 2)
    # Lmv1pv2(p|m)
    local prefix43=$(trim_3ad ${vec2_dir} ${Lmv1p_fq1} ${Lmv1p_fq2} ${VEC2} 2)
    # Lmv1mv2(p|m)
    local prefix44=$(trim_3ad ${vec2_dir} ${Lmv1m_fq1} ${Lmv1m_fq2} ${VEC2} 2)
    anno_dir ${gene_bed} ${vec2_dir}/v2
    echo "[5/x] - wrap all statistics"
    local summary="${out_dir}/summary.txt"
    wrap_stat ${out_dir} ${fq1} > ${summary}
    echo "save stat to file: ${summary}"
    echo "[6/x] - finish"
}
export -f hts_y2h


status=$(check_commands)
[[ ${status} = *failed* ]] && echo "!!! error, check above missing files" && exit 1
hts_y2h $@

