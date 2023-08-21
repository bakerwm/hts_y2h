#!/usr/bin/bash 

# annotation for paired fastq files

################################################################################
## Global variables
SRC_DIR=$(dirname $(realpath -s $0))
# [[ -z ${N_CPU} ]] && N_CPU=8
ANNO_FQ="${SRC_DIR}/anno_fq.sh"
[[ ! -f ${ANNO_FQ} ]] && echo "script not found: ${ANNO_FQ}" && exit 1
################################################################################


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
        awk 'FNR==NR {a[$4]=$0; next} {if($4 in a){print $0,a[$4]}}' ${bed1} ${bed2} > ${out1}
    # # file-2:
    # out1="${out_dir}/$(basename ${bed1/.bed/.read12.bed})"
    # [[ ! -f ${out2} ]] && \
    #     awk 'FNR==NR {a[$4]=$8; next} {if($4 in a){print $4,$8,a[$4]}}' ${bed2} ${bed1} > ${out2}
    ## log
    echo "    ... save $(wc -l ${out1}) records"
}
export -f merge_bed_by_id


function anno_fq_pair() {
    [[ $# -lt 4 ]] && echo "Usage: anno_fq_pair <out_dir> <gene_bed> <fq1> <fq2>" && return 1
    local out_dir=$1
    local gene_bed=$2
    local fq1=$3
    local fq2=$4
    # 1. anno each fq file
    local fq1_name="$(basename ${fq1%.fq.gz})"
    local bed1="${out_dir}/${fq1_name}.anno.bed"
    bash ${ANNO_FQ} ${out_dir} ${gene_bed} ${fq1}
    ## fq2
    local fq2_name="$(basename ${fq2%.fq.gz})"
    local bed2="${out_dir}/${fq2_name}.anno.bed"
    bash ${ANNO_FQ} ${out_dir} ${gene_bed} ${fq2}
    # 2. paring read1,read2
    local bed_pe="${bed1%.bed}.read12.bed"
    merge_bed_by_id ${out_dir} ${bed1} ${bed2}
}
export -f anno_fq_pair


anno_fq_pair $@
