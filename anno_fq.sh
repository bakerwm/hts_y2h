#!/usr/bin/bash 

################################################################################
## Global variables
SRC_DIR=$(dirname $(realpath -s $0))
[[ -z ${N_CPU} ]] && N_CPU=8
HG38_IDX="/data/biodata/genome/hg38/hisat2_index/hg38"
# GENE_GTF="/data/biodata/ensembl/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz"
# PARSE_GTF="${SRC_DIR}/parse_gtf.py"
ANNO_PY="${SRC_DIR}/anno_bed.py"
# [[ ! -f ${GENE_GTF} ]] && echo "hg38 gtf not found: ${GENE_GTF}" && exit 1
# [[ ! -f ${PARSE_GTF} ]] && echo "script not found: ${PARSE_GTF}" && exit 1
[[ ! -f ${ANNO_PY} ]] && echo "script not found: ${ANNO_PY}" && exit 1
[[ ! -f ${HG38_IDX}.1.ht2 ]] && echo "hisat2 index not found: ${HG38_IDX}" && exit 1
################################################################################


function gtf2bed() {
    local out_dir=$1
    local gtf=$2 # Ensembl gtf
    local gname=$(basename ${gtf%.gz})
    gname=${gname%.gtf} # trim suffix
    local bed="${out_dir}/${gname}.bed.gz" # bgzipped compressed
    if [[ ! -f ${bed} ]]
    then
        [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
        local bed_raw="${out_dir}/gene_files/${gname}.gene.f-1.pc.g100.bed" # 
        python ${PARSE_GTF} -o ${out_dir}/gene_files -f gene -g protein_coding -gs 100 -fs -1 -i ${gtf}
        # gzipped, tabix
        bgzip -c ${bed_raw} > ${bed}
        tabix -p bed ${bed}
    fi
    [[ ! -f ${bed} ]] && echo "bed file not exists" && return 1
    # echo ${bed}
}
export -f gtf2bed


# Annotate fq, by gene_name/gene_id/NULL ...
# human library
function align_hg38() {
    local out_dir=$1
    local fq=$2 # fasta
    ## Global variable: HG38_IDX
    local fname=$(basename ${fq%.fq.gz})
    local bam="${out_dir}/${fname}.bam"
    local bed="${bam/.bam/.bed}"
    local log="${out_dir}/${fname}.hisat2.log"
    # [[ -f ${bam} ]] && echo "file exists: ${bam}" && return 1
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    [[ ! -f ${bam} ]] && \
        hisat2 -p ${N_CPU} --very-sensitive --add-chrname -x ${HG38_IDX} -U ${fq} 2> ${log} | \
        samtools view -Sub -F 0x4 -F 2048 -q 10 - | \
        samtools sort -@ ${N_CPU} -o ${bam} - && \
        samtools index ${bam}
    # convert bam to bed
    [[ ! -f ${bed} ]] && \
        bedtools bamtobed -i ${bam} > ${bed} #
    # echo ${bed}
}
export -f align_hg38


# para: <out_dir> <gene_bed> <in_bed>
function anno_bed() {
    [[ $# -lt 3 ]] && echo "Usage: anno_bed <out_dir> <gene_bed> <in_bed>" && return 1
    local out_dir=$1
    local gene_bed=$2
    local in_bed=$3 # query
    local bname=$(basename ${in_bed%.bed}.anno.bed) # new name
    local out_bed="${out_dir}/${bname}" # output 
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # prepare bed
    local gname=$(basename ${gene_bed%.gz})
    gname=${gname%.bed}
    # local gname=$(basename ${GENE_GTF%.gz})
    # gname=${gname%.gtf} # trim suffix
    ## move gene_bed to parameter
    # local gene_bed="${out_dir}/db/${gname}.bed.gz" # bgzipped compressed
    # gtf2bed ${out_dir}/db ${GENE_GTF}
    [[ ! -f ${out_bed} ]] && python ${ANNO_PY} ${gene_bed} ${in_bed} ${out_bed}
    [[ -f ${out_bed} ]] && echo ${out_bed} 
}
export -f anno_bed


# annotate fastq, by gene_bed
# para: <out_dir> <gene_bed> <fq>
function anno_fq() {
    [[ $# -lt 3 ]] && echo "Usage: anno_fq <out_dir> <gene_bed> <fq>" && return 1
    local out_dir=$1
    local gene_bed=$2
    local fq_in=$3 # fq.gz
    local fname="$(basename ${fq_in%.fq.gz})"
    local bed_raw="${out_dir}/${fname}.bed"
    [[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
    # 1. map to hg38 genome
    echo "    [1/2] - mapping reads to hg38 genome"
    [[ ! -f ${bed_raw} ]] && align_hg38 ${out_dir} ${fq_in}
    # 2. get annotation (gene_name)
    echo "    [2/2] - annotating mapped reads by gene_name"
    local bed_anno=$(anno_bed ${out_dir} ${gene_bed} ${bed_raw})
    # 3. log
    echo "    ... save to file: ${bed_anno}" #
}
export -f anno_fq


anno_fq $@
