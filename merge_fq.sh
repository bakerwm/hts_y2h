# Merge two fastq files by read id
# output only common reads (exists in both files)
# Example: 
# bash merge_fq.sh output a.fq b.fq 
## output fasta
[[ $# -lt 3 ]] && echo "Usage: bash merge_fq.sh <out_dir> <fq1> <fq2> " && exit 1
out_dir=$1
fq1=$2
fq2=$3
[[ ! -d ${out_dir} ]] && mkdir -p ${out_dir}
out1="${out_dir}/$(basename ${fq1})"
out2="${out_dir}/$(basename ${fq2})"
out1=${out1/.fq/.fa}
out2=${out2/.fq/.fa}
[[ ! -f ${out1} ]] && \
    bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print ">"$name"\n"a[$name]"\n"}}' ${fq1} ${fq2} | gzip > ${out1}

[[ ! -f ${out2} ]] && \
    bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print ">"$name"\n"$seq"\n"}}' ${fq1} ${fq2} | gzip > ${out2}
# # table
# tab="${out_dir}/$(basename ${fq1/.fq}.table.txt)"
# [[ ! -f ${tab} ]] && \
#     bioawk -cfastx 'FNR==NR {a[$name]=$seq; next} {if($name in a){print $name,a[$name],$seq}}' ${fq1} ${fq2} | gzip > ${tab}
    

#EOF
