#!/bin/sh
#PBS -l walltime=08:00:00
#PBS -e ${1}/.annotate_sam.sh.e
#PBS -o ${1}/.annotate_sam.sh.o

# Designed to run on all the SAM files in a directory

dr=$(dirname "$0")
[ -z "${1}" ] && echo "Specify destination directory as first argument"
dest=${1}
shift

echo "Handling $*"

for fn in $* ; do
    base=$(echo "$fn" | sed 's/\.sam$//')

    python ${dr}/rep.py \
        --basename ${base} \
        --rep ${dr}/hg19.chr1.fa.out \
        --sam-in ${fn} \
        --exome-bed ${dr}/exome_chr1.bed > ${dest}/${base}.csv

    echo "Correctness-by-overlapped-feature table in ${dest}/${base}.csv"
done
