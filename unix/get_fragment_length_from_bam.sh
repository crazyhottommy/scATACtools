#! /bin/bash

set -e
set -u
set -o pipefail

## column 9 is the template length according to the SAM specification https://samtools.github.io/hts-specs/SAMv1.pdf
## fragment length <=30 are artifact
## see https://github.com/GreenleafLab/NucleoATAC/issues/18
## the shell script was adopted from Xi Chen https://dbrg77.wordpress.com/2017/02/10/atac-seq-insert-size-plotting/
show_help() {
cat << EOF
get the fragment length of every read pair from a bam file
The output can be plotted to identify the ~200 bp peak 
see https://github.com/GreenleafLab/NucleoATAC/issues/18

Usage: ${0##*/}  <path to bam file> <output txt file>
Example: ${0##*/} my.bam my.txt
EOF
}

if [[ $# -le 1 ]];then show_help;exit 1;fi

## check if samtools in path

((command -v samtools) > /dev/null) || \
(echo "[get_fragment_length_from_bam.sh] error: samtools is not in your $PATH, \
	install it by 'conda install -c bioconda samtools' on your machine" && exit 1)

bam=$1
out=$2

check_file_exists() {
    if [ ! -f "${bam}" ]; then
        echo "[get_fragment_length_from_bam.sh] error: file '$bam' does not exist." >&2
        exit 1
    fi
}

check_file_exists

#samtools view $bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
#| sed -e 's/^[ \t]*//' | awk '{print $1"\t"$2}' > $out

# this is better for downstream R plotting using the raw data
samtools view $bam | awk '$9>0' | cut -f 9 > $out


