#! /bin/bash

set -e
set -u
set -o pipefail

## column 9 is the template length according to the SAM specification https://samtools.github.io/hts-specs/SAMv1.pdf
## fragment length <=30 are artifact
## see https://github.com/GreenleafLab/NucleoATAC/issues/18
## the shell script was adopted from Xi Chen https://dbrg77.wordpress.com/2017/02/10/atac-seq-insert-size-plotting/
bam=$1
out=$2

show_help() {
cat << EOF
get the counts of fragment length from a bam file
The first column is the count of the fragment length in the second column

Usage: ${0##*/}  <path to bam file> <output txt file>
Example: ${0##*/} my.bam my.txt
EOF
}

if [[ $# -le 1 ]];then show_help;exit 1;fi
# read input

## check if samtools in path

((command -v samtools) > /dev/null) || \
(echo "[get_fragment_length_from_bam.sh] error: samtools is not in your $PATH, \
	install it by 'conda install -c bioconda samtools' on your machine" && exit 1)

check_file_exists() {
    if [ ! -f "${bam}" ]; then
        echo "[get_fragment_length_from_bam.sh] error: file '$bam' does not exist." >&2
        exit 1
    fi
}

check_file_exists


samtools view $bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
| sed -e 's/^[ \t]*//' | awk '{print $1"\t"$2}' > $out