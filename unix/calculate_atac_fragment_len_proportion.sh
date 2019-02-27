#! /bin/bash
set -euo pipefail

show_help() {
cat << EOF
get the proportion of fragment length fall into a range in ATACseq bam
default is 100bp to 300bp. A successful ATACseq experiment has a 
peak around 200bp in the fragment length distribution.

Usage: ${0##*/}  <path to the fragment length file or stdin> 
Example: 
${0##*/} fragment.txt 
${0##*/} fragment.txt 150 350
samtools view my.bam | awk '{ if (\$9>0) print \$9}' | ${0##*/} - 
EOF
}

if [[ $# == 0 ]]; then show_help; exit 1; fi

file=${1:-/dev/stdin}
low=${2:-100}
high=${3:-300}

awk -v low=$low -v high=$high '$1>=low && $1<=high {++a} END {print a/NR}' $file
