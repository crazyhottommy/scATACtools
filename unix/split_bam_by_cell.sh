#! /bin/bash
set -e
set -u
set -o pipefail

## first, use split_bam_by_cluster.py to split cells by clusters
## if one has 30 clusters, this will result in 30 bam files

## then for each cluster level bam, one can further split to cell level bam file
## if one has 500 cells in pbmc_cluster_4.bam, split_bam_by_cell.py will result in 500 bam files.

## when one has thousands of files in a single directory, even ls takes long time
## so better to split the cell level bam files to a seprate directory for each cluster.

bam=$1
outputdir=$(basename $bam .bam)
./split_bam_by_cell.py $bam -prefix $outputdir -outdir $outputdir