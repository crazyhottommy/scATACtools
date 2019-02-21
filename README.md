# scATACutils
utilities for 10x scATAC data 

* `split_scATAC_bam_by_cluster.py` split a 10x scATAC bam file by cluster id.
see a detailed blog post at https://divingintogeneticsandgenomics.rbind.io/post/split-a-10xscatac-bam-file-by-cluster/

for a bam file size of 20G containing 5000 cells and 12 clusters, it takes ~3 hours.

* `split_scATAC_bam_by_cell.py` split a 10x scATAC bam file by cell barcode.

for a bam file size of 1G (one cluster) containing 285 cells, it takes ~7 mins.