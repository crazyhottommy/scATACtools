#! /bin/env python3

import pysam
import csv
import argparse
import os.path
import sys

parser = argparse.ArgumentParser()
parser.add_argument("bam", help="Required. the FULL path to the 10x scATAC bam file generated \
    by cellranger-atac count")
parser.add_argument("-prefix", help="Optional, the prefix of the output bam, default is barcode.bam")
args = parser.parse_args()


if os.path.exists(args.bam):
    pass
else:
    print("10x scATAC bam not found")
    sys.exit(1)


fin = pysam.AlignmentFile(args.bam, "rb")

## dict with keys are barcode, values are outbam
fouts_dict = {}

for read in fin:
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
        if args.prefix:
            fout_name = args.prefix + "_" + cell_barcode + ".bam"
        else:
            fout_name = cell_barcode + ".bam"    
        if cell_barcode not in fouts_dict:
            fouts_dict[cell_barcode] = pysam.AlignmentFile(fout_name, "wb", template = fin)
        fouts_dict[cell_barcode].write(read)
    else: 
        continue
    

## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()
