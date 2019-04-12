#! /usr/bin/env Rscript

'Calculate Frip (Fraction of reads in peaks) from 10xscATACseq bam
needs bioconductor libraries: Rsamtools, GenomicRanges, GenomicAlignments, rtracklayer
and R package dplyr.
It takes ~10 hours for a 20G 10x 5k pbmc bam file.
use calculate_atac_Frip_per_cell_from_fragment.R if you want a faster solution.

Usage:
    calculate_atac_Frip_per_cell.R [ --barcode=<CB> --barcodeList=<FILE>]  <bam> <bed> <output>
    
Options:
    -h --help       Show this screen.
    -v --version    Show version.
    --barcode=<CB>  tag for cellbarcode [default: CB]
    --barcodeList=<FILE>   input filename for the valid barcodes, one column dataframe without header
                    (e.g. from cellranger filtered_peak_bc_matrix/barcodes.tsv), if not provided,
                    all cell barcodes in the bam file will be calculated.
Arguments:
    bam  input filename of 10x merged scATAC bam
    bed  input filename for the peaks in bed format
    output output filename for the Frip score
' -> doc

library(docopt)
arguments <- docopt(doc, version = 'calculate_atac_Frip_per_cell_from_bam.R v1.0\n\n')

## the function was copied and modified from https://github.com/caleblareau/chromVARxx/blob/master/R/getCounts-tweaks.R
## thanks Caleb for sharing! bedtools maybe able to do the same thing, but need to test.

suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))
suppressMessages(library(readr))

get_counts_by_cellbarcode<- function(bamfile, peaks, barcodeTag, barcodeList, validBarcodes, mapqFilter = 0){
        
        # Import alignments and get overlaps with peaks
        GA <- readGAlignments(bamfile, use.names = TRUE, param = ScanBamParam(
                tag = c(barcodeTag), mapqFilter = 0))
        ## some reads in the bam do not have CB tag, filter out
        GA<- GA[!is.na(mcols(GA)[,barcodeTag])]
        
        ## filter only the valid barcodes if barcodeList is given
        if (!is.null(barcodeList)){
            validBarcodes<- read_tsv(barcodeList, col_names = F)
            GA<- GA[mcols(GA)[, barcodeTag] %in% validBarcodes$X1]
        }
        
        peaks_df<- read_tsv(peaks, col_names = F)
        colnames(peaks_df)<- c("chr", "start", "end")
        mat_rownames<- paste(peaks_df$chr, peaks_df$start, peaks_df$end, sep = ":")

        peaks<- import(peaks, format= "BED")
        ovPEAK <- findOverlaps(peaks, GA)
        
        # Determine universe of unique barcodes
        barcodes <- as.character(mcols(GA)[,barcodeTag])
        uniqueBarcodes <- unique(barcodes)
        id <- factor(barcodes, levels = uniqueBarcodes)
        
        # Assemble counts
        countdf <- data.frame(peaks = queryHits(ovPEAK),
                              sample = as.numeric(id)[subjectHits(ovPEAK)],
                              read = names(GA)[subjectHits(ovPEAK)]) %>%
                distinct() %>%  # by filtering on distinct read / peak / sample trios, ensure that PE reads that overlap peak aren't double counted
                select(-one_of("read")) %>% 
                group_by(peaks,sample) %>% summarise(count = n()) %>% data.matrix()
        
        m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks)),
                                  j = c(countdf[,2], length(uniqueBarcodes)),
                                  x = c(countdf[,3],0))
        colnames(m) <- uniqueBarcodes
        
    
        m<- as.matrix(m)
        rownames(m)<- mat_rownames
        # Generate colData
        depth <- data.frame(
                sample = as.numeric(id),
                read = names(GA)) %>%
                distinct() %>% group_by(sample) %>% summarise(depth = n()) %>% data.matrix()
        colData <- data.frame(
                sample = uniqueBarcodes,
                depth = depth[,2],
                FRIP = Matrix::colSums(m)/depth[,2])
        return(list(m = m, colData= colData))
}


df_list<- get_counts_by_cellbarcode(arguments$bam, arguments$bed, 
                                    barcodeTag = arguments$barcode,
                                    barcodeList = arguments$barcodeList)
write.table(df_list$colData, arguments$output, sep="\t", row.names = F, col.names = T, quote =F)
write.table(df_list$m, file = "matrix.txt", col.names = T, row.names = T, sep = "\t", quote =F)