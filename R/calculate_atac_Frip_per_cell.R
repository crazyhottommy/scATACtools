#! /usr/bin/env Rscript

'Calculate Frip (Fraction of reads in peaks) from 10xscATACseq bam
needs bioconductor libraries: Rsamtools, GenomicRanges, GenomicAlignments, rtracklayer
and R package dplyr.
Usage:
    calculate_atac_Frip_per_cell.R <bam> <bed> [ barcode=<CB> ] <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --barcode=<CB>  tag for cellbarcode [default: CB]
Arguments:
    bam  input filename of 10x merged scATAC bam
    bed  input filename for the peaks in bed format
    output output filename for the Frip score
' -> doc

library(docopt)
arguments <- docopt(doc, version = 'calculate_atac_Frip_per_cell.R v1.0\n\n')

## the function was copied and modified from https://github.com/caleblareau/chromVARxx/blob/master/R/getCounts-tweaks.R
## thanks Caleb for sharing! bedtools maybe able to do the same thing, but need to test.

suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))

get_counts_by_cellbarcode<- function(bamfile, peaks, barcodeTag, mapqFilter = 0){
        
        # Import alignments and get overlaps with peaks
        GA <- readGAlignments(bamfile, use.names = TRUE, param = ScanBamParam(
                tag = c(barcodeTag), mapqFilter = 0))
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
        
        # Generate colData
        depth <- data.frame(
                sample = as.numeric(id),
                read = names(GA)) %>%
                distinct() %>% group_by(sample) %>% summarise(depth = n()) %>% data.matrix()
        colData <- data.frame(
                sample = uniqueBarcodes,
                depth = depth[,2],
                FRIP = Matrix::colSums(m)/depth[,2])
        return(colData)
}


df<- get_counts_by_cellbarcode(arguments$bam, arguments$bed, barcodeTag = arguments$barcode)
write.table(df, arguments$output, sep="\t", row.names = F, col.names = T, quote =F)




        