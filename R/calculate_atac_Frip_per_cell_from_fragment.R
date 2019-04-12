#! /usr/bin/env Rscript

'Calculate Frip (Fraction of reads in peaks) from 10xscATACseq fragments.tsv.gz
needs bioconductor libraries: Rsamtools, GenomicRanges, GenomicAlignments, rtracklayer
and R package dplyr, data.table::fread
This is much faster than the calculate_atac_Frip_per_cell_from_bam.R
Takes 10mins for a 1G fragment.tsv.gz file 

Usage:
    calculate_atac_Frip_per_cell_from_fragment.R <fragment> <bed>  <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
Arguments:
    fragment  input filename of 10x scATAC fragment.tsv.gz. zipped.
    bed  input filename for the peaks in bed format
    output output filename for the Frip score per cell
' -> doc

library(docopt)
arguments <- docopt(doc, version = 'calculate_atac_Frip_per_cell_from_fragment.R v1.0\n\n')

suppressMessages(library(Rsamtools))
suppressMessages(library(GenomicRanges))
suppressMessages(library(GenomicAlignments))
suppressMessages(library(rtracklayer))
suppressMessages(library(dplyr))

## this function copied and modified from https://github.com/caleblareau/scATAC_10X_raw_to_kmers/blob/master/example_kmers.R

getCountsFromFrags <- function(frag_gz_file,
                               peaks){
        
        peaks_gr<- import(peaks, format = "BED")
        # Make GRanges of fragments that are solid for the cells that we care about
        frags_valid <- data.table::fread(paste0("zcat < ", frag_gz_file)) %>% 
                data.frame() %>% 
                mutate(V2 = V2 + 1) %>% # make it 1 based for R
                GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
        
        # Get a denominator, per cell
        denom <- table(GenomicRanges::mcols(frags_valid)$V4)
        barcodes_found <- names(denom)
        
        # Get the overlaps with peaks
        ovPEAK <- GenomicRanges::findOverlaps(peaks_gr, frags_valid)
        
        # Establish a numeric index for the barcodes for sparse matrix purposes
        id <- factor(as.character(GenomicRanges::mcols(frags_valid)$V4), levels = barcodes_found)
        
        # Make sparse matrix with counts with peaks by  unique barcode
        countdf <- data.frame(peaks = S4Vectors::queryHits(ovPEAK),
                              sample = as.numeric(id)[S4Vectors::subjectHits(ovPEAK)]) %>%
                dplyr::group_by(peaks,sample) %>% dplyr::summarise(count = n()) %>% data.matrix()
        
        m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks_gr)),
                                  j = c(countdf[,2], length(barcodes_found)),
                                  x = c(countdf[,3],0))
        colnames(m) <- barcodes_found
        
        # Make a polished colData
        colData <- data.frame(
                sample = barcodes_found,
                depth = as.numeric(denom),
                FRIP = Matrix::colSums(m)/as.numeric(denom)
        )
        return(colData)
}

df<- getCountsFromFrags(arguments$fragment, arguments$bed)
write.table(df, arguments$output, sep="\t", row.names = F, col.names = T, quote =F)

