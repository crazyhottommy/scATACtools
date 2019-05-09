#! /usr/bin/env Rscript

'Calculate banding score from scATACseq data
Usage:
    calculate_atac_banding_score_per_cell.R [ --barcodeList=<FILE> ] <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --barcodeList=<FILE>  a file with one cell barcode per line without header. e.g. barcodes.tsv from cellranger
Arguments:
    input  input filename of fragment insert sizes for all cells or stdin(-). This is the dataframe
    with a header from output of get_insert_size_distribution_per_cell.py. Three columns of this dataframe.
    cell column contains the cell barcode; insert_size column contains the insert size, read_count
    column contains the number of paired reads for that insert size.
    
    output  output filename or stdout (-)
' -> doc

suppressMessages(library(docopt))
suppressMessages(library(tidyverse))
suppressMessages(library(purrr))

arguments <- docopt(doc, version = 'calculate_atac_banding_score_per_cell.R v1.0\n\n')

OpenRead <- function(arg) {
        if (arg %in% c("-", "/dev/stdin")) {
                file("stdin", open = "r")
        } else if (grepl("^/dev/fd/", arg)) {
                fifo(arg, open = "r")
        } else {
                file(arg, open = "r")
        }
}



dat.con <- OpenRead(arguments$input)
frags <- read.table(dat.con, header = TRUE)

## this is a rewrite of https://github.com/shendurelab/mouse-atac/blob/master/banding_scores/calculate_nucleosome_banding_scores.R
## using tidyverse

get_banding_score<- function(df){
        df<- df %>% 
                filter(insert_size >0) %>%
                #make complete insert size from 1 to 1000
                tidyr::complete(insert_size = 1:1000, fill = list(read_count  = 0)) %>%
                arrange(insert_size)
        periodogram<- spec.pgram(df$read_count / max(df$read_count), pad=0.3, tap=0.5, span=2, plot=F, fast=T)
        periodogram$freq<- 1/periodogram$freq
        
        banding_score<- sum(periodogram$spec[periodogram$freq >= 100 & periodogram$freq <= 300])
        return(banding_score)
}

##  use furrr to speed up?

if (!is.null(arguments$barcodeList)){
        banding_scores<- frags %>% 
                filter(cell %in% arguments$barcodeList)
                group_by(cell) %>%
                nest() %>%
                mutate(banding_score = map_dbl(data, get_banding_score)) %>% 
                select(-data)
        
} else{
        banding_scores<- frags %>% 
                group_by(cell) %>%
                nest() %>%
                mutate(banding_score = map_dbl(data, get_banding_score)) %>% 
                select(-data)
}



OpenWrite<- function(arg){
        if (arg %in% c("-", "/dev/stdout")){
                # empty to specify stdout
                out.con<- stdout()
                writeLines(format_tsv(banding_scores), out.con)
        } else {
                out.con<- file(arg, open = "w+")
                write_tsv(banding_scores, out.con)
                
        }
}

OpenWrite(arguments$output)