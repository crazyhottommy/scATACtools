#! /usr/bin/env Rscript
'Calculate proportion of fragment length fall in a range from ATACseq data
Usage:
    calculate_atac_fragment_len_proportion.R [--low=<bp> --high=<bp>] <input> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
    --low=<bp>  low cutoff [default: 100]
    --high=<bp>  high cutoff [default: 300]
Arguments:
    input  input filename of fragment length in a one column dataframe without header or stdin (-)
    output  output filename or stdout (-)
' -> doc

library(docopt)
arguments <- docopt(doc, version = 'calculate_atac_fragment_len_proportion.R v1.0\n\n')

OpenRead <- function(arg) {
        if (arg %in% c("-", "/dev/stdin")) {
                file("stdin", open = "r")
        } else if (grepl("^/dev/fd/", arg)) {
                fifo(arg, open = "r")
        } else {
                file(arg, open = "r")
        }
}

OpenWrite<- function(arg){
        if (arg %in% c("-", "/dev/stdout")){
                # empty to specify stdout
                stdout()
        } else {
                file(arg, open = "w+")
        }
}

dat.con <- OpenRead(arguments$input)
out.con<- OpenWrite(arguments$output)
fragment <- read.table(dat.con, header = FALSE)

len<- fragment$V1

ratio = sum(len >=arguments$low & len <= arguments$high)/length(len)
writeLines(as.character(ratio), out.con)


