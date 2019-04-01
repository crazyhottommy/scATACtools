#! /usr/bin/env Rscript
'plot smooth scatter plot for ATACseq data, x-axis is library size, y-axis is
FRIP score. The input of this script can be obtained by using calculate_atac_Frip_per_cell_from_fragment.R
Usage:
    plot_FRIP_scatter.R <FRIP> <barcode> <output>
    
Options:
    -h --help  Show this screen.
    -v --version  Show version.
Arguments:
    FRIP     input filename of FRIP score for each cell. A dataframe of 3 columns with header:
             sample, depth, FRIP. column 1 is the cell barcode, column 2 is the library size, 
             column 3 is the FRIP score.
    barcode  input filename of valid barcode (e.g. from cellranger filtered_peak_bc_matrix/barcodes.tsv). 
             A one column dataframe without header.
    output  output pdf filename
' -> doc

library(docopt)
arguments <- docopt(doc, version = 'plot_FRIP_scatter.R v1.0\n\n')


suppressMessages(library(tidyverse))
# devtools::install_github('VPetukhov/ggrastr')
#library(ggrastr)
library(ggExtra)
frip<- read_tsv(arguments$FRIP)

passed_barcodes<- read_tsv(arguments$barcode, col_names =F)

frip<- frip %>% mutate(cellranger_pass = if_else(sample %in% passed_barcodes$X1, "yes", "no"))

## with too many points, raster them
#ggplot(frip %>% filter(FRIP <=1), aes(x = log10(depth), y = FRIP)) +
#        geom_point_rast(aes(col = cellranger_pass), alpha = 0.2) + 
#        scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
#        theme_bw(base_size = 14)
#ggsave("results/ATAC_FRIP_scattersmooth_nofilter.pdf", width = 8, height = 6)


filtered_frip<- frip %>% filter(FRIP <=1, cellranger_pass == "yes")

d<- densCols(x = log10(filtered_frip$depth), y = filtered_frip$FRIP, colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
p<- ggplot(filtered_frip, aes(x = log10(depth), y = FRIP)) +
        geom_point(alpha = 0.2, col = d) + 
        scale_color_identity() +
        geom_hline(yintercept = 0.6, linetype = 2, col = "red", size = 1) +
        geom_vline(xintercept = 3, linetype = 2, col = "red", size = 1) +
        theme_bw(base_size = 14)

# marginal density
p<- ggMarginal(p, type="density", fill = "blue")
ggsave(arguments$output, plot = p, width = 6, height = 6)
