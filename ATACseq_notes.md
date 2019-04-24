

### TSS enrichment score

from the supplementary method section of this paper: 

> the Library quality was assessed primarily by the transcription start site (TSS) enrichment
score (see below). A cutoff of “4” was used to determine whether a library was of sufficient quality
to deep sequence. We note, however, that the absolute number of this score depends on the set of
transcription start sites used and this metric may not necessarily be comparable across different
studies. Here, a transcription start site enrichment score of 4 corresponds roughly to a fraction of 
reads in peaks of 15%. However, calculation of the fraction of reads in peaks requires predetermination 
of the peak set to be used and this makes the fraction of reads in peaks a less desirable method of 
quality control when studying novel samples. In this paper, we used the transcripts from “TxDb.Hsapiens.UCSC.hg38.knownGene” to define our transcription start sites and found that the enrichment 
scores from these transcripts were reproducible across genome builds hg19 and hg38.

* ATAC-seq data QC – Transcription start site enrichment, fragment length distribution, and fraction
of reads in peaks

> Enrichment of ATAC-seq accessibility at transcription start sites (TSSs) was used to robustly
quantify ATAC-seq data quality without the need for a defined peak set (which is not available
until all samples have been fully sequenced). First, BAM files were read into a Genomic Ranges
object in R using Rsamtools “scanbam” and then corrected by a constant offset to the read start
(“+” stranded +4 bp, “-” stranded -5 bp). To get the fragment length distribution, the width of each
fragment/GRange was plotted. To get the TSS enrichment profile, each TSS from the R package
“TxDb.Hsapiens.UCSC.hg38.knownGene” (accessed by transcripts(TxDb)) was extended 2000
bp in each direction and overlapped with the insertions (each end of a fragment) using
“findOverlaps”. Next, the distance between the insertions and the strand-corrected TSS was
calculated and the number of insertions occurring in each single-base bin was summed. To
normalize this value to the local background, the accessibility at each position +/- 2000 bp from
the TSS was normalized to the mean of the accessibility at positions +/-1900-2000 bp from the
TSS. The final TSS enrichment reported was the maximum enrichment value within +/- 50 bp of
the TSS after smoothing with a rolling mean every 51 bp. To calculate the fraction of reads in
peaks, we calculated the number of insertions that overlapped a peak using “countOverlaps” and
divided by the total number of insertions (number of fragments x 2, as each paired-end read
represents an individual Tn5 insertion).
