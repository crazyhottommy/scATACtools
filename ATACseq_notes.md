

### TSS enrichment score

from the supplementary method section of this paper [The chromatin accessibility landscape of primary human cancers](https://science.sciencemag.org/content/362/6413/eaav1898) 

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

from the method of [Massively parallel single-cell chromatin landscapes of human immune cell development and intratumoral T cell exhaustion](https://www.biorxiv.org/content/10.1101/610550v1)

>Filtering Cells by TSS Enrichment and Unique Fragments
Enrichment of ATAC-seq accessibility at TSS was used to quantify data quality without
the need for a defined peak set. Calculating enrichment at TSS was performed as
previously described48, and TSS positions were acquired from the Bioconductor package
from “TxDb.Hsapiens.UCSC.hg19.knownGene”. Briefly, Tn5 corrected insertions were
aggregated +/- 2,000 bp relative (TSS strand-corrected) to each unique TSS genome
wide. Then this profile was normalized to the mean accessibility +/- 1,900-2,000 bp from
the TSS and smoothed every 51bp in R. The calculated TSS enrichment represents the
max of the smoothed profile at the TSS. We then filtered all single cells that had at least
1,000 unique fragments and a TSS enrichment of 8 for all data sets.

Note that each study uses a different cutoff of the TSS enrichment score.

* [ENCODE definition](https://www.encodeproject.org/data-standards/terms/#enrichment) not so clear to me..

>Transcription Start Site (TSS) Enrichment Score - The TSS enrichment calculation is a signal to noise calculation. The reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 1000 bp in either direction (for a total of 2000bp). This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth. This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome) there should be an increase in signal up to a peak in the middle. We take the signal value at the center of the distribution after this normalization as our TSS enrichment metric. Used to evaluate ATAC-seq. 

