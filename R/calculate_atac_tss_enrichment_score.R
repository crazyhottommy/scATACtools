#' checkClass function
#' 
#' check whether the x object corresponds to the given class
#'
#' @param x object
#' @param class.name class name
#' @param var.name uses x object
#' @keywords internal
checkClass = function(x, class.name, var.name = deparse(substitute(x))){
  
  fun.name = match.call(call=sys.call(sys.parent(n=1)))[[1]]
  if(!class(x) %in% class.name)
    stop(paste(fun.name,': ', 
               var.name, 
               ' is not of class: ', 
               paste(class.name, collapse=' '), 
               '\n', sep=''))
}

### remove the tss that do not have coverage
### I took some code from the ScoreMatrix.R function in the genomation package.
### give the credit due :)
### see https://github.com/BIMSBbioinfo/genomation/blob/master/R/scoreMatrix.R#L113
constrainRanges = function(target, windows){
  
  checkClass(target, c('SimpleRleList','RleList','CompressedRleList'))
  checkClass(windows, 'GRanges')
  
  mcols(windows)$X_rank = 1:length(windows)
  r.chr.len = elementNROWS(target)
  constraint = GRanges(seqnames=names(r.chr.len),
                       IRanges(start=rep(1,length(r.chr.len)),
                               end=as.numeric(r.chr.len)))
  # suppressWarnings is done becuause GenomicRanges function give warnings 
  #if you don't have the same seqnames in both objects
  win.list.chr = suppressWarnings(subsetByOverlaps(windows, 
                                                   constraint,
                                                   type = "within",
                                                   ignore.strand = TRUE))
  
  if(length(win.list.chr) == 0)
    stop('All windows fell have coordinates outside windows boundaries')
  return(win.list.chr)
}


#' Calculate tss enrichment score from 10xscATAC fragment file
#'
#' @param frag_gz_file  fragment.tsv.gz file from 10x cellranger-atac output or 
#' anyother tool but in the same format.
#' @param cut_site whether or not use the the end of the read as a cutting site.
#' @param txs  a txdb object
#' @param flank flanking bp of tss (upstream and downstream)
#' @param endFlank  bp end flanks of flank for local noise control
#'     flank               flank
#'  ---------------|-----------------
#'                tss
#'  ---                           ---
#'  endFlank                     endFlank
#'  
#' @param highest_tss_flank bp flanking tss windown for choosing the highest tss score.
#' The highest tss enrichment score is not always exactly at tss.
#' @param barcodeList valid barcode list, a file with one column 
#' @param smooth window size to smooth 
#' @param strand.aware consider tss strandness when calculating 
#'
#' @return
#' @export
#'
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' scores<- TssEnrichmentFromFrags("fragment.tsv.gz", txs = txs)

TssEnrichmentFromFrags <- function(frag_gz_file,
                               cut_site = TRUE,
                               txs,
                               flank = 1000,
                               endFlank = 100,
                               highest_tss_flank= 50,
                               smooth = 50,
                               strand.aware = TRUE,
                               workers = 1,
                               barcodeList = NULL){
        
        # Make GRanges of fragments that are solid for the cells that we care about
        frags_valid <- data.table::fread(paste0("zcat < ", frag_gz_file)) %>% 
                data.frame() %>% 
                mutate(V2 = V2 + 1) %>% # make it 1 based for R
                GenomicRanges::makeGRangesFromDataFrame(seqnames.field = "V1", start.field = "V2", end.field = "V3", keep.extra.columns = TRUE)
        if (!is.null(barcodeList)){
                validBarcodes<- read_tsv(barcodeList, col_names = F)
                frags_valid<- frags_valid[frags_valid$V4 %in% validBarcodes$X1]
        }
        
        if(cut_site){
          frags_valid <- c(
            GRanges(seqnames = seqnames(frags_valid), ranges = IRanges(start(frags_valid), start(frags_valid)), RG = mcols(fragments)[,"V4"]),
            GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,"V4"])
          ) }
    
        
        # common chromosome names 
        seqlev<- intersect(seqlevels(frags_valid), seqlevels(txs))
        frags_valid<- keepSeqlevels(frags_valid, seqlev, pruning.mode="coarse")
        
        # calculate coverage per cell
        frags_valid_per_cell<- split(frags_valid, frags_valid$V4)
        # can add the chromose length as additional argument for coverage
        # to get 0 coverages if there are no reads there.
        # ...not for now.
        # this step can taek minutes 
        cvgs<- lapply(frags_valid_per_cell, function(x) coverage(x))
        
        ## return only transcripts with a coverage 
        txs<- keepSeqlevels(txs, seqlev, pruning.mode="coarse")
        txs <- unique(txs)
        
        txs.flanks<- promoters(txs, upstream = flank, 
                            downstream = flank)
        txs.length<- length(txs.flanks)
        
        #lapply(cvgs, TssEnrichmentSingleCell(cvg, txs.flanks))
        
        multicoreParam <- BiocParallel::MulticoreParam(workers = workers)
        BiocParallel::register(multicoreParam)
        TssEnrichmentScores<- BiocParallel::bplapply(cvgs, TssEnrichmentSingleCell, txs.flanks)
}    

TssEnrichmentSingleCell<- function(cvg, txs.flanks){
        ## remove tss not in the coverage and assign a unique id for each tss: X_rank
        txs.flanks<- constrainRanges(cvg, txs.flanks)
        
        if(length(txs.flanks)!=txs.length){
              warning(paste0(txs.length-length(txs.flanks),
                             " Tss removed because they fall out of the coverage"))
            }
        
        # convert GRanges to IntergerRangesList does not maintain the order
        # a unique id was given for each Ranges
        myViews<- Views(cvg, as(txs.flanks, "IntegerRangesList"))
        mat = lapply(myViews,function(x) t(viewApply(x,as.vector)) )
        mat = do.call("rbind",mat)
        
        r.list=split(mcols(txs.flanks)[,"X_rank"], as.vector(seqnames(txs.flanks))  )
        r.list=r.list[order(names(r.list))]
        ranks=do.call("c",r.list)
        rownames(mat) = ranks
        
        if(strand.aware == TRUE){
              orig.rows=txs.flanks[strand(txs.flanks) == '-',]$X_rank
              mat[rownames(mat) %in% orig.rows,] = mat[rownames(mat) %in% 
                                                         orig.rows, ncol(mat):1]
        }
        
        # reorder according to the original Granges (txs)
        mat = mat[order(ranks),]
        
  
        ### normlization by the endFlank local noise
        profile <- colSums(mat)
        # Handles low depth cells
        # for each tss
        #profile_mat_norm <- apply(mat, 1, function(x) x/max(mean(x[c(1:endFlank,(flank*2-endFlank):(flank*2))]), 0.5)) 
        profile_norm <- profile/mean(profile[c(1:endFlank,(flank*2-endFlank+1):(flank*2))])

        #smooth
        #profile_mat_norm_smooth <- apply(profile_mat_norm, 1, function(x) zoo::rollmean(x, smooth, fill = 1))
        profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)
        

        #enrichment
        max_finite <- function(x){
        suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
        }
        
        #e_mat <- apply(profile_mat_norm_smooth, 1, function(x)
        #       max_finite(x[(flank-highest_tss_flank):(flank+highest_tss_flank)]))
        #names(e_mat) <- rownames(profile_mat_norm_smooth)
        e <- max_finite(profile_norm_smooth[(flank-highest_tss_flank):(flank+highest_tss_flank)])
        return(e)
}

