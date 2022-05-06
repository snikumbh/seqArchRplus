

#' @title Plot heatmaps of motifs occuring in seqArchR clusters
#'
#' @param sname The sample name
#'
#' @param seqs The sequences as a DNAStringSet object
#'
#' @param flanks Flank size. The same flank is used upstream and downstream.
#' A vector of values is also accepted when more than oone flanks should be
#' visualized.
#'
#' @param clusts List of sequence Ids in each cluster.
#'
#' @param use_colors Specify colors to use
#'
#' @param motifs A vector of motifs to be visualized in the sequence. This can
#' be any words formed by the
#' \href{https://www.bioinformatics.org/sms2/iupac.html}{IUPAC code}.
#' For example, TATAA, CG, WW, SS etc.
#'
#' @param dir_path The path to the directory
#'
#' @param fheight,fwidth,funits Height and width of the PDF file, and the units
#'  in which they are specified
#'
#' @param n_cores Numeric. If you wish to parallelize annotation of peaks,
#' specify the number of cores. Default is 1 (serial)
#'
#'
#' @return Nothing. PNG images are written to disk using the provided filenames.
#'
#' @importFrom Biostrings width
#'
#' @export
#'
plot_motif_heatmaps <- function(sname, seqs, flanks = c(50), clusts,
                                use_colors = NULL, motifs, dir_path,
                                fheight = 500, fwidth = 500, funits = "px",
                                n_cores = 1){
    cli::cli_h1(paste0("Motif heatmaps"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    parallelize <- FALSE
    if(n_cores > 1) parallelize <- TRUE
    bpparam <- .handle_multicore(crs = n_cores, parallelize = parallelize)
    ##
    nClust <- length(clusts)
    if(is.null(use_colors)){
        message("Using default colors")
        nClust_colors <- .get_ncolors(n = nClust, palname = "Set1")
    }
    ##
    all_motifs_str <- paste0(motifs, collapse="_")

    ###
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##
    ## Iterate over different flanks
    sam <- lapply(flanks, function(x){
        cli::cli_h3(paste0("Using flank: ", x))
        fl_up <- fl_down <- x

        maxL <- Biostrings::width(seqs[1])
        midP <- base::ceiling(maxL/2)
        seqs <- Biostrings::subseq(seqs, start = midP - fl_up,
                                   end = midP + fl_down)
        ##
        seq_order <- unlist(clusts)
        #### Using heatmaps pkg
        patt_hm_list500 <- BiocParallel::bplapply(
                                    motifs,
                                    function(x){
                                      hm <- heatmaps::PatternHeatmap(
                                          seq = seqs[seq_order],
                                          pattern = x,
                                          coords = c(-1*fl_up, fl_down),
                                          label = x)
                                      shm <- heatmaps::smoothHeatmap(hm,
                                                            sigma=c(2,2))
                                    },
                                    BPPARAM=bpparam
        )
        ##
        fname_suffix <- paste0(paste(fl_up, "up", fl_down, "down",
                                     "motifHeatmaps", all_motifs_str,
                                     sep = "_"), ".png")
        fname <- file.path(result_dir_path, fname_suffix)

        clust_lens <- unlist(lapply(clusts, length))
        ##
        grDevices::png(fname, height = fheight, width = fwidth, units = funits)
        pl_hms <- heatmaps::plotHeatmapList(patt_hm_list500,
                                        box.width = 1.3,
                                        cex.label = 1.1,
                                        cex.axis = 0.7, scale.lwd = 0.5,
                                        tcl=-0.3, padj=-1.5,
                                        partition = clust_lens,
                                        partition.legend = TRUE,
                                        partition.lines = TRUE,
                                        partition.col = nClust_colors,
                                        legend = TRUE,
                                        legend.width = 0.3, cex.legend = 0.8,
                                        legend.pos = "r")
        grDevices::dev.off()
        cli::cli_alert_success(paste0("Written to: ", fname))
    })
    ##
}
## =============================================================================
