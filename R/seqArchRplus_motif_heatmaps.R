

#' @title Plot heatmaps of motifs occurring in seqArchR clusters
#'
#' @param sname The sample name
#'
#' @param seqs The sequences as a DNAStringSet object
#'
#' @param flanks Flank size. The same flank is used upstream and downstream.
#' A vector of values is also accepted when more than one flanks should be
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
#' @examples
#'
#' library(Biostrings)
#' raw_seqs <- Biostrings::readDNAStringSet(
#'                           filepath = system.file("extdata",
#'                             "example_promoters200.fa.gz",
#'                             package = "seqArchRplus",
#'                             mustWork = TRUE)
#'                         )
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' plot_motif_heatmaps(sname = "sample1", seqs = raw_seqs,
#'                     flanks = c(10, 20, 50),
#'                     clusts = use_clusts,
#'                     motifs = c("WW", "SS", "TATAA", "CG", "Y"),
#'                     dir_path = tempdir(),
#'                     fheight = 800, fwidth = 1600)
#'
#' @author Sarvesh Nikumbh
plot_motif_heatmaps <- function(sname, seqs, flanks = c(50), clusts,
                                use_colors = NULL, motifs, dir_path,
                                fheight = 500, fwidth = 500, funits = "px",
                                n_cores = 1) {
    cli::cli_h1(paste0("Motif heatmaps"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    parallelize <- FALSE
    if (n_cores > 1) parallelize <- TRUE
    bpparam <- .handle_multicore(crs = n_cores, parallelize = parallelize)
    ##
    nClust <- length(clusts)
    if (is.null(use_colors)) {
        message("Using default colors")
        nClust_colors <- .get_ncolors(n = nClust, palname = "Set1")
    }
    ##
    all_motifs_str <- paste0(motifs, collapse = "_")

    ##
    if(!is.null(dir_path))
        result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##
    ## Iterate over different flanks
    sam <- lapply(flanks, function(x) {
        cli::cli_h3(paste0("Using flank: ", x))
        fl_up <- fl_down <- x

        maxL <- Biostrings::width(seqs[1])
        midP <- base::ceiling(maxL / 2)
        seqs <- Biostrings::subseq(seqs,
            start = midP - fl_up,
            end = midP + fl_down
        )

        ##
        seq_order <- unlist(clusts)
        #### Using heatmaps pkg
        patt_hm_list500 <- BiocParallel::bplapply(
            motifs,
            function(x) {
                hm <- heatmaps::PatternHeatmap(
                    seq = seqs[seq_order],
                    pattern = x,
                    coords = c(-1 * fl_up, fl_down),
                    label = x
                )
                shm <- heatmaps::smoothHeatmap(hm,
                    sigma = c(2, 2)
                )
            },
            BPPARAM = bpparam
        )
        ##
        if(!is.null(dir_path)){
            fname_suffix <- paste0(paste(fl_up, "up", fl_down, "down",
                "motifHeatmaps", all_motifs_str,
                sep = "_"
            ), ".png")
            fname <- file.path(result_dir_path, fname_suffix)
        }

        clust_lens <- lengths(clusts)

        ##


        if(!is.null(dir_path)){
            grDevices::png(fname, height = fheight, width = fwidth,
                            units = funits)
        }
        pl_hms <- heatmaps::plotHeatmapList(
            patt_hm_list500,
            box.width = 1.3,
            cex.label = 1.1,
            cex.axis = 0.7, scale.lwd = 0.5,
            tcl = -0.3, padj = -1.5,
            partition = clust_lens,
            partition.legend = TRUE,
            partition.lines = TRUE,
            partition.col = nClust_colors,
            legend = TRUE,
            legend.width = 0.3,
            cex.legend = 0.8,
            legend.pos = "r"
        )
        if(!is.null(dir_path)){
            grDevices::dev.off()
            cli::cli_alert_success(paste0("Written to: ", fname))
        }


    })
    ##
}
## =============================================================================


.get_cluster_legend_plot <- function(clusts, use_colors, ext = ".png",
                                        result_dir_path,
                                        fwidth = 200, fheight = 800){
    # Cluster partitions
    nClust <- length(clusts)
    clust_lens <- rev(lengths(clusts))

    if (is.null(use_colors)) {
        message("Using default colors")
        nClust_colors <- .get_ncolors(n = nClust, palname = "Set1")
    }
    ##

    opts <- heatmaps::heatmapOptions()
    opts$partition <- clust_lens
    opts$colors <- nClust_colors
    clust_pl_fname <- file.path(result_dir_path,
                                paste0("motif_heatmaps_ClusteringLegend", ext))
    png(filename = clust_pl_fname, width = fwidth, height = fheight,
        units = "px")
    heatmaps::plot_clusters(opts)
    dev.off()

    clust_pl_fname
}
## =============================================================================


## Using seqPattern and cowplot
## TODO: different flanks

#' @title Plot heatmaps of motifs occurring in seqArchR clusters
#'
#' @description This function uses the seqPattern package. It is
#'  recommended to use this function rather than `plot_motif_heatmaps`
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
#' @param fheight,fwidth Height and width of the individual heatmap
#' plots in pixels. Default values are 500 and 300
#'
#' @param n_cores Numeric. If you wish to parallelize annotation of peaks,
#' specify the number of cores. Default is 1 (serial)
#'
#' @param type Specify either of "png" or "jpg" to obtain PNG or JPEG files as
#' output
#'
#' @param ... Additional arguments passed to
#' \code{\link[seqPattern]{plotPatternDensityMap}}
#'
#' @return Nothing. Images are written to disk using the provided filenames.
#' In addition, two legends are printed to separate files: the color legend
#' and the clustering legend, which can be then combined with the heatmaps.
#' The heatmaps themselves have the cluster numbers marked on the vertical axis.
#'
#'
#' @importFrom Biostrings width
#' @importFrom seqPattern plotPatternDensityMap
#' @importFrom grDevices png
#' @importFrom ggimage as.ggplot
#' @importFrom magick image_read
#'
#' @export
#'
#' @examples
#'
#' library(Biostrings)
#' raw_seqs <- Biostrings::readDNAStringSet(
#'                           filepath = system.file("extdata",
#'                             "example_promoters200.fa.gz",
#'                             package = "seqArchRplus",
#'                             mustWork = TRUE)
#'                         )
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' plot_motif_heatmaps2(sname = "sample1", seqs = raw_seqs,
#'                     flanks = c(50),
#'                     clusts = use_clusts,
#'                     motifs = c("WW", "SS"),
#'                     dir_path = tempdir())
#'
#' @author Sarvesh Nikumbh

plot_motif_heatmaps2 <- function(sname,
                            seqs, flanks = c(50), clusts,
                            use_colors = NULL, motifs, dir_path,
                            fheight = 1.5*fwidth, fwidth = 2000,
                            hm_scale_factor = 0.75,
                            n_cores = 1, type = c("png", "jpg") , ...){

    cli::cli_h1(paste0("Motif heatmaps"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    use_ext <- paste0(".", match.arg(type))

    parallelize <- FALSE
    if (n_cores > 1) parallelize <- TRUE
    # bpparam <- .handle_multicore(crs = n_cores, parallelize = parallelize)
    ##
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##
    all_motifs_str <- paste0(motifs, collapse = "_")

    ##
    clust_pl_fname <- .get_cluster_legend_plot(clusts, use_colors,
                            ext = use_ext, result_dir_path,
                            fwidth = 200, fheight = 0.75*fheight)


    ######
    ## Iterate over different flanks
    sam <- lapply(flanks, function(x) {
        cli::cli_h3(paste0("Using flank: ", x))
        fl_up <- fl_down <- x

        maxL <- Biostrings::width(seqs[1])
        midP <- base::ceiling(maxL / 2)
        seqs <- Biostrings::subseq(seqs,
            start = midP - fl_up,
            end = midP + fl_down
        )

        ##
        seq_order <- unlist(clusts)
        #### Using seqPattern pkg
        fname_suffix <- paste(fl_up, "up", fl_down, "down",
            "motifHeatmaps", all_motifs_str,
            sep = "_"
        )

        fname <- file.path(result_dir_path, fname_suffix)
        fname_w_ext <- file.path(result_dir_path, paste0(fname_suffix, use_ext))

        clust_lens <- lengths(clusts)

        use_ticks <- cumsum(clust_lens)

        ##
        use_tempdir <- tempdir()
        use_outFile <- file.path(
            use_tempdir,
            paste0("TStamp_", format(Sys.time(), "%Y-%m-%d_%H_%M_%S"),
                "_PatternDensityMap"))

        use_xticksAt <- c(1, Biostrings::width(seqs[1])/2,
            Biostrings::width(seqs[1])-1)
        use_yticksAt <- rev(cumsum(clust_lens))

        ##
        seqPattern::plotPatternDensityMap(
            regionsSeq = seqs,
            seqOrder = seq_order,
            flankUp = fl_up,
            flankDown = fl_down,
            patterns = motifs,
            addReferenceLine = TRUE,
            xTicksAt = use_xticksAt,
            xTicks = c(-x, 0, x),
            yTicksAt = use_yticksAt,
            yTicks = seq_along(use_yticksAt),
            useMulticore = parallelize,
            nrCores = n_cores,
            outFile = use_outFile,
            plotWidth = fwidth,
            plotHeight = fheight,
            ...)
        ##


        ##
        fname_patt <- basename(use_outFile)
        files <- list.files(use_tempdir, pattern = fname_patt,
            full.names = TRUE)
        ##
        legend_file <- grep("ColorLegend", files)
        final_legend_fname <- paste0(fname, ".ColorLegend", use_ext)
        file.copy(from = files[legend_file], to = final_legend_fname)

        ##
        ## order the files by the sequence in motifs TODO
        ## files <- c(files[-legend_file], files[legend_file])
        files <- files[-legend_file]
        ## For the moment, keep legend file separate
        stopifnot(length(files) > 0)

        files <- rev(files)
        #files <- c(files, clust_pl_fname)
        pick_pls <- lapply(seq_along(files), function(f){
            pl <- ggimage::as.ggplot(magick::image_read(files[f]))
        })


        grid_pl <- cowplot::plot_grid(plotlist = pick_pls,
                                        nrow = 1, ncol = length(files),
                                        align = c("hv")
        )

        grid_pl2 <- cowplot::plot_grid(
            ggimage::as.ggplot(magick::image_read(final_legend_fname)), grid_pl,
            scale = c(hm_scale_factor, 1), rel_widths = c(0.1, 1))

        ##
        cowplot::save_plot(filename = fname_w_ext, plot = grid_pl2)
        # cowplot::ggsave2(filename = fname_w_ext, plot = grid_pl2,
        #     scale = 1.5, width = 4*length(files), height = 2.5*length(files),
        #     units = "in")
        ##
        cli::cli_alert_success(paste0("Motif heatmaps written to: ",
            fname_w_ext))
        cli::cli_alert_success(paste0("Clustering legend in file:",
            clust_pl_fname))
        cli::cli_alert_success(paste0("Color legend in file:",
            final_legend_fname))
    })


}



