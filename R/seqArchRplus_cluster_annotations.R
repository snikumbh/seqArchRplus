## seqArchRplus_cluster_annotations

#' @title per_cluster_annotations
#'
#' @param sname Sample name
#' @param clusts List of sequence ids in each cluster.
#' @param tc_gr Tag clusters as \code{\link[GenomicRanges]{GRanges}}. If
#' `cager_obj` is not provided (is NULL), this argument is required.
#' It will be ignored only if cager_obj is provided. Default is NULL
#' @param cager_obj A CAGEexp object obtained from the CAGEr package, if and
#' when CAGEr was used to process the raw CAGE data
#' @param qLow,qUp The interquantile boundaries to be considered for obtaining
#' tag clusters from the CAGEexp object. See \code{\link[CAGEr]{tagClusters}}
#' @param txdb_obj A TxDb object storing transcript metadata
#' @param tss_region For ChIPseeker
#' @param orgdb_obj Organism-level annotation package
#' @param one_plot Logical. Default is TRUE. If set to FALSE the barplots of
#' annotations per cluster are returned as a list, else all are condensed into
#' single plot
#' @param dir_path Specify the /path/to/directory to store results
#' @param txt_size Adjust text size for the plots
#' @param use_suffix,use_prefix Character. Specify any suffix and/or prefix
#' you wish to add to the cluster labels
#' @param n_cores Numeric. If you wish to parallelize annotation of peaks,
#' specify the number of cores. Default is 1 (serial)
#'
#' @export
#'
#' @return
#' When `one_plot = TRUE`, a single plot where annotation barplots for each
#' cluster are put together (ordered as per the clusters in `clusts`).
#' Otherwise, a list of annotation barplots is returned (again ordered by
#' the clusters in `clusts`).
#'
#'
#' @importFrom BiocParallel bplapply MulticoreParam SerialParam SnowParam
#' multicoreWorkers ipclock ipcunlock ipcremove
#'
per_cluster_annotations <- function(sname, clusts, tc_gr = NULL,
                                    cager_obj = NULL,
                                    qLow = 0.1,
                                    qUp = 0.9,
                                    txdb_obj = NULL,
                                    tss_region = NULL,
                                    orgdb_obj = NULL,
                                    one_plot = TRUE,
                                    dir_path = NULL,
                                    txt_size = 12,
                                    use_suffix = NULL, use_prefix = "C",
                                    n_cores = 1) {
    message("SAMARTH--------")
    cli::cli_h1(paste0("All clusters' genomic annotations"))
    cli::cli_h2(paste0("Sample: ", sname))
    ## Check all needed arguments supplied
    tc_gr <- .handle_tc_cager(tc_gr, cager_obj, sname, qLow, qUp)
    stopifnot(!is.null(tc_gr))
    ## clusts should be a list
    if(!is.list(clusts)) clusts <- list(clusts)
    ## as many records in tc_gr as number of sequence IDs in clusts
    if(length(tc_gr) != sum(lengths(clusts))){
        stop("Nb. of records in `tc_gr` should match the nb. of sequence IDs
            in `clusts`")
    }
    ##
    parallelize <- FALSE
    if (n_cores > 1) parallelize <- TRUE
    if(parallelize)
        bpparam <- .handle_multicore(crs = n_cores, parallelize = parallelize)
    ##
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    fname <- file.path(result_dir_path, paste0("Clusterwise_annotations.pdf"))
    ##
    clust_labels <- make_cluster_labels(clust = clusts, use_prefix, use_suffix)


    ## Using BiocParallel
    ## Based on: https://support.bioconductor.org/p/110570/
    # id <- BiocParallel::ipcid()
    # clustwise_anno <- BiocParallel::bplapply(clusts, function(x, id) {
    #     BiocParallel::ipclock(id)
    #     foo_anno <- ChIPseeker::annotatePeak(tc_gr[x, ],
    #         tssRegion = tss_region,
    #         TxDb = txdb_obj,
    #         annoDb = orgdb_obj
    #     )
    #     BiocParallel::ipcunlock(id)
    #     foo_anno
    # }, id, BPPARAM = bpparam)
    # BiocParallel::ipcremove(id)

    ## Without BiocParallel
    clustwise_anno <- lapply(clusts, function(x) {
        foo_anno <- ChIPseeker::annotatePeak(tc_gr[x, ],
            tssRegion = tss_region,
            TxDb = txdb_obj,
            annoDb = orgdb_obj
        )
        foo_anno
    })

    names(clustwise_anno) <- seq(1, length(clusts))

    colrs <- RColorBrewer::brewer.pal(n = 9, name = "Paired")
    ## Note: Generally, speed is not an issue when the number of clusters is
    ## around 40-50 or even up to 100.
    ## So we can use the base::rbind instead of the dplyr::bind_rows which is
    ## generally faster.
    ##

    ## without dplyr?
    df_list <- lapply(clustwise_anno, function(x) x@annoStat)

    sam <- do.call("rbind", df_list)
    per_df_nrows <- unlist(lapply(df_list, nrow))
    ## Add a column clust that we need downstream
    sam$clust <- rep(seq_along(df_list), times = per_df_nrows)

    sam <- sam[, c(3, 1, 2)]
    rownames(sam) <- seq_len(sum(per_df_nrows))
    ##

    sam$Feature <- factor(sam$Feature,
        levels = c(
            "Promoter", "5' UTR", "1st Exon", "Other Exon",
            "1st Intron", "Other Intron", "Downstream (<=300)",
            "3' UTR", "Distal Intergenic"
        )
    )
    ##
    names(colrs) <- levels(sam$Feature)
    ##
    if (one_plot) {
        ## Solution using custom plotting
        clustwise_annobar <- .get_prop_anno_oneplot(sam,
            txt_size = txt_size,
            colrs = colrs
        )

        ## Put in the cluster names when printing independently
        ind_print_plot <- clustwise_annobar +
            ggplot2::scale_y_discrete(
                labels = rev(clust_labels),
                expand = expansion(add = c(0.5, 0.5))
            ) +
            ggplot2::theme(axis.text.y = element_text(
                size = txt_size,
                hjust = 0.85
            ))
        ##
        ggplot2::ggsave(
            filename = fname, plot = ind_print_plot,
            device = "pdf", width = 15, height = 20, units = "in",
            dpi = 300
        )

        return(clustwise_annobar)
    } else {
        ## return individual plots as a list
        sam_split <- split(sam, f = factor(sam$clust,
            levels = seq_along(clusts)
        ))

        ## Using BiocParallel
        # annobar_list <- BiocParallel::bplapply(sam_split, function(anno_df) {
        #     ##
        #     anno_df$Feature <- factor(anno_df$Feature,
        #         levels = levels(sam$Feature)
        #     )
        #     anno_pl <- .get_prop_anno_listplot(anno_df,
        #         txt_size = txt_size,
        #         colrs = colrs
        #     )
        # }, BPPARAM = bpparam)

        ## Without BiocParallel
        annobar_list <- lapply(sam_split, function(anno_df) {
            ##
            anno_df$Feature <- factor(anno_df$Feature,
                levels = levels(sam$Feature)
            )
            anno_pl <- .get_prop_anno_listplot(anno_df,
                txt_size = txt_size,
                colrs = colrs
            )
        })

        return(annobar_list)
    }
}
## =============================================================================


## anno_df holds information on the annoStats from the ChIPSeeker annotatePeak
## return object
## Expected columns are 'Frequency', 'Feature', 'clust'.
##
##
.get_prop_anno_oneplot <- function(anno_df, txt_size, colrs) {
    clustwise_annobar <- ggplot(
        anno_df,
        aes(
            y = forcats::fct_rev(
                stats::reorder(clust, as.numeric(clust))
            ),
            x = Frequency, fill = Feature
        )
    ) +
        geom_bar(
            stat = "identity", width = 0.6,
            position = position_fill(reverse = TRUE)
        ) +
        ggplot2::scale_fill_manual(name = "", values = colrs) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.text.y = element_blank(),
            axis.ticks.length.y.left = unit(0.1, units = "cm"),
            axis.ticks.length.x.bottom = unit(0.1, units = "cm"),
            axis.text.x.bottom = element_text(size = txt_size),
            axis.title.x.bottom = element_text(size = txt_size),
            legend.position = "bottom",
            legend.text = element_text(size = txt_size)
        ) +
        ggplot2::scale_y_discrete(expand = expansion(add = c(0.75, 0))) +
        ggplot2::scale_x_continuous(expand = expansion(mult = 0.01)) +
        ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE)) +
        ggplot2::ylab(NULL) +
        ggplot2::xlab("Percentage (%)")
    ##
    clustwise_annobar
}
## =============================================================================


.get_prop_anno_listplot <- function(anno_df, txt_size, colrs) {
    ##
    pl <- ggplot(anno_df, aes(y = clust, x = Frequency, fill = Feature)) +
        # ggplot2::theme_bw() +
        ggplot2::geom_bar(
            stat = "identity", width = 0.9,
            position = position_fill(reverse = TRUE)
        ) +
        ggplot2::scale_fill_manual(values = colrs) +
        ggplot2::scale_y_discrete(expand = expansion(add = c(0.1, 0.1))) +
        ggplot2::scale_x_continuous(expand = expansion(mult = 0.05)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title.y = element_text(size = txt_size),
            axis.title = element_text(size = txt_size),
            axis.text = element_text(size = txt_size),
            axis.text.y = element_text(hjust = 1.0, size = txt_size),
            axis.ticks.length.y.left = unit(0.2, units = "cm"),
            panel.grid = element_blank(),
            legend.text = element_text(size = txt_size),
            legend.title = element_text(size = txt_size),
            legend.position = "bottom",
            legend.direction = "horizontal"
        ) +
        ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE)) +
        ggplot2::labs(x = "Percentage (%)", y = "Cluster") +
        NULL
    return(pl)
}
## =============================================================================


.get_fixed_anno_ord <- function() {
    anno_terms_ord <- c(
        "Promoter", "5' UTR", "3' UTR", "1st Exon",
        "Other Exon", "1st Intron", "Intron",
        "Downstream",
        "Distal Intergenic"
    )

    anno_terms_ord
}
## =============================================================================


.get_named_colors <- function(anno_terms_ord, palname = "Set1") {
    use_colors <- .get_ncolors(n = length(anno_terms_ord), palname = palname)
    if (palname == "Set1") use_colors <- use_colors[c(2:length(use_colors), 1)]
    names(use_colors) <- anno_terms_ord
    use_colors
}
## =============================================================================
