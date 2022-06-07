## iqw-tpm-conservation


#' @title IQW, TPM plots
#'
#' @param sname Sample name
#'
#' @param dir_path The path to the directory where files will be written.
#'
#' @param info_df A DataFrame object holding information on the tag clusters.
#' Expected columns (names) are 'chr', 'start', 'end', 'IQW', 'domTPM', and
#' 'strand'.
#'
#' @param clusts List of sequence ids in each cluster.
#'
#' @param iqw Logical. Specify TRUE if boxplots of interquantile widths (IQW)
#' of the tag clusters corresponding to the promoters in each cluster
#' are to be plotted.
#'
#' @param tpm Logical. Specify TRUE if boxplots of TPM values for all clusters
#' are to be plotted.
#'
#' @param cons Logical. Specify TRUE if boxplots of conservation scores
#' (PhastCons scores) for all clusters. If this is TRUE, an additional column
#' named 'cons' is expected in the `info_df`.
#'
#' @param txt_size Specify text size to be used in the plots.
#'
#' @param wTitle Logical. If TRUE, the returned plot will contain a default
#'  title, which is the same as the filename. See details.
#'
#' @param use_suffix,use_prefix Character. Specify any suffix and/or prefix
#' you wish to add to the filename.
#'
#' @details The plots are written to a file named
#' "Sample_<sample_name>_IQW_TPM_Cons_plot.pdf" if all of `iqw`, `tpm`,
#' and `cons` are set to TRUE. This is also set as the plot title if `wTitle`
#' is set to TRUE.
#'
#' All plots are arranged by the IQWs (smallest on top, largest at the bottom),
#' even iff `iqw` is set to FALSE.
#'
#' @return The plot(s) as a ggplot2 object. The order of the plots is IQW,
#' followed by TPM values (when specified), followed by conservation scores
#' (when specified).
#'
#' @importFrom ggplot2 aes theme theme_void element_text element_blank
#' ggplot scale_y_continuous scale_y_discrete geom_boxplot dup_axis unit
#' ggtitle margin position_fill guide_legend geom_bar expansion xlab ylab
#' scale_x_log10 theme_bw
#' @importFrom cli cli_h1 cli_h1 cli_alert_warning cli_alert_info
#' @importFrom stats reorder
#' @importFrom forcats fct_reorder
#'
#' @export
#'
#' @examples
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname,
#'          sep = "\t", header = TRUE,
#'          col.names = c("chr", "start", "end", "width",
#'                  "dominant_ctss", "domTPM",
#'                  "strand",	"score", "nr_ctss",
#'                  "q_0.1", "q_0.9", "IQW", "tpm"))
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#'
#' iqw_tpm_pl <- iqw_tpm_plots(sname = "sample1",
#'                             dir_path = tempdir(),
#'                             info_df = info_df,
#'                             iqw = TRUE,
#'                             tpm = TRUE,
#'                             cons = FALSE,
#'                             clusts = use_clusts,
#'                             txt_size = 14)
#'
#'
#' @author Sarvesh Nikumbh
iqw_tpm_plots <- function(sname, dir_path, info_df, clusts, iqw = TRUE,
                            tpm = TRUE, cons = TRUE, txt_size = 12,
                            use_suffix = NULL, use_prefix = "C",
                            wTitle = TRUE) {
    cli::cli_h1(paste0("IQW-ordered boxplots"))
    cli::cli_h2(paste0("Sample: ", sname))

    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)

    info_df$clust_ID <- .get_clust_id_column(info_df, clusts)
    ##
    if (!any(iqw, tpm, cons)) {
        cli::cli_alert_warning(paste0(
            "None of `iqw`, `tpm` or `cons` are ",
            "set to TRUE. Setting `iqw` to `TRUE`"
        ))
        iqw <- TRUE
    }

    ##
    iqw_str <- ifelse(iqw, "_IQW", "")
    tpm_str <- ifelse(tpm, "_TPM", "")
    cons_str <- ifelse(cons, "_Cons", "")

    cli::cli_h2(paste("Plotting:", iqw_str, tpm_str, cons_str, "plots"))

    title_str <- paste0("Sample", sname, iqw_str, tpm_str, cons_str, "_plot")
    fname <- file.path(result_dir_path, paste0(title_str, ".pdf"))
    ##
    plots <- c(iqw, tpm, cons)
    which_plots <- which(plots)
    ## Can store as many plots as asked for (by setting TRUE)
    ## begin by all entries in plot_list to NULL, will later filter out NULLs
    plot_list <- vector("list", length(which_plots))
    set_yax_txt <- FALSE

    ##
    plot_list <- lapply(which_plots, function(x) {
        set_iqw <- ifelse(x == 1, TRUE, FALSE)
        set_tpm <- ifelse(x == 2, TRUE, FALSE)
        set_cons <- ifelse(x == 3, TRUE, FALSE)

        ## When we enter the loop for the first time, that is going to be for
        ## the first plot. This is when we need the y-axis text. From the
        ## second time onwards, there is a plot on the left to what is being
        ## drawn now. Therefore, no need for y-axis text.
        set_yax_txt <- ifelse(which(which_plots == x) == 1, TRUE, FALSE)
        ##
        temp_pl <- .get_iqw_ord_plot(
            iqw = set_iqw, tpm = set_tpm,
            phast = set_cons,
            y_axis_text = set_yax_txt,
            use_notch = FALSE, txt_size = txt_size,
            info_df = info_df, seqs_clust = clusts,
            use_suffix = use_suffix
        )
        temp_pl
    })

    ## Keep only non-null entries
    plot_list <- lapply(plot_list, function(x) if (!is.null(x)) x)

    ## IQW ordered plots
    if (length(plot_list) > 1) {
        rel_width_vals <- c(1, rep(0.85, times = length(plot_list) - 1))
    } else {
        rel_width_vals <- 1
    }

    comb_pl <- cowplot::plot_grid(
        plotlist = plot_list, nrow = 1,
        ncol = length(plot_list), align = "h",
        rel_widths = rel_width_vals
    )
    ##
    pl_w_title <- comb_pl + ggplot2::ggtitle(label = title_str)
    ##
    if (!is.null(fname)) {
        cowplot::save_plot(
            filename = fname, plot = pl_w_title,
            base_height = 14, base_width = 5 * length(plot_list)
        )
    }
    ##
    if (wTitle) {
        return(pl_w_title)
    } else {
        return(comb_pl)
    }
}
## =============================================================================


.get_iqw_ord_plot <- function(iqw = FALSE, tpm = FALSE, phast = FALSE,
                                use_notch = FALSE, y_axis_text = TRUE,
                                txt_size = 12, info_df, seqs_clust,
                                use_suffix = "X", use_prefix = "C") {
    clust_labels <- make_cluster_labels(
        clust = seqs_clust,
        use_prefix = use_prefix,
        use_suffix = use_suffix
    )
    clr <- RColorBrewer::brewer.pal(3, "Dark2")
    ##
    if (iqw) {
        pl <- ggplot(
            info_df,
            aes(
                y = forcats::fct_reorder(clust_ID, IQW,
                    .fun = median,
                    .desc = TRUE
                ),
                x = IQW
            )
        ) +
            ggplot2::geom_boxplot(
                outlier.size = 1, width = 0.5,
                notch = use_notch, color = "black",
                fill = clr[1]
            ) +
            xlab("IQW")
    }
    ##
    if (tpm) {
        pl <- ggplot(
            info_df,
            aes(
                y = forcats::fct_reorder(clust_ID, IQW,
                    .fun = median,
                    .desc = TRUE
                ),
                x = domTPM
            )
        ) +
            geom_boxplot(
                outlier.size = 1, width = 0.5, notch = use_notch,
                color = "black", fill = clr[2]
            ) +
            xlab("TPM") +
            ylab("Clusters")
    }
    if (phast) {
        pl <- ggplot(
            info_df,
            aes(
                y = forcats::fct_reorder(clust_ID, IQW,
                    .fun = median,
                    .desc = TRUE
                ),
                x = phast
            )
        ) +
            geom_boxplot(
                outlier.size = 1, width = 0.5, notch = use_notch,
                color = "black", fill = clr[3]
            ) +
            xlab("PhastCons score") +
            ylab("Clusters")
    }
    ##
    if (tpm || iqw) {
        pl <- pl + scale_x_log10() + ggplot2::annotation_logticks(
            sides = "b",
            short = unit(0.1, "cm"),
            mid = unit(0.2, "cm"),
            long = unit(0.3, "cm")
        )
    }
    pl <- pl + ylab("Clusters") +
        scale_y_discrete(
            labels = rev(clust_labels),
            expand = expansion(add = c(0.55, 0.55))
        ) +
        theme_bw() +
        theme(
            axis.title = element_text(colour = "black", size = txt_size),
            axis.text.x = element_text(
                colour = "black", angle = 0,
                size = txt_size,
                vjust = 0.5, hjust = 0.5
            ),
            axis.text.y = element_text(
                colour = "black",
                size = txt_size, hjust = 1
            )
        )
    ##

    if (!y_axis_text) {
        pl <- pl + theme(
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
        )
    }

    return(pl)
}
## =============================================================================


#' @title Order clusters by median or mean interquantile widths
#'
#' @param sname The sample name
#' @param clusts List of sequence ids in each cluster.
#' @param info_df The data.frame with all tag clusters information. The
#' following columns are expected in the data.frame:"chr", "start", "end",
#' "width", "strand", "score", "nr_ctss", "dominant_ctss", "domTPM",
#' "IQW", "tpm" and two additional columns based on qLow and qUp used.
#'
#' @param order_by_median Logical. Whether to order to clusters by their
#' median (when TRUE) or mean (when FALSE) interquantile widths.
#'
#' @return
#' The list of clusters ordered by their mean/median interquantile widths
#' (shortest first).
#'
#' @importFrom stats median
#'
#' @export
#'
#' @examples
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname,
#'          sep = "\t", header = TRUE,
#'          col.names = c("chr", "start", "end", "width",
#'                  "dominant_ctss", "domTPM",
#'                  "strand",	"score", "nr_ctss",
#'                  "q_0.1", "q_0.9", "IQW", "tpm"))
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#'
#' ordered_clusts <- seqArchRplus::order_clusters_iqw(
#'                                  sname = "sample1",
#'                                  clusts = use_clusts,
#'                                  info_df = info_df,
#'                                  order_by_median = TRUE)
#'
order_clusters_iqw <- function(sname, clusts, info_df,
                                order_by_median = TRUE) {
    cli::cli_h1(paste0("Order clusters by IQW"))
    cli::cli_h2(paste0("Sample: ", sname))
    cluster_medians_IQW <- unlist(lapply(clusts, function(x) {
        stats::median(info_df$IQW[x])
    }))
    cluster_means_IQW <- unlist(lapply(clusts, function(x) {
        base::mean(info_df$IQW[x])
    }))

    if (order_by_median) {
        ascending_order_IQW <- sort(cluster_medians_IQW,
            decreasing = FALSE,
            index.return = TRUE
        )
    } else {
        ascending_order_IQW <- sort(cluster_means_IQW,
            decreasing = FALSE,
            index.return = TRUE
        )
    }
    ##
    clusts_list_ordered <- lapply(ascending_order_IQW$ix, function(x) {
        clusts[[x]]
    })

    clusts_list_ordered
}
## =============================================================================
