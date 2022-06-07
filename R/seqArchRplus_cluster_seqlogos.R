## seqArchRplus_cluster_seqlogos


#' @title Plot per cluster sequence logos
#'
#' @param sname The sample name
#' @param seqs The raw sequences as a DNAStringSet. These are also available
#' as part of the seqArchR result object.
#' @param clusts List of sequence ids in each cluster.
#' @param pos_lab The position labels
#' @param bits_yax The yaxis limits. See argument in seqArchR
#' @param strand_sep Logical. Whether sequences are to be separated by strand.
#' @param one_plot Logical. Whether all sequence logos should be combined
#' into one grid (with ncol = 1)?
#'
#' @param info_df The information data.frame
#'
#' @param txt_size The font size for text
#'
#' @param save_png Logical. Set TRUE if you would like to save the
#' architectures sequence logos as PNG files
#'
#' @param dir_path The path to the directory where plot will be saved.
#'
#'
#' @details Plots the sequence logos of all clusters
#'
#' @return If `one_plot` is TRUE, one plot as a grid of all sequence logos is
#' returned.
#'
#' If `one_plot` is FALSE, a set of ggplot2-based sequence logos for all
#' clusters is saved to disk with the default filename
#' `Architectures_0-max.pdf`. This is a multi-page PDF document with the
#' sequence logo for each cluster on a separate page. Also, the list of plots
#' is returned.
#'
#' @importFrom seqArchR plot_ggseqlogo_of_seqs plot_arch_for_clusters
#' collate_seqArchR_result get_seqs_clust_list seqs_str collate_clusters
#' viz_seqs_acgt_mat
#'
#' @export
#'
#' @examples
#'
#' library(Biostrings)
#' raw_seqs <- Biostrings::readDNAStringSet(
#'                           filepath = system.file("extdata",
#'                             "example_promoters45.fa.gz",
#'                             package = "seqArchRplus",
#'                             mustWork = TRUE)
#'                         )
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' seqlogo_pl <- per_cluster_seqlogos(sname = "sample1",
#'                                    seqs = raw_seqs,
#'                                    clusts = use_clusts,
#'                                    pos_lab = -45:45,
#'                                    bits_yax = "max",
#'                                    strand_sep = FALSE,
#'                                    one_plot = TRUE,
#'                                    dir_path = tempdir(),
#'                                    txt_size = 14)
#'
#' @author Sarvesh Nikumbh
per_cluster_seqlogos <- function(sname, seqs = NULL, clusts,
                                    pos_lab = NULL, bits_yax = "max",
                                    strand_sep = FALSE, one_plot = TRUE,
                                    info_df = NULL, txt_size = 12,
                                    save_png = FALSE, dir_path) {
    if (strand_sep) {
        if (is.null(info_df)) {
            stop(
                "`info_df` cannot be NULL when ",
                "`strand_sep` is TRUE"
            )
        }
        strand_sep_pl <- strand_sep_seqlogos(
            sname = sname, seqs = seqs,
            clusts = clusts, info_df = info_df,
            pos_lab = pos_lab, bits_yax = bits_yax,
            dir_path = dir_path, txt_size = txt_size,
            save_png = save_png
        )
        return(strand_sep_pl)
    }
    cli::cli_h1(paste0("All clusters' sequence logos"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##
    message("Generating architectures for clusters of sequences...")
    fname <- file.path(result_dir_path, paste0("Architectures_0-max.pdf"))
    ##
    T_F_titles <- TRUE
    if (one_plot) T_F_titles <- FALSE
    ## arch_list plots directly using seqArchR function
    arch_pl <- seqArchR::plot_arch_for_clusters(
        seqs = seqs,
        clust_list = clusts,
        pos_lab = pos_lab,
        set_titles = T_F_titles,
        method = "bits",
        pdf_name = NULL, show = FALSE
    )

    arch_list <- lapply(seq_along(arch_pl), function(y) {
        pl <- arch_pl[[y]] + ggplot2::theme(
            axis.text = ggplot2::element_text(size = txt_size),
            axis.text.y = ggplot2::element_text(vjust = 0.5),
            axis.title.y = ggplot2::element_text(size = txt_size),
            axis.ticks.length = ggplot2::unit(0.04, "cm"),
            plot.margin = ggplot2::unit(c(0.2, 0.1, -0.4, 0.1), "cm")
        ) +
            ggplot2::scale_y_continuous(
                sec.axis = ggplot2::dup_axis(
                    name = paste0("C", y),
                    labels = NULL
                )
            )
        pl
    })
    ## Making a plot_grid and plotting is tedious; the PDF height has to
    ## be set by trial and error to a very large value. Instead, we decide
    ## to plot this a pair on one page
    if (one_plot) {
        comb_one_pl <- cowplot::plot_grid(plotlist = arch_list, ncol = 1)
        return(comb_one_pl)
    } else {
        grDevices::pdf(file = fname, width = 20, height = 2, onefile = TRUE)
        lapply(arch_list, print)
        grDevices::dev.off()

        ## save PNGs
        if (save_png) {
            .save_PNGs(
                dir_path = result_dir_path, plot_list = arch_list,
                txt_size = txt_size
            )
        }
        ##
        return(arch_list)
    }
    ##
}
## =============================================================================

.save_PNGs <- function(dir_path, plot_list, strand_sep = FALSE, txt_size = 10) {
    use_res_path <- file.path(dir_path, "arch_png")
    if (strand_sep) use_res_path <- file.path(dir_path, "strand_sep_arch_png")

    use_ht <- 3
    if (strand_sep) use_ht <- 6

    stopifnot(.check_and_create_dir(use_res_path))
    for (p in seq_along(plot_list)) {
        fname <- file.path(
            use_res_path,
            paste0("Architecture_clust", p, "_0-max.png")
        )
        pl <- plot_list[[p]] +
            ggplot2::theme(
                axis.text = ggplot2::element_text(size = txt_size),
                axis.title.y = ggplot2::element_text(size = txt_size),
                axis.text.y = ggplot2::element_text(size = txt_size)
            )
        cowplot::ggsave2(fname,
            plot = pl, device = "svg",
            width = 20, height = use_ht, units = "cm", dpi = 300
        )
    }
}
## =============================================================================





strand_sep_seqlogos <- function(sname, seqs, clusts, info_df, pos_lab,
                                bits_yax, dir_path, txt_size = 12,
                                save_png = FALSE) {
    cli::cli_h1(paste0("All clusters' strand-separated sequence logos"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##

    ## STRAND SEPARTATED SEQUENCE LOGOS (AUTO LIMITS)
    fname <- file.path(
        result_dir_path,
        paste0("Architectures_0-max_strand_separated.pdf")
    )


    per_cl_idx_by_strand <- lapply(clusts, function(x) {
        idx_p_n <- lapply(c("+", "-"), function(y) {
            .get_strand_specific_indices(
                df = info_df,
                seq_ids_in_cl = x,
                strand_val = y
            )
        })
        names(idx_p_n) <- c("+", "-")
        idx_p_n
    })


    per_cl_title_by_strand <- lapply(seq_along(clusts), function(x) {
        title_p_n <- lapply(c("+", "-"), function(y) {
            .get_strand_plot_title(x,
                nclust = length(clusts),
                this_n = length(per_cl_idx_by_strand[[x]][[y]]),
                tot_n = length(clusts[[x]]),
                strand_val = y
            )
        })
        names(title_p_n) <- c("+", "-")
        title_p_n
    })

    plots_p_n <- lapply(seq_along(clusts), function(x) {
        if (length(per_cl_idx_by_strand[[x]]$`+`) < 1) {
            print("plus is zero")
            pl_p <- ggplot2::ggplot() +
                ggplot2::theme_void()
        } else {
            pl_p <- seqArchR::plot_ggseqlogo_of_seqs(
                seqs = seqs[per_cl_idx_by_strand[[x]]$`+`],
                pos_lab = pos_lab,
                xt_freq = 5,
                title = per_cl_title_by_strand[[x]]$`+`,
                bits_yax = bits_yax
            )
            pl_p <- pl_p + theme(
                axis.text = element_text(size = txt_size),
                axis.text.y = element_text(vjust = 0.5),
                axis.title.y = element_text(size = txt_size),
                axis.ticks.length = unit(0.04, "cm"),
                plot.margin = unit(c(0.2, 0.1, -0.4, 0.1), "cm")
            ) +
                ggplot2::scale_y_continuous(
                    sec.axis = dup_axis(
                        name = paste0("C", x),
                        labels = NULL
                    )
                )
        }
        if (length(per_cl_idx_by_strand[[x]]$`-`) < 1) {
            print("negative is zero")
            pl_n <- ggplot() +
                theme_void()
        } else {
            pl_n <- seqArchR::plot_ggseqlogo_of_seqs(
                seqs = seqs[per_cl_idx_by_strand[[x]]$`-`],
                pos_lab = pos_lab,
                xt_freq = 5,
                title = per_cl_title_by_strand[[x]]$`-`,
                bits_yax = bits_yax
            )
            pl_n <- pl_n + theme(
                axis.text = element_text(size = txt_size),
                axis.text.y = element_text(vjust = 0.5),
                axis.title.y = element_text(size = txt_size),
                axis.ticks.length = unit(0.04, "cm"),
                plot.margin = unit(c(0.2, 0.1, -0.4, 0.1), "cm")
            ) +
                ggplot2::scale_y_continuous(
                    sec.axis = dup_axis(
                        name = paste0("C", x),
                        labels = NULL
                    )
                )
        }
        ##
        pl_p
        ##
        comb_pl_p_n <- cowplot::plot_grid(pl_p, pl_n, nrow = 2, align = "v")
        comb_pl_p_n
    })
    ## Making a plot_grid and plotting is tedious; the PDF height has to be set
    ## by trial and error to a very large value. Instead, we decide to plot this
    ## a pair on one page
    grDevices::pdf(file = fname, width = 20, height = 4.5, onefile = TRUE)
    lapply(plots_p_n, print)
    grDevices::dev.off()

    ## save PNGs
    if (save_png) {
        .save_PNGs(
            dir_path = result_dir_path, plot_list = plots_p_n,
            strand_sep = TRUE, txt_size = 10
        )
    }


    return(plots_p_n)
}
## =============================================================================
