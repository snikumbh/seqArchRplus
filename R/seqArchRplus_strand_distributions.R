## strand_distributions

#' @title Visualize as barplots how promoters in each cluster are distributed
#' on different chromosomes and strands
#'
#' @param sname The sample name
#' @param clusts List of sequence ids in each cluster.
#' @param info_df The information data.frame
#' @param dir_path Specify the path to the directory on disk where plots
#' will be saved
#' @param colrs Specify colors used for two strands. By default, the fourth
#' and fifth color from the palette 'Set3' are used
#' @param txt_size Size of the text in the plots (includes plot title,
#' y-axis title, axis texts, legend title and text)
#' @param fwidth Width of the individual plots in file
#' @param fheight Height of the plots in file
#'
#' @return A list of plots showing the per cluster division of promoters on
#' chromosomes and strands. These plots are also written to disk in file named
#' "Per_cluster_strand_distributions.pdf"
#'
#' @importFrom ggplot2 theme_classic scale_fill_manual element_line
#'
#' @export
#'
#' @examples
#'
#' library(RColorBrewer)
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname,
#'          sep = "\t", header = TRUE,
#'          col.names = c("chr", "start", "end", "width",
#'                  "strand", "score", "nr_ctss",
#'                  "dominant_ctss", "domTPM",
#'                  "q_0.1", "q_0.9", "IQW", "tpm"))
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' pair_colrs <- RColorBrewer::brewer.pal(n = 5, name = "Set3")[4:5]
#' per_cl_strand_pl <- per_cluster_strand_dist(sname = "sample1",
#'                                             clusts = use_clusts,
#'                                             info_df = info_df,
#'                                             dir_path = tempdir(),
#'                                             colrs = pair_colrs
#'                                             )
#'
#' @author Sarvesh Nikumbh
per_cluster_strand_dist <- function(sname, clusts, info_df, dir_path = NULL,
                                    colrs = c("#FB8072", "#80B1D3"),
                                    txt_size = 14,
                                    fwidth = 20, fheight = 3) {
    cli::cli_h1(paste0("All clusters' strand distributions"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    if(!is.null(dir_path)){
        result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
        fname <- file.path(
            result_dir_path,
            paste0("Per_cluster_strand_distributions.pdf")
        )
    }


    chr_set <- unique(info_df$chr)
    nChr <- length(chr_set)
    empty_df <- data.frame(
        "Chromosomes" = rep(chr_set, each = 2),
        "Strand" = rep(c("-", "+"), times = nChr),
        "Frequency" = rep(0, times = nChr)
    )
    ## Using 'times' instead of 'each' for Strand column and
    ## 'each' instead of 'times in Chromosomes column
    ## in empty_df makes
    ## it easier when iterating to pick each chromosome for pie_list

    ## Used later when + and - frequencies are zero
    # dummy_df <- empty_df
    # dummy_df[, "Frequency"] <- 5
    # ##
    chrStr <- paste0(empty_df$Chromosomes, empty_df$Strand)

    plot_list <- lapply(seq_along(clusts), function(x) {
        df_strand <- data.frame(
            "chr" = info_df$chr[clusts[[x]]],
            "strand" = info_df$strand[clusts[[x]]]
        )
        df_strand_tab <- table(df_strand)
        df_strand_df <- data.frame(df_strand_tab)
        RchrStr <- paste0(df_strand_df$chr, df_strand_df$strand)

        filled_df <- empty_df

        filled_df$chrStr <- chrStr
        fill_idx <- unlist(lapply(RchrStr, function(a) {
            which(filled_df$chrStr == a)
        }))

        filled_df$Frequency[fill_idx] <- df_strand_df$Freq

        ## Barplots
        pl <- ggplot(filled_df) +
            ggplot2::geom_col(aes(
                x = Chromosomes, y = Frequency,
                fill = Strand
            ),
            position = ggplot2::position_dodge(width = 0.5)
            ) +
            scale_fill_manual(values = colrs) +
            xlab("") +
            ggtitle(label = .get_strand_plot_title(
                this_id = x,
                nclust = length(clusts),
                tot_n = length(clusts[[x]]),
                strand_val = NULL
            )) +
            theme_classic() +
            theme(
                plot.title = element_text(size = txt_size),
                axis.text = element_text(size = txt_size),
                axis.title.y = element_text(size = txt_size),
                legend.title = element_text(size = txt_size),
                legend.text = element_text(
                    size = txt_size,
                    face = "plain",
                    color = "black"
                ),
                legend.position = "right",
                panel.grid.major.y = element_line()
            )
        return(pl)
    })

    if(!is.null(dir_path)){
        grDevices::pdf(file = fname, width = fwidth, height = fheight,
            onefile = TRUE)
        lapply(plot_list, gridExtra::grid.arrange)
        grDevices::dev.off()
    }
    return(plot_list)
}
## =============================================================================
