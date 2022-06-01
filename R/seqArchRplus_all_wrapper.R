## wrapper to call all together


#' @title Generate all plots for a given sample
#'
#' @param sname The sample name
#' @param bed_info_fname The BED filename with information on tag clusters. See
#' details for expected columns (column names)/information
#'
#' @param custom_colnames Specify custom column/header names to be used
#' with the BED file information
#'
#' @param seqArchR_clusts The seqArchR clusters' list
#'
#' @param raw_seqs The sequences corresponding to the cluster elements (also
#' available from the seqArchR result object)
#'
#' @param cager_obj The CAGEr object. This expects that
#' \code{\link[CAGEr]{clusterCTSS}} has been run beforehand. Default is NULL
#'
#' @param tc_gr The tag clusters as a \code{\link[GenomicRanges]{GRanges}}.
#' Default is NULL
#'
#' @param use_q_bound Logical. Write the lower and upper quantiles as tag
#' cluster boundaries in BED track files with tag clusters. Default is TRUE
#'
#' @param order_by_iqw Logical. Set TRUE to order clusters by the IQW median
#' or mean. Set argument `use_median_iqw`.to TRUE to use the median, else will
#' use mean if FALSE
#'
#' @param use_median_iqw Logical. Set TRUE if the median IQW for each cluster
#' is used to order the clusters. Otherwise, the mean IQW will be used when set
#' to FALSE
#'
#' @param iqw,tpm,cons Logical. Specify TRUE when the corresponding plots
#' should be included
#'
#' @param pos_lab The position labels
#'
#' @param txdb_obj The TranscriptsDB object
#'
#' @param org_name The organism name. This is used in the tracknames for tracks
#' writte nas BED files
#'
#' @param qLow,qUp Numeric values between 0 and 1. These are required when
#' cager_obj is provided instead of the tag clusters `tc_gr`
#'
#' @param tss_region For ChIPseeker, "region range of TSS"
#'
#' @param raw_seqs_mh Specify the sequences to be used for motif heatmaps, if
#' they are different than the sequences clustered by seqArchr. Default is NULL,
#' when `raw_seqs`` are used
#'
#' @param motifs Specify a character vector of motif words (using IUPAC
#' notation) to be visualized as a heatmap
#'
#' @param motif_heatmaps_flanks Specify a vector of different flank values to
#' be considered for visualization. Same size flanks are considered upstream
#' as well as downstream, hence one value suffices for eac hvisualization.
#' When a vector `c(50, 100, 200)` is specified, three motif heatmap files
#' (three separate PNG files) are created, each with one flank size. The
#' motif heatmap file will contain separate heatmaps for each of the specified
#' motifs in the `motifs` argument
#'
#'
#' @param dir_path The path to the directory where files are saved
#'
#' @param txt_size The text size to be used for the plots. This is
#' some high value because the plots are written with `dpi=300` and are often
#' large in size, especially the combined panel plots
#'
#' @details The expected columns (and column names) in the BED file are
#' "chr", "start", "end", "width", "strand", "score", "nr_ctss",
#' "dominant_ctss", "domTPM", "q_<qLow>", "q_<qUp>", "IQW", "tpm". Depending on
#' the values for arguments qLow and qUp, the corresponding column names are
#' formed. For example, if `qLow` and `qUp` are 0.1 and 0.9, the column names
#' are "q_0.1" and "q_0.9". These columns are mostly present by default in the
#' CAGEr tag clusters.
#'
#' The supplied clusters are ordered by their mean/median interquantile widths
#' before proceeding to generate the visualizations.
#'
#' @return A list holding generated plots; some which are directly written to
#' disk are not included in this list.
#'
#' The included plots are:
#'
#' - Boxplots of IQW (and TPM and conservation score when available)
#' distributions for each cluster (as a single combined plot)
#' - Annotation percentages per cluster as stacked barplots (as a single
#' combined plot)
#' - Annotation percentages per cluster as stacked barplots as a list
#' - Sequence logos of all cluster architectures (as a single combined plot)
#' - Sequence logos of all cluster architectures as a list
#' - Strand-separated sequence logos of all cluster architectures as a list
#' - Per cluster distribution of tag clusters on chromosomes and strands
#'
#' In addition, the following plots are written to disk:
#' - Visualization of all clustered sequences as a matrix
#' - Visualization of motif occurrences (for specified motifs) in all
#' clustered sequences
#'
#' In addition, the individual clusters from seqArchR are written to disk as
#' BED track files that can be viewed in the genome browser/IGV.
#'
#' @export
#'
#' @examples
#'
#' library(GenomicRanges)
#' library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' library(ChIPseeker)
#' library(Biostrings)
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' ## info_df <- read.delim(file = bed_fname,
#' ##         sep = "\t", header = TRUE)
#'
#'
#' tc_gr <- readRDS(system.file("extdata", "example_tc_gr.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#'
#'
#' raw_seqs <- Biostrings::readDNAStringSet(
#'                           filepath = system.file("extdata",
#'                             "example_promoters45.fa.gz",
#'                             package = "seqArchRplus",
#'                             mustWork = TRUE)
#'                         )
#'
#' raw_seqs_mh <- Biostrings::readDNAStringSet(
#'                           filepath = system.file("extdata",
#'                             "example_promoters200.fa.gz",
#'                             package = "seqArchRplus",
#'                             mustWork = TRUE)
#'                         )
#'
#'
#' generate_all_plots(sname = "sample1",
#'                    bed_info_fname = bed_fname,
#'                    seqArchR_clusts = use_clusts,
#'                    raw_seqs = raw_seqs,
#'                    tc_gr = tc_gr,
#'                    use_q_bound = FALSE,
#'                    order_by_iqw = TRUE,
#'                    use_median_iqw = TRUE,
#'                    iqw = TRUE, tpm = TRUE, cons = FALSE,
#'                    pos_lab = -45:45,
#'                    txdb_obj = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                    org_name = "Dmelanogaster22",
#'                    qLow = 0.1, qUp = 0.9,
#'                    tss_region = c(-500, 100),
#'                    raw_seqs_mh = raw_seqs_mh,
#'                    motifs = c("WW", "SS", "TATAA", "CG"),
#'                    motif_heatmaps_flanks = c(50, 100, 200),
#'                    dir_path = tempdir(),
#'                    txt_size = 25)
#'
#'
generate_all_plots <- function(sname, bed_info_fname,
                                custom_colnames = NULL,
                                seqArchR_clusts,
                                raw_seqs,
                                cager_obj = NULL,
                                tc_gr = NULL,
                                use_q_bound = FALSE,
                                order_by_iqw = TRUE,
                                use_median_iqw = TRUE,
                                iqw = TRUE, tpm = TRUE, cons = FALSE,
                                pos_lab = NULL,
                                txdb_obj = NULL,
                                org_name = NULL,
                                qLow = 0.1, qUp = 0.9,
                                tss_region = c(-500, 100),
                                raw_seqs_mh = NULL,
                                motifs = c("WW", "SS", "TATAA", "CG"),
                                motif_heatmaps_flanks = c(50, 100, 200),
                                dir_path,
                                txt_size = 25) {
    iqw_tpm_pl <-
        annotations_oneplot_pl <-
        annotations_list_pl <-
        seqlogos_oneplot_pl <-
        seqlogos_list_pl <-
        stranded_seqlogos_pl <-
        per_cl_strand_pl <- NULL
    ##
    # ## Prepare tc_gr
    # tc_gr2 <- .handle_tc_cager(tc_gr, cager_obj, sname, qLow, qUp)
    # stopifnot(!is.null(tc_gr2))
    # ## If tc_gr was prepared from a BED file, populate the clusts arg
    #
    # if(!is.null(tc_gr2[[2]]) && tc_gr2[[2]] == "bed"){
    #     clusts <- seq(length(tc_gr2[[1]]))
    # }
    #
    # ## clusts should be a list
    # if(!is.list(clusts)) clusts <- list(clusts)
    # ##
    # tc_gr <- tc_gr2[[1]]
    ##
    if (is.null(raw_seqs_mh)) {
        raw_seqs_mh <- raw_seqs
    }
    ##
    if(is.null(custom_colnames)){
        message("NO custom colnames")
        info_df <- utils::read.delim(file = bed_info_fname,
                                        sep = "\t", header = TRUE)
    }else{
        info_df <- utils::read.delim(file = bed_info_fname,
            sep = "\t", header = TRUE, col.names = custom_colnames
                )
    }
    # info_df <- utils::read.delim(
    #     file = bed_info_fname,
    #     sep = "\t", header = TRUE,
    #     col.names = c(
    #         "chr", "start", "end", "width",
    #         "strand", "score", "nr_ctss",
    #         "dominant_ctss", "domTPM",
    #         paste0("q_", c(qLow, qUp)),
    #         # "q_0.1", "q_0.9",
    #         "IQW", "tpm"
    #     )
    # )

    ##
    seqArchR_clusts_ord <- seqArchR_clusts
    if (order_by_iqw) {
        ## get clusters ordered by median IQW values
        seqArchR_clusts_ord <- order_clusters_iqw(
            sname = sname, clusts = seqArchR_clusts,
            info_df = info_df, order_by_median = TRUE
        )
    }
    ##
    use_clusts <- seqArchR_clusts_ord
    ##
    seqs_acgt_image(
        sname = sname,
        seqs = raw_seqs,
        seqs_ord = unlist(use_clusts),
        pos_lab = pos_lab, dir_path = dir_path
    )
    ##
    write_seqArchR_cluster_track_bed(
        sname = sname,
        clusts = use_clusts,
        info_df = info_df,
        use_q_bound = use_q_bound,
        one_zip_all = FALSE,
        org_name = org_name,
        dir_path = dir_path,
        include_in_report = FALSE,
        strand_sep = FALSE
    )
    ##
    iqw_tpm_pl <- iqw_tpm_plots(
        sname = sname,
        dir_path = dir_path,
        info_df = info_df,
        iqw = iqw, tpm = tpm, cons = cons,
        clusts = use_clusts,
        txt_size = txt_size
    )

    annotations_oneplot_pl <- per_cluster_annotations(
        sname = sname,
        clusts = use_clusts,
        tc_gr = tc_gr,
        cager_obj = NULL,
        qLow = qLow, qUp = qUp,
        txdb_obj = txdb_obj,
        tss_region = tss_region,
        orgdb_obj = NULL, dir_path = dir_path,
        one_plot = TRUE, txt_size = txt_size
    )

    annotations_list_pl <-
        per_cluster_annotations(
            sname = sname,
            clusts = use_clusts,
            tc_gr = tc_gr,
            cager_obj = NULL,
            qLow = qLow, qUp = qUp,
            txdb_obj = txdb_obj,
            tss_region = tss_region,
            orgdb_obj = NULL, dir_path = dir_path,
            one_plot = FALSE,
            txt_size = 12
        )

    seqlogos_oneplot_pl <-
        per_cluster_seqlogos(
            sname = sname,
            seqs = raw_seqs,
            clusts = use_clusts,
            pos_lab = pos_lab, bits_yax = "max",
            strand_sep = FALSE, one_plot = TRUE,
            dir_path = dir_path,
            txt_size = txt_size
        )
    ## for combining later
    seqlogos_list_pl <-
        per_cluster_seqlogos(
            sname = sname,
            seqs = raw_seqs,
            clusts = use_clusts,
            pos_lab = pos_lab, bits_yax = "max",
            strand_sep = FALSE, one_plot = FALSE,
            dir_path = dir_path,
            txt_size = 12
        )

    stranded_seqlogos_pl <-
        per_cluster_seqlogos(
            sname = sname,
            seqs = raw_seqs,
            clusts = use_clusts,
            pos_lab = pos_lab,
            bits_yax = "max",
            info_df = info_df,
            strand_sep = TRUE, one_plot = FALSE,
            dir_path = dir_path,
            txt_size = 12
        )

    plot_motif_heatmaps(
        sname = sname, seqs = raw_seqs_mh,
        flanks = motif_heatmaps_flanks,
        clusts = use_clusts,
        motifs = motifs,
        dir_path = dir_path,
        fheight = 800, fwidth = 1600
    )

    pair_colrs <- RColorBrewer::brewer.pal(n = 5, name = "Set3")
    per_cl_strand_pl <- per_cluster_strand_dist(
        sname = sname,
        clusts = use_clusts,
        info_df = info_df,
        dir_path = dir_path,
        colrs = pair_colrs[4:5]
    )

    return(list(
        iqw_tpm_pl,
        annotations_oneplot_pl,
        annotations_list_pl,
        seqlogos_oneplot_pl,
        seqlogos_list_pl,
        stranded_seqlogos_pl,
        per_cl_strand_pl
    ))
}
## =============================================================================
