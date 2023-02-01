## seqArchRplus_cluster_annotations

#' @title per_cluster_annotations
#'
#' @description This function helps annotate the genomic regions specified in
#' `tc_gr` with features, namely, promoter-TSS (transcription start site),
#' exons, 5'UTR, 3'UTR, introns and (distal) intergenic regions. This requires
#' that the annotations are available as a TxDb object. The selected genomic
#' regions can be specified as a single GenomicRanges object. These regions
#' can be specified directly as a BED file (when available) or select specific
#' regions from a larger set of regions based on some clustering.
#'
#' When working with CAGE data, if the CAGEr package was used and the
#' corresponding CAGEexp object is available, it can also be used -- see
#' `cager_obj` argument.
#'
#' @param sname Sample name. Default is NULL. This is a required argument
#' if the CAGEexp object is provided. See `cager_obj` argument
#' @param clusts List of sequence IDs in each cluster. This can be NULL only
#' when a BED file is passed to the argument `tc_gr`
#' @param tc_gr Tag clusters as \code{\link[GenomicRanges]{GRanges}} or a
#' BED file (specify filename with path). If `cager_obj` is not provided (i.e.,
#'  it is NULL), this argument is required. It will be ignored only if
#'  `cager_obj` is provided. Default is NULL
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
#' @details
#' When annotations for only selected clusters are required, alter
#' the `clusts` argument to specify only those selected clusters. Because
#' the `clusts` list holds the IDs of sequences belonging to each cluster, the
#' corresponding records are selected from the `tc_gr` GRanges object. This
#' approach requires that sequence IDs in `clusts` are directly associated with
#' the ranges in `tc_gr`. Also, see examples.
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
#' @examples
#'
#' ## Need the TxDb object to run these examples
#' if(require("TxDb.Dmelanogaster.UCSC.dm6.ensGene")){
#'
#' library(GenomicRanges)
#' library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' library(ChIPseeker)
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname,
#'          sep = "\t", header = TRUE)
#'
#' tc_gr_from_df <- GenomicRanges::makeGRangesFromDataFrame(info_df,
#'                                                   keep.extra.columns = TRUE)
#'
#' tc_gr <- readRDS(system.file("extdata", "example_tc_gr.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' tdir <- tempdir()
#'
#' # Get annotations for all clusters in use_clusts
#' annotations_pl <- per_cluster_annotations(sname = "sample1",
#'                          clusts = use_clusts,
#'                          tc_gr = tc_gr_from_df,
#'                          txdb_obj = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                          one_plot = FALSE,
#'                          dir_path = tdir,
#'                          tss_region = c(-500,100))
#'
#' # Get annotations for selected clusters in use_clusts
#' # -- First two clusters
#' selected_clusts <- lapply(seq(2), function(x) use_clusts[[x]])
#' # OR
#' # -- Mixed set of clusters, say 1 and 3 out of total 3
#' selected_clusts <- lapply(c(1,3), function(x) use_clusts[[x]])
#' #
#' annotations_pl <- per_cluster_annotations(sname = "sample1",
#'                          clusts = selected_clusts,
#'                          tc_gr = tc_gr,
#'                          txdb_obj = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                          one_plot = FALSE,
#'                          dir_path = tdir,
#'                          tss_region = c(-500,100))
#'
#' # Alternatively, you can also directly specify a BED file to the `tc_gr`
#' # argument. This is useful when one may not have access to the CAGEexp
#' # object, but only clusters' information is available in a BED file.
#' #
#'
#' annotations_pl <- per_cluster_annotations(sname = "sample1",
#'                          clusts = NULL,
#'                          tc_gr = bed_fname,
#'                          txdb_obj = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                          one_plot = FALSE,
#'                          dir_path = tdir,
#'                          tss_region = c(-500,100))
#' }
#'
#' @author Sarvesh Nikumbh
per_cluster_annotations <- function(sname = NULL, clusts = NULL,
                                    tc_gr = NULL,
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

    cli::cli_h1(paste0("All clusters' genomic annotations"))
    cli::cli_h2(paste0("Sample: ", sname))
    ## Check all needed arguments supplied


    ## Prepare tc_gr
    tc_gr2 <- .handle_tc_cager(tc_gr, cager_obj, sname, qLow, qUp)
    stopifnot(!is.null(tc_gr2))
    ## If tc_gr was prepared from a BED file, populate the clusts arg

    if(!is.null(tc_gr2[[2]]) && tc_gr2[[2]] == "bed"){
        clusts <- seq(length(tc_gr2[[1]]))
    }

    ## clusts should be a list
    if(!is.list(clusts)) clusts <- list(clusts)
    ##
    tc_gr <- tc_gr2[[1]]

    ##
    ## as many records in tc_gr as number of sequence IDs in clusts
    # if(length(tc_gr) != sum(lengths(clusts))){
    #     stop("Nb. of records in `tc_gr` should match the nb. of sequence IDs
    #         in `clusts`")
    # }
    ##
    parallelize <- FALSE
    if (n_cores > 1) parallelize <- TRUE
    if(parallelize)
        bpparam <- .handle_multicore(crs = n_cores, parallelize = parallelize)
    ##
    if(!is.null(dir_path)){
        result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
        fname <- file.path(result_dir_path,
                        paste0("Clusterwise_annotations.pdf"))
    }
    ##
    clust_labels <- .make_cluster_labels(clust = clusts, use_prefix, use_suffix)

    ## Without BiocParallel
    # clustwise_anno <- lapply(clusts, function(x) {
    #     foo_anno <- ChIPseeker::annotatePeak(tc_gr[x, ],
    #         tssRegion = tss_region,
    #         TxDb = txdb_obj,
    #         annoDb = orgdb_obj
    #     )
    #     foo_anno
    # })
    # names(clustwise_anno) <- seq(1, length(clusts))
    # return(clustwise_anno) ## to check
    # ## without dplyr?
    # df_list <- lapply(clustwise_anno, function(x) x@annoStat)

    ##
    peakAnno <- ChIPseeker::annotatePeak(tc_gr,
        tssRegion = tss_region,
        TxDb = txdb_obj,
        annoDb = orgdb_obj
    )
    ##

    clustwise_anno <- lapply(clusts, function(x){
        .get_anno_stat(peakAnno, x)
    })
    names(clustwise_anno) <- seq(1, length(clusts))
    ##

    colrs <- RColorBrewer::brewer.pal(n = 9, name = "Paired")
    ## Note: Generally, speed is not an issue when the number of clusters is
    ## around 40-50 or even up to 100.
    ## So we can use the base::rbind instead of the dplyr::bind_rows which is
    ## generally faster.
    ##

    sam <- do.call("rbind", clustwise_anno)
    per_df_nrows <- unlist(lapply(clustwise_anno, nrow))
    ## Add a column clust that we need downstream
    # sam$clust <- rep(seq_along(clustwise_anno), times = per_df_nrows)

    sam$clust <- rep(clust_labels, times = per_df_nrows)
    sam$clust <- factor(sam$clust, levels = rev(clust_labels))

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
        ind_print_plot <- ind_print_plot +
            ggplot2::labs(x = "Proportion", y = "Cluster")
        ##
        if(!is.null(dir_path)){
            ggplot2::ggsave(
                filename = fname, plot = ind_print_plot,
                device = "pdf", width = 15, height = 20, units = "in",
                dpi = 300
            )
        }
        return(ind_print_plot)
    } else {
        ## return individual plots as a list
        sam_split <- split(sam, f = sam$clust)

        ## Without BiocParallel
        annobar_list <- lapply(seq_along(sam_split), function(x) {
            ##
            anno_df <- sam_split[[x]]
            use_clust_label <- clust_labels[x]
            anno_df$Feature <- factor(anno_df$Feature,
                levels = levels(sam$Feature)
            )
            anno_pl <- .get_prop_anno_listplot(anno_df,
                txt_size = txt_size,
                colrs = colrs
            )
            anno_pl <- anno_pl +
                ggplot2::labs(x = "Proportion", y = use_clust_label)
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
        aes(y = clust, x = Frequency, fill = Feature)
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
            axis.title.y.left = element_text(size = txt_size),
            legend.position = "bottom",
            legend.text = element_text(size = txt_size)
        ) +
        ggplot2::scale_y_discrete(expand = expansion(add = c(0.75, 0))) +
        ggplot2::scale_x_continuous(expand = expansion(mult = 0.01)) +
        ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE)) +
        NULL
    ##
    clustwise_annobar
}
## =============================================================================

.get_anno_stat <- function(peakAnno, idx){

    anno_terms1 <- as.data.frame(peakAnno)[idx,]$annotation
    anno_terms_splits <- strsplit(anno_terms1, " ")
    anno_terms <- unlist(lapply(anno_terms_splits, function(x){
        if(is.na(x[1])){ "NA";
        }else if(grepl("Intron", x[1]) || grepl("Exon", x[1])){
            if(x[4] == 1){
                paste("1st", x[1])
            }else{
                paste("Other", x[1])
            }
        }else if(grepl("Promoter", x[1])){
            x[1]
        }else{ paste(x[1], x[2]) }
    }))
    ##
    wordFreqTable <- table(anno_terms)
    ##
    wordFreqDF <- data.frame(Feature = rownames(wordFreqTable),
        Frequency = as.vector(100*wordFreqTable/sum(wordFreqTable)))
    wordFreqDF
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
        # # ggplot2::labs(x = "Percentage (%)", y = "Cluster") +
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

