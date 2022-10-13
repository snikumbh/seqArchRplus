

#' @title Perform per cluster GO term enrichment analysis
#'
#' @description This function helps identify GO terms enriched per cluster.
#' This requires that the annotations are available as a TxDb object. The
#' selected genomic regions can be specified as a single GenomicRanges object.
#' These regions can be specified directly as a BED file (when available) or
#' select specific regions from a larger set of regions based on some
#' clustering.
#'
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
#' @param use_keytype Either of "ENTREZID" or "ENSEMBL". Required for use with
#' \code{\link[clusterProfiler]{enrichGO}}
#' @param one_file Logical. Default is TRUE. If set to FALSE the plots of
#' GO terms enriched per cluster are returned as a list, else all are
#' written to a single file as separate pages
#' @param bar_or_dot Specify "dot" for dotplot (default), or "bar" for barplot
#' @param dir_path Specify the /path/to/directory to store results
#' @param txt_size Adjust text size for the plots
#' @param n_cores For future use
#'
#' @details
#' Both `txdb_obj` and `orgdb_obj` are required.
#'
#'
#'
#' Per cluster, the enriched GO terms are visualized as a dot plot which shows
#' the enriched terms on the vertical axis and the ratio of genes that are
#' enriched for the given GO term vs. the total genes in the cluster.
#'
#'
#' @return
#' The list of dot plots showing GO term enrichments per cluster. This is a
#' list of ggplot2 plots.
#'
#' When `one_file` is set to TRUE (default), in addition to returning the list
#' of dot plots, these plots are also written to disk as a PDF, with one plot
#' per page.
#'
#' @importFrom clusterProfiler enrichGO dotplot
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics barplot
#'
#' @export
#'
#' @examples
#' library(GenomicRanges)
#' library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
#' library(ChIPseeker) ## important to load this package
#' library(org.Dm.eg.db)
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname, sep = "\t", header = TRUE)
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
#' # Get GO term enrichments for all clusters in use_clusts
#' go_pl <- per_cluster_go_term_enrichments(sname = "sample1",
#'                          clusts = use_clusts[1:2],
#'                          tc_gr = tc_gr_from_df,
#'                          txdb_obj = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                          dir_path = tdir,
#'                          one_file = FALSE,
#'                          tss_region = c(-500,100),
#'                          orgdb_obj = "org.Dm.eg.db")
#'
#'
#'
#' @author Sarvesh Nikumbh
per_cluster_go_term_enrichments <- function(sname = NULL, clusts = NULL,
                                    tc_gr = NULL,
                                    cager_obj = NULL,
                                    qLow = 0.1,
                                    qUp = 0.9,
                                    txdb_obj = NULL,
                                    tss_region = NULL,
                                    orgdb_obj = NULL,
                                    use_keytype = "ENTREZID",
                                    one_file = TRUE,
                                    bar_or_dot = "dot",
                                    dir_path = NULL,
                                    txt_size = 12,
                                    n_cores = 1){

    cli::cli_h1(paste0("All clusters' GO term enrichments"))
    cli::cli_h2(paste0("Sample: ", sname))
    ## Check all needed arguments supplied

    if(!is.null(dir_path)){
        result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
        fname <- file.path(result_dir_path,
            paste0("Clusterwise_GO_term_enrichments.pdf"))
    }

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

    # tc_gr <- .handle_tc_cager(tc_gr, cager_obj, sname, qLow, qUp)
    #
    peakAnno <- ChIPseeker::annotatePeak(tc_gr,
        tssRegion = tss_region,
        TxDb = txdb_obj,
        annoDb = orgdb_obj
    )

    go_pl <- lapply(seq_along(clusts), function(x){
        plot_title_text <- paste0("Arch: ", x,
                        " (n = ", length(clusts[[x]]), ")")
        pl <- .get_go_term_dot_plot(peakAnno, useOrgDb = orgdb_obj,
            useKeyType = use_keytype, choose_idx = clusts[[x]],
            plot_title_text = plot_title_text, bar_or_dot = bar_or_dot,
            font.size = txt_size)
    })

    if(!one_file) go_pl
    else{
        grDevices::pdf(file = fname, width = 10, height = 5, onefile = TRUE)
        lapply(go_pl, gridExtra::grid.arrange)
        grDevices::dev.off()
    }
    return(go_pl)
}


## peakAnno is the result obj obtained from ChIPseeker::annotatePeak function
## choose_idx is the indices from the df corresponding to the cluster
.get_go_term_dot_plot <- function(peakAnno, useOrgDb, useKeyType,
    choose_idx, plot_title_text,
    bar_or_dot = "dot", font.size = 12){
    # any_NAs <- which(is.na(as.data.frame(peakAnno)[choose_idx, ]$ENSEMBL))
    # plot_title_text <- paste(plot_title_text, length(any_NAs), "NAs")
    ##
    foo_go <- tryCatch(clusterProfiler::enrichGO(gene =
            as.data.frame(peakAnno)[choose_idx, useKeyType],
        keyType = useKeyType,
        OrgDb = useOrgDb, ont = "all"),
        error = function(e) e)

    if(is(foo_go, 'simpleError')){
        message("ontology ALL errored, hence using BP")
        foo_go <- tryCatch(clusterProfiler::enrichGO(gene =
                as.data.frame(peakAnno)[choose_idx, useKeyType],
            keyType = useKeyType,
            OrgDb = useOrgDb, ont = "BP"),
            error = function(e) e)
    }

    use_label <- "N/A"
    if(is(foo_go, 'simpleError')){
        message("ontology ALL and BP both errored, hence N/A")
        foo_go <- NULL
        use_label <- "All and BP errored, hence N/A"
    }

    ##
    if(!is.null(foo_go) && nrow(foo_go) > 0){
        if(bar_or_dot == "dot"){
            p1 <- clusterProfiler::dotplot(foo_go,
                title = plot_title_text,
                font.size = font.size
            )
        }else{
            p1 <- barplot(foo_go,
                title = plot_title_text,
                font.size = font.size
            )
        }

        p1 <- p1 + #DOSE::theme_dose(font.size=font.size) +
            ggplot2::scale_y_discrete(labels = scales::label_wrap(30)) +
            ggplot2::theme(title = element_text(size=18),
                legend.text = element_text(size=14),
                legend.title = element_text(size=14))
        return(p1)
    }else{
        message("Null/Zero rows")
        p1 <- ggplot() +
            ggplot2::geom_blank() + theme_bw() +
            theme(panel.grid = element_blank(),
                axis.title = element_blank(),
                axis.text = element_blank(),
                title = element_text(size=18)) +
            ggplot2::geom_text(aes(0,0,label=use_label)) +
            ggplot2::labs(title = plot_title_text)
        return(p1)
    }
}


