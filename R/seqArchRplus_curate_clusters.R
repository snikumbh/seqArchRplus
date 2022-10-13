

#' @title Curate clusters from seqArchR result
#'
#' @description seqArchR result stores the clusters obtained at every iteration.
#' It is possible that the previously chosen agglomeration and/or distance
#' method used with hierarchical clustering does not yield reasonable
#' clusters. This function enables minor curation of the clusters obtained
#' from the hierarchical clustering step.
#'
#' @param sname The sample name
#' @param use_aggl The agglomeration method to be used for hierarchical
#' clustering. This is passed on to `seqArchR::collate_seqArchR_result`.
#' See argument `aggl_method` in \link[seqArchR]{collate_seqArchR_result}.
#' @param use_dist The distance method to be used for hierarchical clustering.
#' This is passed on to \link[seqArchR]{collate_seqArchR_result}. See argument
#' `dist_method` in \link[seqArchR]{collate_seqArchR_result}.
#' @param seqArchR_result The seqArchR result object.
#' @param iter Specify which iteration of the seqArchR result should be used
#' for obtaining clusters.
#' @param pos_lab The position labels
#' @param regularize Logical. Specify TRUE if the basis vector comparison
#' is to be regularized. Requires you to set `topn` which is set to 50 as
#' default. See argument `regularize` in
#' \link[seqArchR]{collate_seqArchR_result}.
#' @param topn Numeric. The top N features (nucleotide-position pairs)
#' that will be used for distance computation, rest will be ignored.
#' See argument `topn` in \link[seqArchR]{collate_seqArchR_result}.
#' @param use_cutk Value of K (number of clusters) for cutting the hierarchical
#' clustering tree.
#' @param need_change A list of elements (cluster IDs in the input clusters)
#' that need re-assignment. Elements
#' @param change_to A list of elements (cluster IDs in the input clusters)
#' to be assigned to those that need re-assignment. In case there is a candidate
#' that needs to be put into a new, independent cluster of itself, use
#' 0 (as numeric). Both `need_change` and `change_to` should be empty lists if
#' no re-assignment is to be performed.
#'
#' @param final Logical, set to TRUE or FALSE
#'
#' @param dir_path The /path/to/the/directory where files will be written.
#' Default is NULL.
#'
#' @details
#' This function helps the user work through the curation in at most three
#' steps.
#'
#' 1. This function performs hierarchical clustering on the seqArchR clusters
#' (of the specified iteration) to obtain a clustering result.
#' The resulting clustering is visualized by a dendrogram, color-coded cluster
#' assignments, and corresponding sequence logos. Using this visualization,
#' the user can identify/estimate the (nearly right) number of clusters to cut
#' the dendrogram. The **first call** uses K = 2.
#'
#' 2. Visually examine and count the tentative number of clusters (K) deemed
#' right. Because these architectures are now arranged by the hierarchical
#' clustering order, identifying this tentative value of K is much easier.
#' Call the function this **second** time with the identified value of K.
#' Look at the visualization now generated to determine if it is good enough,
#' i.e., it requires only minor re-assignments.
#'
#' 3. Identify cases of cluster assignments that you wish to re-assign to
#' different clusters. These can be noted as a list and supplied in a
#' subsequent call to the function.
#'
#' 4. In the **final call** to the function, set `final = TRUE`, supply the re-
#' assignments as two lists `need_change` and `change_to`.
#'
#' More on re-assignments using arguments \code{need_change} and
#' \code{change_to}: If any element is to be put into a new cluster, use a
#' numeral 0 in change_to. This can be done in either scenario: when any
#' element is re-assigned as a singleton cluster of itself, or as clustered
#' with some other element coming from some other existing cluster. Consider,
#' for instance, among some 33 clusters identified by seqArchR, the following
#' re-assignments are executed.
#'
#' Original_clustering <- list(c(1), c(2,3), c(4,5), c(6, 7))
#' need_change <- list(c(5), c(2), c(3, 6))
#' change_to <- list(1, 0, 0)
#'
#' In the above, element 5 is re-assigned to the cluster containing element 1.
#' Element 2 is re-assigned to a new, singleton cluster of itself, while
#' elements 3 and 6 (which initially can belong to same/any two clusters) are
#' collated together. Note that it is important to use \code{c()}.
#'
#' Also see examples below.
#'
#'
#' @return
#' This function returns a list holding (a) 'curation_plot': plot showing the
#' dendrogram + sequence logos, and (b) 'clusters_list': the sequence clusters
#' as a list.
#'
#' to help
#' perform curation and document it.
#'
#'
#' When `final = FALSE`, the `curation_plot' shows the dendrogram + sequence
#' logos of clusters (ordered by \code{hclust} ordering).
#' The 'clusters_list' holds the hclust ordered clusters.
#' If the `dir_path`
#' is specified, a PDF file showing the same figure is also written at the
#' location using the default filename
#' `<Sample_name>_dend_arch_list_reg_top50_euclid_complete_<K>clusters.pdf`.
#'
#' When `final = TRUE`, 'clusters_list' holds the clusters with
#' re-assignments executed, and 'curation_plot' of dendrogram + sequence logos
#' now has an additional panel showing the sequence logos upon collation of
#' clusters as specified by the re-assignments.
#' Also the cluster IDs in the dendrograms have colors showing the
#' re-assignments, i.e., elements that were re-assigned to different clusters,
#' also have appropriate color changes reflected.
#'
#' When `dir_path` is provided, the curation plot is written to disk using
#' the same filename as before except for suffix 'final' attched to it.
#'
#'
#'
#' @importFrom stats cutree median
#'
#' @export
#'
#' @examples
#'
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
#' seqArchR_result <- readRDS(system.file("extdata", "seqArchR_result.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' use_aggl <- "complete"
#' use_dist <- "euclid"
#'
#' ## get seqArchR clusters custom curated
#' seqArchR_clusts <- seqArchRplus::curate_clusters(sname = "sample1",
#'     use_aggl = use_aggl, use_dist = use_dist,
#'     seqArchR_result = seqArchR_result, iter = 5,
#'     pos_lab = NULL, regularize = FALSE, topn = 50,
#'     use_cutk = 5, final = FALSE, dir_path = tempdir())
#'
#' ## Form the lists need_change and change_to for re-assignments
#' need_change <- list(c(2))
#' change_to <- list(c(1))
#'
#' ## This fetches us clusters with custom/curated collation in _arbitrary_
#' ## order. See the next function call to order_clusters_iqw that orders
#' ## these clusters by their median/mean IQW
#'
#' seqArchR_clusts <- seqArchRplus::curate_clusters(sname = "sample1",
#'     use_aggl = use_aggl, use_dist = use_dist,
#'     seqArchR_result = seqArchR_result, iter = 5,
#'     pos_lab = NULL, regularize = FALSE, topn = 50,
#'     use_cutk = 5,
#'     need_change = need_change, change_to = change_to,
#'     final = TRUE, dir_path = tempdir())
#'
#'
#'
#' @author Sarvesh Nikumbh
curate_clusters <- function(sname, use_aggl = "ward.D", use_dist = "euclid",
                            seqArchR_result, iter, pos_lab = NULL,
                            regularize = TRUE, topn = 50, use_cutk = 2,
                            need_change = NULL, change_to = NULL,
                            final = FALSE, dir_path = NULL) {
    cli::cli_h1(paste0("seqArchR result clusters curation"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    fname_suffix <- ""
    if (final) {
        if (is.null(need_change) || is.null(change_to)) {
            stop("Both `need_change` and `change_to` should be specified when ",
                "`final` is TRUE")
        } else {
            fname_suffix <- "_final"
        }
    }

    if(!is.null(dir_path)){
        result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    }

    reg_suffix <- ""
    reg_suffix <- ifelse(regularize, paste0("reg_top", topn, "_"), "")
    aggl_suffix <- paste0(use_aggl, "_")
    dist_suffix <- paste0(use_dist, "_")

    if(!is.null(dir_path)){
        fname <- file.path(result_dir_path, paste0(
            sname, "_dend_arch_list_",
            reg_suffix, dist_suffix, aggl_suffix, use_cutk,
            "clusters", fname_suffix
        ))
        ## extension .pdf is added in the .plot_dend_arch function downstream
    }else{
        fname <- NULL
    }

    clusts_reord <- seqArchR::collate_seqArchR_result(
        result = seqArchR_result, iter = iter, clust_method = "hc",
        aggl_method = use_aggl, dist_method = use_dist,
        regularize = regularize, topn = topn,
        collate = FALSE, return_order = TRUE,
        flag = list(debugFlag = FALSE, verboseFlag = TRUE)
    )

    ## these are in default seqArchR ordering
    clust_seqArchR_ord_list <- seqArchR::get_seqs_clust_list(
        seqs_clust_lab = seqArchR_result$seqsClustLabels[[iter]]
    )

    ## these are now ordered by the hc ordering
    clust_hc_ord_list <- lapply(clusts_reord$order, function(x) {
        clust_seqArchR_ord_list[[x]]
    })

    ordered_seqArchR_pl <- seqArchR::plot_arch_for_clusters(
        seqs = seqArchR::seqs_str(seqArchR_result),
        clust_list = clust_hc_ord_list, pos_lab = pos_lab,
        xt_freq = 5, set_titles = FALSE, show = FALSE, bits_yax = "auto"
    )

    ## reverse the ordering to fit the dendrogram order from visual top
    ordered_arch_pl2 <- lapply(rev(ordered_seqArchR_pl), function(pl) {
        pl <- pl +
            ggplot2::theme(
                axis.text = ggplot2::element_text(size = 0),
                axis.text.x = ggplot2::element_text(
                    angle = 0, vjust = 2, hjust = 0.5
                ),
                axis.text.y = ggplot2::element_text(vjust = 0.5),
                axis.title.y = ggplot2::element_text(size = 0),
                axis.ticks.length = ggplot2::unit(0.00, "cm"),
                plot.margin = ggplot2::unit(c(-0.1, 0, -0.4, -0.4), "cm")
            )
    })

    sam_foo <- cowplot::plot_grid(plotlist = ordered_arch_pl2, ncol = 1)

    ## This function plots the grid of dendrogram + seqlogos
    ## & also writes it to a PDF with provided fname
    if (!final) {
        sam_foo2 <- .plot_dend_arch(
            arch_plot = sam_foo, fname = fname,
            use_cutk = use_cutk,
            clusts = clusts_reord,
            use_ht = 60, plot_png = FALSE, lwd = 0.4,
            repel = TRUE, show_labels = TRUE,
            labels_track_height = 0.25, rect = TRUE,
            rect_fill = TRUE, color_labels_by_k = TRUE
        )
        return(list(curation_plot = sam_foo2,
                    clusters_list = clust_hc_ord_list))
    } else {
        temp_clusts <- stats::cutree(clusts_reord, k = use_cutk)
        names(temp_clusts) <- NULL

        temp_clusts <- .handle_cl_reassignment(
            need_change = need_change,
            change_to = change_to,
            temp_clusts = temp_clusts, ret_cl_assignments = TRUE
        )
        clust_list <- .handle_cl_reassignment(
            need_change = need_change,
            change_to = change_to,
            temp_clusts = temp_clusts, ret_cl_assignments = FALSE
        )


        seqs_clusters_as_list <- seqArchR::collate_clusters(
            to_clust = clust_list,
            orig_clust = seqArchR::get_seqs_clust_list(
                seqArchR_result$seqsClustLabels[[iter]]
            )
        )
        ##
        use_color <- scales::hue_pal()(length(unique(temp_clusts)))

        sam_foo2 <- .plot_dend_arch(
            arch_plot = sam_foo, fname = fname,
            use_ht = 60,
            use_cutk = use_cutk, # length(unique(temp_clusts)),
            clusts = clusts_reord, rect = TRUE, rect_fill = TRUE,
            label_cols = use_color[temp_clusts[clusts_reord$order]],
            k_colors = use_color,
            clust_assignment = clust_list,
            new_clusts = seqs_clusters_as_list,
            rawSeqs = seqArchR::seqs_str(seqArchR_result),
            palette = FALSE, plot_png = FALSE
        )
        return(list(curation_plot = sam_foo2,
                    clusters_list = seqs_clusters_as_list))
    }
}
## =============================================================================


## temp_clusts[c(21,2,29)] <- temp_clusts[c(19)]
## temp_clusts[c(8)] <- temp_clusts[c(11)]
##
## need_change <- list(c(21,2,29), c(8))
## change_to <- list(c(19), c(11))
##
## when no change is to be made, set both need_change and change_to to an empty
## list like so list()
##
## temp_clusts
##
## need_change is a list of clust IDs that will be reassigned
## change_to is a list of clusters they will be assigned to
## Both the lists have a one-to-one mapping, meaning that element 1 in list
## need_change is assigned element 1 in change_to list
##
## Use 0 for new, independent clusters
## ??Use negative numbers for independent but same clusters for two different
## elements
##
## If ret_cl_assignments is TRUE (i.e., return cluster assignments), changed
## temp_clusts is returned, Otherwise, a clust_list with clusterIDs is returned
##
.handle_cl_reassignment <- function(need_change, change_to,
                                    temp_clusts, ret_cl_assignments = TRUE) {
    if (!all(lengths(list(need_change, change_to)) == 0)) {
        stopifnot(length(need_change) == length(change_to))
        sec_list_lens <- lengths(change_to)
        if (!all(sec_list_lens == 1)) {
            stop("All elements of the `change_to` list should be length 1")
        }
        ###
        #
        ## Re-assign here
        ##
        ## Following scenarios can arise:
        ## -- An existing cluster becoming empty/null set due to re-assignments
        ##      -- members can be re-assigned to another existing cluster
        ##      -- members can be re-assigned to a new singleton cluster
        ##         containing only itself (use numeral 0 in change_to)
        ##      -- members can be re-assigned to a new cluster (non-existent
        ##         earlier), but this is not a singleton. This can be achieved
        ##         by using combine in the need_change list entry.
        ##
        ## -- Normal re-assignments: those that do not lead to empty clusters
        ##
        ##
        ## Two things can happen when new clusters get created:
        ## A. Some existing cluster can have all its elements re-assigned, then
        ## this cluster is empty
        ## B. Some elements can be re-assigned to completely independent
        ## clusters, i.e., new clusters. These are marked by numeral 0  or
        ## negative numbers in change_to.
        ## Instead of letting them create confusions, we first
        ## implement all non-zero and non-negative  re-assignments, and
        ## handle resultant null sets as below.
        ##
        ## Finally, implement type B re-assignments.
        ##

        ##
        ##
        ## Re-assignments here change_to is an existing cluster

        pos_idx <- which(unlist(change_to) > 0)
        if(length(pos_idx) > 0){
            for (i in seq_along(pos_idx)) {
                temp_clusts[need_change[[pos_idx[i]]]] <-
                    temp_clusts[change_to[[pos_idx[i]]]]
            }
        }
        # for (x in seq_along(need_change)) {
        #     if(change_to[[x]] > 0){
        #         temp_clusts[need_change[[x]]] <- temp_clusts[change_to[[x]]]
        #     }else{
        #
        #     }
        # }
        ## When reassigning, some clusters may become null sets. They should be
        ## omitted. The below code does that
        existing_clust <- sort(unique(temp_clusts))
        for (i in seq_along(existing_clust)) {
            idx <- which(temp_clusts == existing_clust[i])
            temp_clusts[idx] <- i
        }

        ##
        ## Type B re-assignments
        ##
        zero_idx <- which(unlist(change_to) == 0)
        if(length(zero_idx) > 0){
            ## Make further few clusters
            nCl <- length(unique(temp_clusts))
            for (i in seq_along(zero_idx)) {
                temp_clusts[need_change[[zero_idx[i]]]] <- nCl + i
            }
        }

        ## When reassigning, some clusters may become null sets. They should be
        ## omitted. The below code does that
        existing_clust <- sort(unique(temp_clusts))
        for (i in seq_along(existing_clust)) {
            idx <- which(temp_clusts == existing_clust[i])
            temp_clusts[idx] <- i
        }

        neg_idx <- which(unlist(change_to) < 0)
        if(length(neg_idx) > 0){
            ## Make further few clusters
            nCl <- length(unique(temp_clusts))
            for (i in seq_along(neg_idx)) {
                temp_clusts[need_change[[neg_idx[i]]]] <- nCl + i
            }
        }

        ## When reassigning, some clusters may become null sets. They should be
        ## omitted. The below code does that
        existing_clust <- sort(unique(temp_clusts))
        for (i in seq_along(existing_clust)) {
            idx <- which(temp_clusts == existing_clust[i])
            temp_clusts[idx] <- i
        }
    }
    if(!ret_cl_assignments){
        ## Also need to show alongside, how the final clusters' seqlogos look
        clust_list <- lapply(unique(temp_clusts), function(x) {
            which(temp_clusts == x)
        })
        return(clust_list)
    }else{
        return(temp_clusts)
    }
}
## =============================================================================


.plot_dend_arch <- function(arch_plot, fname, use_ht = 40, use_wd = 50,
                            use_cutk = 2, use_cuth = NULL, clusts,
                            clust_assignment = NULL,
                            new_clusts = NULL, rawSeqs = NULL, pos_lab = NULL,
                            plot_png = TRUE, ...) {
    dend_pl2 <- factoextra::fviz_dend(clusts,
        horiz = TRUE, main = "",
        k = use_cutk,
        h = use_cuth,
        ...
    )

    dend_pl2 <- dend_pl2 + ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(add = c(0.6, 0.2))
    ) +
        NULL

    if (is.null(new_clusts)) {
        sam_foo2 <- cowplot::plot_grid(dend_pl2, arch_plot,
            ncol = 2, align = "hv",
            rel_widths = c(0.40, 1), rel_heights = c(1, 1)
        )
        if(!is.null(fname)){
            ## pdf
            cowplot::ggsave2(paste0(fname, ".pdf"),
                plot = sam_foo2,
                width = use_wd, height = use_ht,
                units = "cm", dpi = 300
            )
            ## png
            if (plot_png) {
                cowplot::ggsave2(paste0(fname, ".png"),
                    plot = sam_foo2,
                    width = 40, height = use_ht,
                    units = "cm", dpi = 300
                )
            }
        }
    } else {

        ## new_clusts info is provided
        new_arch_pl <- seqArchR::plot_arch_for_clusters(
            seqs = rawSeqs,
            clust_list = new_clusts, pos_lab = pos_lab,
            xt_freq = 5, set_titles = FALSE, show = FALSE, bits_yax = "auto"
        )

        new_arch_pl2 <- lapply(seq_along(new_arch_pl), function(x) {
            pl <- new_arch_pl[[x]] +
                ggtitle(paste0(
                    "Obtained by combining: ",
                    paste(clust_assignment[[x]], collapse = ", ")
                )) +
                ggplot2::theme(
                    axis.text = ggplot2::element_text(size = 0),
                    axis.text.x = ggplot2::element_text(
                        angle = 0, vjust = 2, hjust = 0.5
                    ),
                    axis.text.y = ggplot2::element_text(vjust = 0.5),
                    axis.title.y = ggplot2::element_text(size = 0),
                    axis.ticks.length = ggplot2::unit(0.00, "cm"),
                    plot.title = element_text(margin = margin(1, 0, 0, 0)),
                    plot.margin = ggplot2::unit(c(-0.1, 0, -0.4, -0.4), "cm")
                )
        })
        sam_foo <- cowplot::plot_grid(plotlist = new_arch_pl2, ncol = 1)
        ##
        sam_foo2 <- cowplot::plot_grid(dend_pl2, arch_plot, sam_foo,
            ncol = 3,
            rel_widths = c(0.40, 1, 1),
            rel_heights = c(1, 1, 1),
            align = "hv"
        )
        if(!is.null(fname)){
            ## pdf
            cowplot::ggsave2(paste0(fname, ".pdf"),
                plot = sam_foo2,
                width = use_wd + 20, height = use_ht,
                units = "cm", dpi = 300
            )
            ## png
            if (plot_png) {
                cowplot::ggsave2(paste0(fname, ".png"),
                    plot = sam_foo2,
                    width = 60, height = use_ht,
                    units = "cm", dpi = 300
                )
            }
        }
    }

    return(sam_foo2)
}
## =============================================================================
