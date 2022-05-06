

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
#' This is passed on to `seqArchR::collate_seqArchR_result`. See argument
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
#' If any element is to be put into a new cluster, use a numeral 0 in change_to
#'
#' @param final Logical, set to TRUE or FALSE
#'
#' @param dir_path The path to the directory where files will be written
#'
#' @details
#' This function helps the user work through the curation in at least three
#' steps.
#'
#' 1. This function uses hierarchical clustering to obtain a clustering result.
#' The resulting clustering is visualized by a dendrogram, color-coded cluster
#' assignments, and corresponding sequence logos. Using this visualization,
#' the user can identify/estimate the (nearly right) number of clusters to cut
#' the dendrogram. The first call uses K = 1.
#'
#' 2. Re-call the function with the identified value of K. Look at the
#' visualization to determine if it is good enough, i.e., it requires minor
#' re-assignments).
#'
#' 3. Identify cases of cluster assignments that you wish to re-assign to
#' different clusters. These can be noted as a list and supplied in a
#' subsequent call to the function.
#'
#' 2. In the final call to the function, set `final = TRUE`, supply the re-
#' assignments as two lists `need_change` and `change_to`.
#'
#' @return
#' When `final = TRUE`, the return value of the function is the clusters
#' with re-assignments executed. Otherwise, the clusters are ordered by
#' according to the hclust result/dendrogram.
#'
#' @importFrom stats cutree median
#'
#' @export
#'
curate_clusters <- function(sname, use_aggl = "ward.D", use_dist = "euclid",
                            seqArchR_result, iter, pos_lab = NULL,
                            regularize = TRUE, topn = 50, use_cutk = 2,
                            need_change = NULL, change_to = NULL,
                            final = FALSE, dir_path){
    cli::cli_h1(paste0("seqArchR result clusters curation"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    fname_suffix <- ""
    if(final){
        if(is.null(need_change) || is.null(change_to)){
            stop("Both `need_change` and `change_to` should be specified when
                 `final` is TRUE")
        }else{
            fname_suffix <- "_final"
        }
    }

    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)

    reg_suffix <- ""
    reg_suffix <- ifelse(regularize, paste0("reg_top", topn, "_"), "")
    aggl_suffix <- paste0(use_aggl, "_")
    dist_suffix <- paste0(use_dist, "_")
    fname <- file.path(result_dir_path, paste0(sname, "_dend_arch_list_",
                               reg_suffix, dist_suffix, aggl_suffix, use_cutk,
                               "clusters", fname_suffix))
    ## extension .pdf is added in the .plot_dend_arch function downstream

    clusts_reord <- seqArchR::collate_seqArchR_result(
        result = seqArchR_result, iter = iter, clust_method = 'hc',
        aggl_method = use_aggl, dist_method = use_dist,
        regularize = regularize, topn = topn,
        collate = FALSE, return_order = TRUE,
        flag = list(debugFlag = FALSE, verboseFlag = TRUE))

    ## these are in default seqArchR ordering
    clust_seqArchR_ord_list <- seqArchR::get_seqs_clust_list(
        seqs_clust_lab = seqArchR_result$seqsClustLabels[[iter]])

    ## these are now ordered by the hc ordering
    clust_hc_ord_list <- lapply(clusts_reord$order, function(x){
        clust_seqArchR_ord_list[[x]]
    })

    ordered_seqArchR_pl <- seqArchR::plot_arch_for_clusters(
        seqs = seqArchR::seqs_str(seqArchR_result),
        clust_list = clust_hc_ord_list, pos_lab = pos_lab,
        xt_freq = 5, set_titles = FALSE, show = FALSE, bits_yax = "auto")

    ## reverse the ordering to fit the dendrogram order from visual top
    ordered_arch_pl2 <- lapply(rev(ordered_seqArchR_pl), function(pl){
        pl <- pl +
            ggplot2::theme(axis.text = ggplot2::element_text(size = 0),
                       axis.text.x = ggplot2::element_text(
                           angle = 0, vjust = 2, hjust = 0.5),
                       axis.text.y = ggplot2::element_text(vjust = 0.5),
                       axis.title.y = ggplot2::element_text(size = 0),
                       axis.ticks.length = ggplot2::unit(0.00, "cm"),
                       plot.margin = ggplot2::unit(c(-0.1,0,-0.4,-0.4), "cm"))
    })

    sam_foo <- cowplot::plot_grid(plotlist = ordered_arch_pl2, ncol = 1)

    ## This function plots the grid of dendrogram + seqlogos
    ## & also writes it to a PDF with provided fname
    if(!final){
        sam_foo2 <- .plot_dend_arch(arch_plot = sam_foo, fname = fname,
                                    use_cutk = use_cutk,
                                    clusts = clusts_reord,
                                    use_ht = 60, plot_png = FALSE, lwd = 0.4,
                                    repel = TRUE, show_labels = TRUE,
                                    labels_track_height = 0.25, rect = TRUE,
                                    rect_fill = TRUE, color_labels_by_k = TRUE)

        return(ordered_arch_pl2)
    }else{
        temp_clusts <- stats::cutree(clusts_reord, k = use_cutk)
        names(temp_clusts) <- NULL
        clust_list <- .handle_cl_reassignment(need_change = need_change,
                                              change_to = change_to,
                                              temp_clusts = temp_clusts)

        seqs_clusters_as_list <- seqArchR::collate_clusters(
            to_clust = clust_list,
            orig_clust = seqArchR::get_seqs_clust_list(
                seqArchR_result$seqsClustLabels[[iter]])
        )
        ##
        use_color <- scales::hue_pal()(length(unique(temp_clusts)))
        sam_foo2 <- .plot_dend_arch(arch_plot = sam_foo, fname = fname,
                        use_ht = 60,
                        use_cutk = use_cutk,#length(unique(temp_clusts)),
                        clusts = clusts_reord, rect = TRUE, rect_fill = TRUE,
                        label_cols = use_color[temp_clusts[clusts_reord$order]],
                        k_colors = use_color,
                        clust_assignment = clust_list,
                        new_clusts = seqs_clusters_as_list,
                        rawSeqs = seqArchR::seqs_str(seqArchR_result),
                        palette = FALSE, plot_png = FALSE)

        return(seqs_clusters_as_list)
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
## Both the lists have a one-to-one mapping, meaning that element 1 is list
## need_change is assigned to element 1 in changes_to list
##
## Use 0 for new, independent clusters
## ??Use negative numbers for independent but same clusters for two different
## elements
##
.handle_cl_reassignment <- function(need_change, change_to,
                                    temp_clusts){
    if(!all(lengths(list(need_change, change_to)) == 0)){
        stopifnot(length(need_change) == length(change_to))
        sec_list_lens <- lengths(change_to)
        if(!all(sec_list_lens == 1)){
            stop("All elements of the `change_to` list should be length 1")
        }
        ###
        #
        ## Re-assign here
        ## Two things can happen:
        ## A. Some existing cluster can have all its elements re-assigned, then
        ## this cluster is empty
        ## B. Some elements can be re-assigned to completely independent
        ## clusters, i.e., new clusters. These are marked by numeral 0 in
        ## change_to. Instead of letting them create confusions, we first
        ## implement all non-zero re-assignments, and handle resultant
        ## null sets as below. Only then implement type B re-assignments.
        ##
        ## Type A re-assignments
        for( x in seq_along(need_change)){
            temp_clusts[ need_change[[x]] ] <- temp_clusts[ change_to[[x]] ]
        }
        ## When reassigning, some clusters may become null sets. They should be
        ## omitted. The below code does that
        existing_clust <- sort(unique(temp_clusts))
        for(i in seq_along(existing_clust)){
            idx <- which(temp_clusts == existing_clust[i])
            temp_clusts[idx] <- i
        }

        ##
        ## Type B re-assignments
        ##
        zero_idx <- which(unlist(change_to) == 0)
        ## Make further few clusters
        nCl <- length(unique(temp_clusts))
        for(i in seq_along(zero_idx)){
            temp_clusts[ need_change[[ zero_idx[i] ]] ] <- nCl + i
        }

    }
    ## Also need to show alongside, how the final clusters' seqlogos look
    clust_list <- lapply(unique(temp_clusts), function(x){
        which(temp_clusts == x)
    })
    clust_list
}
## =============================================================================

.plot_dend_arch <- function(arch_plot, fname, use_ht = 40, use_wd = 50,
                            use_cutk = 2, use_cuth = NULL, clusts,
                            clust_assignment = NULL,
                            new_clusts = NULL, rawSeqs = NULL, pos_lab = NULL,
                            plot_png = TRUE, ...
){

    dend_pl2 <- factoextra::fviz_dend(clusts, horiz = TRUE, main = "",
                                      k = use_cutk,
                                      h = use_cuth,
                                      ...
    )

    dend_pl2 <- dend_pl2 + ggplot2::scale_x_continuous(
        expand = ggplot2::expansion(add = c(0.6,0.2))) +
        NULL

    if(is.null(new_clusts)){
        sam_foo2 <- cowplot::plot_grid(dend_pl2, arch_plot,
                                   ncol = 2, align = "hv",
                                   rel_widths = c(0.40,1), rel_heights = c(1,1)
        )
        ## pdf
        cowplot::ggsave2(paste0(fname, ".pdf"), plot = sam_foo2,
                         width = use_wd, height = use_ht,
                         units = "cm", dpi = 300)
        ## png
        if(plot_png){
            cowplot::ggsave2(paste0(fname, ".png"), plot = sam_foo2,
                             width = 40, height = use_ht,
                             units = "cm", dpi = 300)
        }
    }else{

        ## new_clusts info is provided
        new_arch_pl <- seqArchR::plot_arch_for_clusters(
            seqs = rawSeqs,
            clust_list = new_clusts, pos_lab = pos_lab,
            xt_freq = 5, set_titles = FALSE, show = FALSE, bits_yax = "auto")

        new_arch_pl2 <- lapply(seq_along(new_arch_pl), function(x){
            pl <- new_arch_pl[[x]] +
                ggtitle(paste0("Obtained by combining: ",
                               paste(clust_assignment[[x]], collapse= ", "))) +
                ggplot2::theme(axis.text = ggplot2::element_text(size = 0),
                       axis.text.x = ggplot2::element_text(
                           angle = 0, vjust = 2, hjust = 0.5),
                       axis.text.y = ggplot2::element_text(vjust = 0.5),
                       axis.title.y = ggplot2::element_text(size = 0),
                       axis.ticks.length = ggplot2::unit(0.00, "cm"),
                       plot.title = element_text(margin=margin(1,0,0,0)),
                       plot.margin = ggplot2::unit(c(-0.1,0,-0.4,-0.4), "cm"))
        })
        sam_foo <- cowplot::plot_grid(plotlist = new_arch_pl2, ncol = 1)
        ##
        sam_foo2 <- cowplot::plot_grid(dend_pl2, arch_plot, sam_foo,
                                       ncol = 3,
                                       rel_widths = c(0.40,1,1),
                                       rel_heights = c(1,1,1),
                                       align = "hv")
        ## pdf
        cowplot::ggsave2(paste0(fname, ".pdf"), plot = sam_foo2,
                         width = use_wd + 20, height = use_ht,
                         units = "cm", dpi = 300)
        ## png
        if(plot_png){
            cowplot::ggsave2(paste0(fname, ".png"), plot = sam_foo2,
                             width = 60, height = use_ht,
                             units = "cm", dpi = 300)
        }
    }

    return(sam_foo2)

}
## =============================================================================
