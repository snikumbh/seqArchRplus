#' seqArchRplus: A package for downstream analyses of promoter sequence
#' architectures/clusters identified by seqArchR
#'
#' This package facilitates various downstream analyses and visualizations
#' for clusters/architectures identified in promoter sequences. If a seqArchR
#' result object is available, and in addition, the CAGEr object or information
#' on the tag clusters, this package can be used for the following:
#'
#' - Order the sequence architectures by the interquantile widths (IQWs) of the
#' tag clusters they originate from. See `CAGEr` vignette for more information.
#'   \code{\link{order_clusters_iqw}}
#'
#' - Visualize distributions of IQW, TPM and conservation scores per cluster
#'   \code{\link{iqw_tpm_plots}}
#'
#' - Visualize the percentage annotations of promoter sequence per cluster
#'   \code{\link{per_cluster_annotations}}
#'
#' - Visualize heatmaps of occurrences of motifs (as words) in sequences
#'   \code{\link{plot_motif_heatmaps}} and \code{\link{plot_motif_heatmaps2}}
#'
#' - Visualize the above plots as (publication ready) combined panels viewable
#' for different samples as HTML reports
#'
#' - Following per cluster visualizations:
#'
#'     - sequence logos of architectures (including strand-separated ones)
#'       \code{\link{per_cluster_seqlogos}}
#'
#'     - distributions of promoters on different chromosomes/strands
#'       \code{\link{per_cluster_strand_dist}}
#'
#'     - GO term enrichments
#'       \code{\link{per_cluster_go_term_enrichments}}
#'
#' - Produce BED track files of seqArchR clusters for visulization in a genome
#' browser or IGV
#'   \code{\link{write_seqArchR_cluster_track_bed}}
#'
#' - Curate seqArchR clusters
#'   \code{\link{curate_clusters}}
#'
#' - (future) Generate HTML reports that help you navigate this wealth of
#' information with ease, and enable insights and hypotheses generation
#'
#' #' @section Functions for data preparation and manipulation:
#' \itemize{
#' \item \code{\link{prepare_data_from_FASTA}}
#' \item \code{\link{get_one_hot_encoded_seqs}}
#' }
#'
#'
#' @section Functions for visualizations:
#' \itemize{
#' \item \code{\link{plot_arch_for_clusters}}
#' \item \code{\link{plot_ggseqlogo_of_seqs}}
#' \item \code{\link{viz_bas_vec}}
#' \item \code{\link{viz_seqs_acgt_mat}}
#' \item \code{\link{viz_pwm}}
#' }
#'
#' @author Sarvesh Nikumbh
#' @docType package
#' @name seqArchRplus
NULL

