#' seqArchRplus: A package for downstream analyses of promoter sequence
#' architectures/clusters identified by seqArchR
#'
#' This package facilitates various downstream analyses and visualizations
#' for clusters/architectures identified in promoter sequences. If a seqArchR
#' result object is available, and in addition, the CAGEr object or information
#' on the tag clusters, per sequence annotations and GO term information, etc.,
#' this package can be used for the following:
#'
#' - Order the sequence architectures by the interquantile widths (IQWs) of the
#' tag clusters they originate from. See `CAGEr` vignette for more information.
#'
#' - Visualize distributions of IQW, TPM and conservation scores per cluster
#'
#' - Visualize the percentage annotations of promoter sequence per cluster
#'
#' - Visualize the above plots as (publication ready) combined panels viewable
#' for different samples as HTML reports
#'
#' - Following per cluster visualizations:
#'     - sequence logos of architectures
#'     - strand-separated sequence logos of architectures
#'     - distributions of promoters on different chromosomes/strands
#'     - GO term enrichments
#'
#' - Produce BED track files of seqArchR clusters for visulization in a genome
#' browser or IGV
#'
#' - Generate HTML reports that help you navigate this wealth of information
#' with ease, and enable insights and hypotheses generation
#'
#'
#'
#'
#' @docType package
#' @name seqArchRplus
NULL


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# if(getRversion() >= "2.15.1")  utils::globalVariables(c("seqArchRconfig"))
