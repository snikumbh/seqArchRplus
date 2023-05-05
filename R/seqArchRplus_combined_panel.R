
#' @title Form a combined panel of three plots
#'
#' @description For a given sample, this function forms a combined panel of
#' three plots namely, the IQW-TPM boxplots, cluster sequence logos, and
#' annotations per cluster. All of these individual plots can be generated
#' using existing seqArchRplus functions
#'
#' @param iqw_tpm_pl The IQW-TPM plot generated using
#' \code{\link{iqw_tpm_plots}}
#'
#' @param seqlogos_pl The sequence logos oneplot obtained from
#' \code{\link{per_cluster_seqlogos}} by setting the argument
#' \code{`one_plot = TRUE`}
#'
#' @param annot_pl The annotations oneplot obtained from
#' \code{\link{per_cluster_annotations}} by setting the argument
#' \code{`one_plot = TRUE`}
#'
#' @details This functionality requires cowplot. The combined
#' plot panels are arranged horizontally.
#'
#' @return This function returns a \code{ggplot} plot object.
#'
#'
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 ylab theme element_blank
#'
#' @export
#'
#'
#' @author Sarvesh Nikumbh
form_combined_panel <- function(iqw_tpm_pl, seqlogos_pl, annot_pl){

    ## Remove y-axis text/tick labels
    annot_new_pl <- annot_pl + ggplot2::ylab(NULL) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank())

    panel_pl <- cowplot::plot_grid(iqw_tpm_pl, seqlogos_pl, annot_new_pl,
        ncol = 3, axis = "tblh",
        align = "hv")

    return(panel_pl)
}
