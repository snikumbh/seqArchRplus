

#' @title Visualize all sequences as an image
#'
#' @param sname The sample name
#' @param seqs The sequences
#' @param seqs_ord The order of sequences
#' @param pos_lab The position labels
#' @param xt_freq The frequency of xticks
#' @param yt_freq The frequency of yticks
#' @param f_height,f_width The height and width for the PNG image.
#' @param dir_path Specify the /path/to/directory to store results
#'
#' @return Nothing. PNG images are written to disk using the provided filename.
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
#' seqs_acgt_image(sname = "sample1",
#'                 seqs = raw_seqs,
#'                 seqs_ord = unlist(use_clusts),
#'                 pos_lab = -45:45,
#'                 dir_path = tempdir())
#'
#' @author Sarvesh Nikumbh
seqs_acgt_image <- function(sname, seqs, seqs_ord, pos_lab, xt_freq = 5,
                            yt_freq = 500, f_height = 1200,
                            f_width = 600, dir_path) {
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    fname <- file.path(
        result_dir_path,
        paste0(sname, "_ClusteringImage.png")
    )
    seqArchR::viz_seqs_acgt_mat(as.character(seqs[unlist(seqs_ord)]),
        pos_lab = pos_lab,
        xt_freq = 5,
        yt_freq = 500,
        f_height = 1200, f_width = 600,
        save_fname = fname,
        file_type = "png"
    )
}
## =============================================================================
