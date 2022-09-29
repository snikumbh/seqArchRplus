

## 1. Files are written to disk
## 2. These on-disk files are embedded within the HTML and made available
## for download
## 3. If a BED file for any particular cluster is missing on disk, there is a
## good chance that the cluster was an empty set
##
#' @title Write seqArchR cluster information in BED files viewable as browser
#' tracks
#'
#' @description Writes the seqArchR clusters are BED tracks for viewing in IGV
#' or any genome browser
#'
#' @param sname Sample name
#'
#' @param clusts List of sequence ids in each cluster.
#'
#' @param info_df The data.frame holding information to be written to the BED
#' file. Expected columns are "chr", "start", "end", "names", "strand",
#' "domTPM", and "dominant_ctss". Out of these, only "dominant_ctss" column is
#' optional; when present this is visualised as thickStart and thickEnd.
#' The "domTPM" column is treated as scores for visualizing by default. One
#' can also visualize PhastCons score
#'
#' @param use_q_bound Logical. Write the lower and upper quantiles as tag
#' cluster boundaries. Default is TRUE
#'
#' @param use_as_names Specify the column name from info_df which you would like
#' to use as names for display with the track. By default, `use_names` is NULL,
#' and the sequence/tag cluster IDs are used as names.
#'
#' @param one_zip_all Logical. Specify TRUE when the facility to download BED
#' files for all clusters with one click is desired, otherwise FALSE. This is
#' relevant only when include_in_report is TRUE. Default is FALSE
#'
#' @param org_name Organism name
#'
#' @param dir_path The path where a special folder named `Cluster_BED_tracks`
#' (default) is created and all BED files are written inside it
#'
#' @param include_in_report Logical. Specify TRUE when this function is invoked
#' to write BED files to disk *and* provide downloadable links in the HTML
#' report. The corresponding chunk in the Rmarkdown report should set the
#' parameter `result='asis'`. By setting this to FALSE, BED files are written
#' to disk, but no downloadable links are provided. Note: This should be TRUE,
#' when `one_zip_all` argument is set to TRUE (see below). Requires the package
#' `xfun`.
#'
#' @param strand_sep Logical. Specify TRUE if records for each strand are to
#' be written in separate BED files.
#'
#' @details
#'
#' Note on links in HTML:
#' For providing downloadable links in the HTML report, the
#' complete BED files are encoded into base64 strings and embedded with the
#' HTML report itself. This considerably increases the size of the HTML file,
#' and can slow down loading of the HTML file in your browser.
#'
#' Note on BED files:
#' The output BED files have selected columns provided in the `info_df`.
#' These are "chr", "start", "end", "name", "score" (see more info below),
#'  "strand", "dominant_ctss".
#' By default, the sequence/tag cluster IDs are used as names.
#' If `use_as_names` is specified, information from that column in the
#' `info_df` is used as "name".
#'
#' If conservation score (e.g., PhastCons) is available, it is used as the
#' score, otherwise the TPM value of the dominant CTSS (domTPM) is used.
#' The final two columns (when dominantCTSS column is present), are the
#'  'thickStart' and 'thickEnd' values corresponding to the BED format.
#'  The 'thickEnd' column is the dominant_ctss position.
#'
#' Importantly, the lower and upper quantile boundaries are used as the start
#' and end coordinates of the cluster when `use_q_bound` is set to TRUE
#' (the default).
#'
#'
#' @return
#' When `include_in_report = FALSE`, the cluster information is written to disk
#' as BED track files that can be viewed in the genome browser or IGV.
#' Otherwise, a str object holding HTML text is returned that can be included
#' in the report as downloadable links for each cluster BED file (use `cat`).
#' When `one_zip_all = TRUE`, a link to download all files zipped into one is
#' also provided to enable convenience.
#'
#' @export
#'
#' @examples
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname,
#'          sep = "\t", header = TRUE)
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#'
#' # Write seqArchR clusters of promoters/CAGE tag clusters as BED track files
#' # All possible variations are enlisted here
#'
#' # Using quantiles information as the tag cluster boundaries
#' # via arg `use_q_bound = TRUE` and using custom names for each.
#' # Create a new/custom column, and specify the new column name for argument
#' # `use_as_names`. Notice that any custom names can be obtained by this
#' # approach.
#'
#' info_df$use_names <- paste(rownames(info_df), info_df$domTPM, sep = "_")
#'
#' write_seqArchR_cluster_track_bed(sname = "sample1",
#'                                  clusts = use_clusts,
#'                                  info_df = info_df,
#'                                  use_q_bound = FALSE,
#'                                  use_as_names = "use_names",
#'                                  dir_path = tempdir()
#'                                  )
#'
#' # Generating textual output that can be included in HTML reports.
#' # This requires package xfun.
#' \donttest{
#' cat_str <- write_seqArchR_cluster_track_bed(sname = "sample1",
#'                                  clusts = use_clusts,
#'                                  info_df = info_df,
#'                                  use_q_bound = FALSE,
#'                                  use_as_names = "use_names",
#'                                  dir_path = tempdir(),
#'                                  include_in_report = TRUE,
#'                                  one_zip_all = TRUE
#'                                  )
#' }
#'
#'
#' @author Sarvesh Nikumbh
write_seqArchR_cluster_track_bed <- function(sname, clusts = NULL, info_df,
                                                use_q_bound = TRUE,
                                                use_as_names = NULL,
                                                one_zip_all = FALSE,
                                                org_name = NULL,
                                                dir_path = NULL,
                                                include_in_report = FALSE,
                                                strand_sep = FALSE) {
    cli::cli_h1(paste0("Writing cluster track BED files"))
    cli::cli_h2(paste0("Sample: ", sname))

    cat_str <- ""
    if (include_in_report) {
        if (!requireNamespace("xfun", quietly = TRUE)) {
            stop(
                "Please install R package 'xfun' for ability to ",
                "embed downloadable links for cluster BED files ",
                "in your report."
            )
        }
    }

    ##
    cli::cli_alert_info(paste0("Preparing cluster-wise BED for ", sname))
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    bedFilesPath <- file.path(result_dir_path, "Cluster_BED_tracks")
    stopifnot(.check_and_create_dir(bedFilesPath))
    cli::cli_alert_info(paste0(
        "Writing cluster BED track files at: ",
        bedFilesPath
    ))

    ## Write bed to disk -- each cluster in a separate bed file
    if (include_in_report) {
        prefix_str <- "Individual"
        if (strand_sep) prefix_str <- "Strand-separated individual"
        cat_str <- paste0(cat_str, paste0(
            "\n\n### ", prefix_str, " cluster track BED files [",
            sname, "]\n\n"
        ))
    }

    if(use_q_bound)
        cli::cli_alert_info(paste0("Using quantiles as tag ",
            "cluster boundaries"))

    for (lo in seq_along(clusts)) {
        fname_prefix <- paste0(org_name, "_TC_sample_", sname, "_cluster")
        track.name_prefix <- paste(sname, org_name, "clust", sep = "_")
        pos_fname <- paste0(fname_prefix, lo, "_plus_strand.bed")
        neg_fname <- paste0(fname_prefix, lo, "_minus_strand.bed")
        bedfname <- file.path(bedFilesPath, paste0(fname_prefix, lo, ".bed"))
        bedfname_m <- file.path(bedFilesPath, neg_fname)
        bedfname_p <- file.path(bedFilesPath, pos_fname)
        ##
        strands <- c("+", "-", "na")
        strands <- strands[ifelse(strand_sep, c(1, 2), 3)]
        for (strn in strands) {
            chosen_idx <- clusts[[lo]]
            bedFilename <- switch(strn,
                "+" = bedfname_p,
                "-" = bedfname_m,
                "na" = bedfname
            )
            ##
            if (strn != "na" && strand_sep) {
                chosen_idx <- .get_strand_specific_indices(
                    df = info_df,
                    seq_ids_in_cl = clusts[[lo]],
                    strand_val = strn
                )
            }
            strand_str <- ""
            if (strand_sep) strand_str <- paste0("(", strn, "strand)")
            ## TEST: no records in a cluster, or on one or both strand
            if (length(chosen_idx) < 1) {
                cli::cli_alert_warning(paste(
                    "Cluster ", lo, strand_str,
                    ": No records"
                ))
                if (include_in_report) {
                    cat_str <- paste0(cat_str, .write_empty_string())
                }
            } else {
                limit_df <- info_df[chosen_idx, ]
                strand_track_str <- ""
                if (strn != "na") {
                    strand_track_str <- ifelse(strn == "+", "_plus_strand",
                        "_negative_strand"
                    )
                }
                track.name <- paste0(track.name_prefix, lo, strand_track_str)


                if(is.null(use_as_names)){
                    set_names <- chosen_idx
                }else{
                    ## make sure column exists
                    if(use_as_names %in% colnames(limit_df)){
                        set_names <- as.vector(limit_df[, use_as_names])
                    }else{
                        stop("Specified column name in `use_names` does not ",
                            "exist in `info_df`")
                    }
                }


                .write_as_track_bed(given_df = limit_df,
                    use_names = set_names,
                    track_name = track.name,
                    bedFilename = bedFilename,
                    use_q_bound = use_q_bound
                )

                if (include_in_report) {
                    .create_dload_text(
                        embedFile = TRUE,
                        use_path = bedFilename,
                        use_text = paste(
                            "Download coordinates for cluster",
                            lo, strand_str, " as browser track"
                        )
                    )
                }
                ##
            }
        }
    }
    ## This facility to have all files zipped together for download makes sense
    ## , atleast for now, when it is to be included in the report
    if (one_zip_all && include_in_report) {
        cat_str < paste0(cat_str, paste0(
            "### All cluster track files (one zipped folder of",
            "all BED files)\n\n"
        ))
        cat_str <- paste0(cat_str, .create_dload_text(
            embedFile = FALSE, use_path = bedFilesPath,
            use_text = "Download all track files (zip)"
        ))
        return(cat_str)
    }
} ## function Ends
## =============================================================================


# @title
#
# @param embedFile Logical. Set to TRUE if an on-disk file is to be ebedded. By
# setting FALSE, a directory can be embedded.
#
#
.create_dload_text <- function(embedFile = TRUE, use_path = NULL,
                                use_text = NULL) {
    if (!embedFile) {
        dload_this <- xfun::embed_dir(use_path, text = use_text)
    }else{
        dload_this <- xfun::embed_file(path = use_path, text = use_text)
    }

    cat_str <- paste0(
        "\n<", dload_this$name,
        " href=\"", dload_this$attribs$href,
        "\" download=\"", dload_this$attribs$download,
        "\">", dload_this$children[[1]][1],
        "</a>\n"
    )
    return(cat_str)
}
## =============================================================================

.write_as_track_bed <- function(given_df, use_names, track_name, bedFilename,
                                use_q_bound = TRUE) {
    ##
    ## given_df should have separate start and end columns
    ##
    old_colnames <- colnames(given_df)

    if (use_q_bound) {

        df_gr <- GenomicRanges::makeGRangesFromDataFrame(given_df,
            keep.extra.columns = TRUE
        )
        ## Because lower and upper boundaries of the quantiles can be anything
        ## as chosen by the user, we should match and find out what columns
        ## are these
        # colIdx <- grep("q_", names(S4Vectors::mcols(df_gr)))
        colIdx <- grep("q_", colnames(given_df))
        ## If use_q_bound is TRUE, colIdx should not be empty
        if(length(colIdx) == 0){
            cli::cli_alert_warning("`use_q_bound` is set to TRUE, but no
                columns with quantile information found")
        }else{
            ## We assume that the qLow appears first and then qUp

            ## This narrowing can fail when IQWs are larger than the range
            ## widths
            new_gr <- GenomicRanges::narrow(df_gr,
                # start = as.integer( GenomicRanges::mcols(df_gr)[, colIdx[1]] ),
                # end = as.integer( GenomicRanges::mcols(df_gr)[, colIdx[2]] )
                start = as.vector(given_df[, colIdx[1]]),
                end = as.vector(given_df[, colIdx[2]])
                )

            ## Check, they should be identical. If not, something is wrong!
            # stopifnot(identical(df_gr$IQW, GenomicRanges::width(new_gr)))
            stopifnot(identical(given_df$IQW, GenomicRanges::width(new_gr)))
            given_df <- as.data.frame(new_gr)

            ## make sure the first colnames is chr and not seqnames
            if (!identical(colnames(given_df), old_colnames)) {
                colnames(given_df) <- c("chr", old_colnames[-1])
            }

        }
    }

    track.description <- track_name
    write(paste('track name="', track_name, '" description="',
        track.description, '" visibility="pack"', ' itemRgb="On"',
        ' colorByStrand="255,0,0 0,0,255"',
        sep = ""
    ),
    file = bedFilename, append = FALSE
    )
    use_score <- ifelse("phast" %in% colnames(given_df),
        given_df$phast, given_df$domTPM
    )

    use_domctss <- ifelse("dominant_ctss" %in% colnames(given_df), TRUE, FALSE)
    ## We do not assume any dominant_ctss positions.
    ## The input DF can be appropriately modified before calling:
    ## add a new column by manuipulating the start and end coordinates

    if(use_domctss){
        utils::write.table(data.frame(given_df$chr,
            formatC(as.integer(given_df$start-1), format = "f", digits = 0),
            formatC(as.integer(given_df$end), format = "f", digits = 0),
            use_names, ## Name column
            score = use_score, ## Score column
            given_df$strand, ## Strand column
            as.integer(given_df$dominant_ctss)-1, ## thickStart
            as.integer(given_df$dominant_ctss) ## thickEnd
        ),
        file = bedFilename,
        append = TRUE, col.names = FALSE, row.names = FALSE,
        quote = FALSE, sep = "\t"
        )
    }else{
        utils::write.table(data.frame(given_df$chr,
            formatC(as.integer(given_df$start-1), format = "f", digits = 0),
            formatC(as.integer(given_df$end), format = "f", digits = 0),
            use_names, ## Name column
            score = use_score, ## Score column
            given_df$strand ## Strand column
        ),
            file = bedFilename,
            append = TRUE, col.names = FALSE, row.names = FALSE,
            quote = FALSE, sep = "\t"
        )
    }
}
## =============================================================================
