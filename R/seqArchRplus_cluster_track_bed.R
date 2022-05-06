

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
#' file
#'
#' @param one_zip_all Logical. Specify TRUE when the facility to download BED
#' files for all clusters with one click is desired, otherwise FALSE. This is
#' relevant only when include_in_report is TRUE
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
#' @details Note: For providing downloadable links in the HTML report, the
#' complete BED files are encoded into base64 strings and embedded with the
#' HTML report itself. This considerably increases the size of the HTML file,
#' and can slow down loading of the HTML file in your browser.
#'
#' @return
#' When `include_in_report = FALSE`, the cluster information is written to disk
#' as BED track files that can be viewed in the genome browser or IGV.
#' Otherwise, HTML text is returned that can be included in the report as
#' downloadable links for each cluster BED file.
#' When `one_zip_all = TRUE`, a link to download all files zipped into one is
#' also provided to enable convenience.
#'
#' @export
#'
write_seqArchR_cluster_track_bed <- function(sname, clusts = NULL, info_df,
                                     one_zip_all = TRUE, org_name = NULL,
                                     dir_path = NULL, include_in_report = FALSE,
                                     strand_sep = FALSE){
    if(include_in_report){
        if(!requireNamespace("xfun", quietly = TRUE)){
            stop("Please install R package 'xfun' for ability to ",
                 "embed downloadable links for cluster BED files ",
                 "in your report.")
        }
    }

    ##
    cli::cli_alert_info("Preparing cluster-wise BED for ", sname)
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    bedFilesPath <- file.path(result_dir_path, "Cluster_BED_tracks")
    stopifnot(.check_and_create_dir(bedFilesPath))
    cli::cli_alert_info(paste0("Writing cluster BED track files at: ",
                               bedFilesPath))

    ## Write bed to disk -- each cluster in a separate bed file
    if(include_in_report){
        prefix_str <- "Individual"
        if(strand_sep) prefix_str <- "Strand-separated individual"
        cat(paste0("\n\n### ", prefix_str, " cluster track BED files [",
                   sname, "]\n\n"))
    }

    for(lo in seq_along(clusts)){
        fname_prefix <- paste0(org_name, "_TC_sample_", sname, "_cluster")
        track.name_prefix <- paste(sname, org_name, "clust", sep = "_")
        pos_fname <- paste0(fname_prefix, lo, "_plus_strand.bed")
        neg_fname <- paste0(fname_prefix, lo, "_minus_strand.bed")
        bedfname <- file.path(bedFilesPath, paste0(fname_prefix, lo, ".bed"))
        bedfname_m <- file.path(bedFilesPath, neg_fname)
        bedfname_p <- file.path(bedFilesPath, pos_fname)
        ##
        strands <- c("+", "-", "na")
        strands <- strands[ifelse(strand_sep, c(1,2), 3)]
        for(strn in strands){
            chosen_idx <- clusts[[lo]]
            bedFilename <- switch(strn,
                                  "+" = bedfname_p,
                                  "-" = bedfname_m,
                                  "na" = bedfname)
            ##
            if(strn != "na" && strand_sep){
                chosen_idx <- .get_strand_specific_indices(
                    df = info_df,
                    seq_ids_in_cl = clusts[[lo]],
                    strand_val = strn)
            }
            strand_str <- ""
            if(strand_sep) strand_str <- paste0("(", strn, "strand)")
            ## TEST: no records in a cluster, or on one or both strand
            if(length(chosen_idx) < 1){
                cli::cli_alert_warning(paste("Cluster ", lo, strand_str,
                                             ": No records"))
                if(include_in_report) .write_empty_string()
            }else{
                limit_df <- info_df[chosen_idx,]
                strand_track_str <- ""
                if(strn != "na"){
                    strand_track_str <- ifelse(strn == "+", "_plus_strand",
                                               "_negative_strand")
                }
                track.name <- paste0(track.name_prefix, lo, strand_track_str)
                .write_as_track_bed(limit_df, chosen_idx,
                                    track_name = track.name,
                                    bedFilename = bedFilename)
                if(include_in_report){
                    .create_dload_text(embedFile = TRUE,
                       use_path = bedFilename,
                       use_text = paste("Download coordinates for cluster",
                                        lo, strand_str, " as browser track"),
                       include_in_report = include_in_report)
                }
                ##
            }
        }
    }
    ## This facility to have all files zipped together for download makes sense
    ## , atleast for now, when it is to be included in the report
    if(one_zip_all && include_in_report){
        cat(paste0("### All cluster track files (one zipped folder of",
                   "all BED files)\n\n"))
        .create_dload_text(embedFile = FALSE, use_path = bedFilesPath,
                           use_text = "Download all track files (zip)",
                           include_in_report = include_in_report)
    }
} ## function Ends
##==============================================================================


# @title
#
# @param embedFile Logical. Set to TRUE if an on-disk file is to be ebedded. By
# setting FALSE, a directory can be embedded.
#
#
.create_dload_text <- function(embedFile = TRUE, use_path = NULL,
                               use_text = NULL, include_in_report = TRUE){
    if(include_in_report){
        dload_this <- xfun::embed_dir(use_path,
                            text = paste("Download all track files (zip)"))

        cat(paste0("\n<", dload_this$name,
                   " href=\"", dload_this$attribs$href,
                   "\" download=\"", dload_this$attribs$download,
                   "\">", dload_this$children[[1]][1],
                   "</a>\n"))
    }
}
##==============================================================================

.write_as_track_bed <- function(given_df, seq_ids, track_name, bedFilename,
                                write_tc_bound = TRUE){
    ##
    old_colnames <- colnames(given_df)
    if(write_tc_bound){
        df_gr <- GenomicRanges::makeGRangesFromDataFrame(given_df,
                                                     keep.extra.columns = TRUE,
                                                     seqnames.field = "chr")
        ## Because lower and upper boundaries of the quantiles can be anything
        ## as chosen by the user, we should match and find out what columns
        ## are these
        colIdx <- grep("q_", names(S4Vectors::mcols(df_gr)))
        ##
        new_gr <- GenomicRanges::narrow(df_gr,
                                start = as.integer(mcols(df_gr)[, colIdx[1]]),
                                end = as.integer(mcols(df_gr)[, colIdx[2]])
        )
        ## Check, they should be identical. If not, something is wrong!
        stopifnot(identical(df_gr$IQW, GenomicRanges::width(new_gr)))
        given_df <- as.data.frame(new_gr)
        ## make sure the first colnames is chr and not seqnames
        if(!identical(colnames(given_df), old_colnames)){
            colnames(given_df) <- c("chr", old_colnames[-1])
        }
    }

    track.description <- track_name
    write(paste('track name="', track_name,'" description="',
            track.description,'" visibility="pack"', ' itemRgb="On"', sep = ''),
          file = bedFilename, append = FALSE)
    use_score <- ifelse("phast" %in% colnames(given_df),
                        given_df$phast, given_df$domTPM)
    utils::write.table(data.frame(given_df$chr,
                              formatC(given_df$start, format = 'f', digits = 0),
                              formatC(given_df$end, format = 'f', digits = 0),
                              seq_ids,
                              score = use_score,
                              given_df$strand),
    file = bedFilename,
    append = TRUE, col.names = FALSE, row.names = FALSE,
    quote = FALSE, sep = '\t')
}
## =============================================================================
