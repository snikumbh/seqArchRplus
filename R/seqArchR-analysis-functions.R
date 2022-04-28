## Functions for downstream analysis of seqArchR clusters
##
## Each of these functions performs the task for one specified sample
## 1. Tag clusters as bed files for viewing in Genome browser/IGV DONE
## 2. Per cluster annotation percentages DONE
## 3. Arranging clusters by IQW, TPM for comparison of architectures DONE
## 4. Seqlogos, usual and strand-separated.
## 5. How to visualize prevalence of clusters on different chromosomes
## 6. Heatmaps arranged by IQW
## 7. Manual curation of clusters

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
#' @param strand_sep Logical. Specify TRUE if records for each strand are to
#' be written in separate BED files.
#'
#' @details Note: For providing downloadable links in the HTML report, the
#' complete BED files are encoded into base64 strings and embedded with the
#' HTML report itself. This considerably increases the size of the HTML file,
#' and can slow down loading of the HTML file in your browser.
#'
#' @export
#'
write_seqArchR_cluster_track_bed <- function(sname, clusts = NULL, info_df,
                                one_zip_all = TRUE, org_name = NULL,
                                dir_path = NULL, include_in_report = TRUE,
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
        strands <- strands[ifelse(strand_sep, 1:2, 3)]
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


.get_strand_plot_title <- function(this_id, nclust, this_n = NULL,
                                    tot_n, strand_val = "+"){
    if(!is.null(strand_val)){
        if(is.null(this_n))
            stop("'this_n' cannot be NULL when strand_val is not NULL")
    }

    clust_names <- sort(as.character(seq(nclust)))

    if(is.null(strand_val)){
        title_str <-  paste0("(", this_id , "/", nclust,
                             ") Arch `", clust_names[this_id], "': ",
                             tot_n, " sequences")
    }else{
        if(strand_val == "+"){
            title_str <-  paste0("(", this_id , "/", nclust,
                                 ") Arch `", clust_names[this_id], "': ",
                                 this_n, " (+ strand) /", tot_n)
        }else if(strand_val == "-"){
            title_str <-  paste0("(", this_id , "/", nclust,
                                 ") Arch `", clust_names[this_id], "': ",
                                 this_n, " (- strand) /", tot_n)
        }
    }


    title_str
}
## =============================================================================

.get_strand_specific_indices <- function(df, seq_ids_in_cl, strand_val = "+"){
    return(seq_ids_in_cl[which(df$strand[seq_ids_in_cl] == strand_val)])
}
## =============================================================================

## Handle per sample result directory
.handle_per_sample_result_dir <- function(sname, dir_path){
    result_dir_path <- file.path(dir_path, paste0(sname, "_results"))
    stopifnot(.check_and_create_dir(result_dir_path))
    result_dir_path
}
## =============================================================================

## check_and_create_dir
.check_and_create_dir <- function(dir_path){
    creation_ok <- FALSE
    if(!dir.exists(dir_path)){
        cli::cli_alert_warning(paste0("Creating directory: ", dir_path))
        creation_ok <- dir.create(dir_path)
    }else{
        cli::cli_alert_warning(paste0("Directory exists: ", dir_path))
        creation_ok <- TRUE
    }
    ## TRUE: success in creation; FALSE: otherwise
    return(creation_ok)
}
## =============================================================================


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
          file = bedFilename, append = F)
    use_score <- ifelse("phast" %in% colnames(given_df),
                        given_df$phast, given_df$domTPM)
    utils::write.table(data.frame(given_df$chr,
                           formatC(given_df$start, format = 'f', digits = 0),
                           formatC(given_df$end, format = 'f', digits = 0),
                           seq_ids,
                           score = use_score,
                           given_df$strand
    ),
    file = bedFilename,
    append = T, col.names = F, row.names = F,
    quote = F, sep = '\t')
}
## =============================================================================


.write_empty_string <- function(){
    cat(paste0("<a href= >Empty",   "</a>\n"))
}
##==============================================================================


#' @title IQW, TPM plots
#'
#' @param sname Sample name
#'
#' @param dir_path The path to the directory where files will be written.
#'
#' @param info_df A DataFrame object holding information on the tag clusters.
#' Expected columns (names) are 'chr', 'start', 'end', 'IQW', 'domTPM', and
#' 'strand'.
#'
#' @param clusts List of sequence ids in each cluster.
#'
#' @param iqw Logical. Specify TRUE if boxplots of interquantile widths (IQW)
#' of the tag clusters corresponding to the promoters in each cluster
#' are to be plotted.
#'
#' @param tpm Logical. Specify TRUE if boxplots of TPM values for all clusters
#' are to be plotted.
#'
#' @param cons Logical. Specify TRUE if boxplots of conservation scores
#' (PhastCons scores) for all clusters. If this is TRUE, an additional column
#' named 'cons' is expected in the `info_df`.
#'
#' @param txt_size Specify text size to be used in the plots.
#'
#' @param wTitle Logical. If TRUE, the returned plot will contain a default
#'  title, which is the same as the filename. See details.
#'
#' @param use_suffix,use_prefix Character. Specify any suffix and/or prefix
#' you wish to add to the filename.
#'
#' @details The plots are written to a file named
#' "Sample_<sample_name>_IQW_TPM_Cons_plot.pdf" if all of `iqw`, `tpm`,
#' and `cons` are set to TRUE. This is also set as the plot title if `wTitle`
#' is set to TRUE.
#'
#' All plots are arranged by the IQWs (smallest on top, largest at the bottom),
#' even iff `iqw` is set to FALSE.
#'
#' @return The plot(s) as a ggplot2 object. The order of the plots is IQW,
#' followed by TPM values, followed by conservation scores.
#'
#' @importFrom ggplot2 aes theme theme_void element_text element_blank
#' ggplot scale_y_continuous scale_y_discrete geom_boxplot dup_axis unit
#' ggtitle margin position_fill guide_legend geom_bar expansion xlab ylab
#' scale_x_log10 theme_bw
#' @importFrom cli cli_h1 cli_h1 cli_alert_warning cli_alert_info
#' @importFrom stats reorder
#' @importFrom forcats fct_reorder
#'
#' @export
#'
iqw_tpm_plots <- function(sname, dir_path, info_df, clusts, iqw = TRUE,
                          tpm = TRUE, cons = TRUE, txt_size = 12,
                          use_suffix = NULL, use_prefix = "C",
                          wTitle = TRUE){
    cli::cli_h1(paste0("IQW-ordered boxplots"))
    cli::cli_h2(paste0("Sample: ", sname))

    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)

    info_df$clust_ID <- .get_clust_id_column(info_df, clusts)
    ##
    if(!any(iqw, tpm, cons)){
        cli::cli_alert_warning(paste0("None of `iqw`, `tpm` or `cons` are ",
                                      "set to TRUE. Setting `iqw` to `TRUE`"))
        iqw <- TRUE
    }

    ##
    iqw_str <- ifelse(iqw, "_IQW", "")
    tpm_str <- ifelse(tpm, "_TPM", "")
    cons_str <- ifelse(cons, "_Cons", "")

    cli::cli_h2(paste("Plotting:", iqw_str, tpm_str, cons_str, "plots"))

    title_str <- paste0("Sample", sname, iqw_str, tpm_str, cons_str, "_plot")
    fname <- file.path(result_dir_path, paste0(title_str, ".pdf"))
    ##
    plots <- c(iqw, tpm, cons)
    which_plots <- which(plots)
    ## Can store as many plots as asked for (by setting TRUE)
    ## begin by all entries in plot_list to NULL, will later filter out NULLs
    plot_list <- vector("list", length(which_plots))
    set_yax_txt <- FALSE

    ##
    plot_list <- lapply(which_plots, function(x){
        set_iqw <- ifelse(x == 1, TRUE, FALSE)
        set_tpm <- ifelse(x == 2, TRUE, FALSE)
        set_cons <- ifelse(x == 3, TRUE, FALSE)

        ## When we enter the loop for the first time, that is going to be for
        ## the first plot. This is when we need the y-axis text. From the
        ## second time onwards, there is a plot on the left to what is being
        ## drawn now. Therefore, no need for y-axis text.
        set_yax_txt <- ifelse(which(which_plots == x) == 1, TRUE, FALSE)
        ##
        temp_pl <- .get_iqw_ord_plot(iqw = set_iqw, tpm = set_tpm,
                            phast = set_cons,
                            y_axis_text = set_yax_txt,
                            use_notch = FALSE, txt_size = txt_size,
                            info_df = info_df, seqs_clust = clusts,
                            use_suffix = use_suffix)
        temp_pl
    })

    ## Keep only non-null entries
    plot_list <- lapply(plot_list, function(x) if(!is.null(x)) x)

    ## IQW ordered plots
    if(length(plot_list) > 1){
        rel_width_vals <- c(1, rep(0.85, times = length(plot_list)-1))
    }else{
        rel_width_vals <- 1
    }

    comb_pl <- cowplot::plot_grid(plotlist = plot_list, nrow = 1,
                                  ncol = length(plot_list), align = "h",
                                  rel_widths = rel_width_vals)
    ##
    pl_w_title <- comb_pl + ggplot2::ggtitle(label = title_str)
    ##
    if(!is.null(fname))
        cowplot::save_plot(filename = fname, plot = pl_w_title,
                       base_height = 14, base_width = 5*length(plot_list))
    ##
    if(wTitle){
        return(pl_w_title)
    }else{
        return(comb_pl)
    }
}
##==============================================================================

.get_clust_id_column <- function(info_df, clusts){
    ## Add new column noting the cluster IDs from seqArchR result
    clust_lab <- rep("0", length(unlist(clusts)))
    clust_names <- sort(as.character(1:length(clusts)))

    for(i in seq_along(clusts)){
        clust_lab[ clusts[[i]] ] <- clust_names[i]
    }
    clust_lab
}
##==============================================================================

make_cluster_labels <- function(clust, use_prefix, use_suffix){
    clust_lens <- lengths(clust)
    clust_labels <- paste(paste0(" (", clust_lens, ") "),
                          use_prefix, 1:length(clust), use_suffix,
                          sep="")
    clust_labels
}
##==============================================================================


.get_iqw_ord_plot <- function(iqw = FALSE, tpm = FALSE, phast = FALSE,
                             use_notch = FALSE, y_axis_text = TRUE,
                             txt_size = 12, info_df, seqs_clust,
                             use_suffix = "X", use_prefix = "C"){

    # clust_lens <- unlist(lapply(seqs_clust, length))
    # clust_labels <- paste(use_prefix, 1:length(seqs_clust), use_suffix,
    #                       paste0(" (", clust_lens, ")"), sep="")
    clust_labels <- make_cluster_labels(clust = seqs_clust,
                                        use_prefix = use_prefix,
                                        use_suffix = use_suffix)
    clr <- RColorBrewer::brewer.pal(3, "Dark2")
    ##
    if(iqw){
        pl <- ggplot(info_df,
                     aes(y=forcats::fct_reorder(clust_ID, IQW,
                                       .fun = median,
                                       .desc = TRUE),
                         x=IQW)) +
            ggplot2::geom_boxplot(outlier.size = 1, width = 0.5,
                                  notch = use_notch, color = "black",
                                  fill = clr[1]) +
            # ggplot2::scale_y_discrete(expand = expansion(add = c(0, 0))) +
            # scale_x_log10(#trans = 'log10',
            # breaks = scales::trans_format("log10", function(x) 10^x),
            # labels = scales::trans_format("log10", scales::math_format(10^.x))
            # ) +
            # coord_cartesian(xlim = c(0,25)) +
            xlab("IQW")
    }
    ##
    if(tpm){
        pl <- ggplot(info_df,
                     aes(y=forcats::fct_reorder(clust_ID, IQW,
                                       .fun = median,
                                       .desc = TRUE),
                         x=domTPM)) +
            geom_boxplot(outlier.size = 1, width = 0.5, notch = use_notch,
                         color = "black", fill = clr[2]) +
            # ggplot2::scale_y_discrete(expand = expansion(add = c(0.75, 0))) +
            # theme_linedraw() +
            # scale_x_log10(#trans = 'log10',
            # breaks = scales::trans_format("log10", function(x) 10^x),
            # labels = scales::trans_format("log10", scales::math_format(10^.x))
            #     ) +
            xlab("TPM") +
            ylab("Clusters")
        # +
        #     scale_y_discrete(labels = rev(clust_labels))
    }
    if(phast){
        pl <- ggplot(info_df,
                     aes(y=forcats::fct_reorder(clust_ID, IQW,
                                       .fun = median,
                                       .desc = TRUE),
                         x=phast)) +
            geom_boxplot(outlier.size = 1, width = 0.5, notch = use_notch,
                         color = "black", fill = clr[3]) +
            xlab("PhastCons score") +
            ylab("Clusters")
    }
    ##
    if(tpm || iqw){
        pl <- pl + scale_x_log10() +
            ggplot2::annotation_logticks(sides="b",
                                         short = unit(0.1, "cm"),
                                         mid = unit(0.2, "cm"),
                                         long = unit(0.3, "cm"))
    }
    pl <- pl  +
        ylab("Clusters") +
        scale_y_discrete(labels = rev(clust_labels),
                         expand = expansion(add = c(0.55, 0.55))) +
        theme_bw() +
        theme(#panel.grid = element_line(color = "grey90"),
              axis.title = element_text(colour = "black", size = txt_size),
              axis.text.x = element_text(colour = "black", angle = 0,
                                    size = txt_size, vjust = 0.5, hjust = 0.5),
              axis.text.y = element_text(colour = "black",
                                         size = txt_size, hjust = 1)
        )
    ##

    if(!y_axis_text){
        pl <- pl  +
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank())
    }

    return(pl)
    ####
}
##==============================================================================


#' @title per_cluster_annotations
#'
#' @param sname Sample name
#' @param clusts List of sequence ids in each cluster.
#' @param tc_gr Tag clusters as \code{\link[GenomicRanges]{GRanges}}. If
#' `cage_obj` is not provided (is NULL), this argument is required.
#' It will be ignored only if cage_obj is provided.
#' @param cage_obj A CAGEexp object obtained from the CAGEr package, if and
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
#' you wish to add to the filename.
#'
#' @export
#'
per_cluster_annotations <- function(sname, clusts, tc_gr,
                                    cage_obj = NULL,
                                    qLow = 0.1,
                                    qUp = 0.9,
                                    txdb_obj = NULL,
                                    tss_region = NULL,
                                    orgdb_obj = NULL,
                                    one_plot = TRUE,
                                    dir_path = NULL,
                                    txt_size = 12,
                                    use_suffix = NULL, use_prefix = "C"){
    message("Text size is: ", txt_size)
    cli::cli_h1(paste0("All clusters' genomic annotations"))
    cli::cli_h2(paste0("Sample: ", sname))
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    fname <- file.path(result_dir_path, paste0("Clusterwise_annotations.pdf"))
    ##
    # clust_labs <- unlist(lapply(seq_along(clusts),
    #                       function(x){
    #                           paste0(x, "_(n=", length(clusts[[x]]), ")")
    #                       }))
    clust_labels <- make_cluster_labels(clust = clusts, use_prefix, use_suffix)
    if(is.null(tc_gr)){
        if(is.null(cage_obj)){
            stop("If tc_gr is NULL, maybe you forgot to supply the `cage_obj`.",
                 "Please also specify `qLow` and `qUp` with the `cage_obj`")
        }else{
            any_null <-any(unlist(lapply(list(qLow, qUp), is.null)))
            if(any_null) stop("Please specify both `qLow` and `qUp`.")
            message("Using qLow = ", qLow, " and qUp = ", qUp)
            tc_gr <- CAGEr::tagClustersGR(cage_obj, sample = sname,
                                returnInterquantileWidth = TRUE,
                                qLow = qLow, qUp = qUp)


        }
    }
    stopifnot(!is.null(tc_gr))

    clustwise_anno <- lapply(clusts, function(x){
        foo_anno <- ChIPseeker::annotatePeak(tc_gr[x,],
                                             tssRegion=tss_region,
                                             TxDb = txdb_obj,
                                             annoDb = orgdb_obj)
        foo_anno
    })
    names(clustwise_anno) <- seq(1, length(clusts))

    colrs <- RColorBrewer::brewer.pal(n = 9, name = "Paired")
    ## Note: Generally, speed is not an issue when the number of clusters is
    ## around 40-50 or even up to 100.
    ## So we can use the base::rbind instead of the dplyr::bind_rows which is
    ## generally faster.
    ##

    ## without dplyr?
    df_list <- lapply(clustwise_anno, function(x) x@annoStat)
    # print(df_list)
    sam <- do.call("rbind", df_list)
    per_df_nrows <- unlist(lapply(df_list, nrow))
    ## Add a column clust that we need downstream
    sam$clust <- rep(seq_along(df_list), times = per_df_nrows)

    sam <- sam[, c(3,1,2)]
    rownames(sam) <- seq_len(sum(per_df_nrows))
    ##

    sam$Feature <- factor(sam$Feature,
                  levels = c("Promoter", "5' UTR", "1st Exon", "Other Exon",
                             "1st Intron", "Other Intron", "Downstream (<=300)",
                             "3' UTR", "Distal Intergenic"))
    ##
    names(colrs) <- levels(sam$Feature)
    ##
    if(one_plot){
    ## Solution using custom plotting
    clustwise_annobar <- get_prop_anno_oneplot(sam, txt_size = txt_size,
                                               colrs = colrs)

    ## Put in the cluster names when printing independently
    ind_print_plot <- clustwise_annobar +
        ggplot2::scale_y_discrete(labels = rev(clust_labels),
                                  expand = expansion(add = c(0.5, 0.5))) +
        ggplot2::theme(axis.text.y = element_text(size = txt_size,
                                                  hjust = 0.85))
    ##
    ggplot2::ggsave(filename = fname, plot = ind_print_plot,
                    device = "pdf", width = 15, height = 20, units = "in",
                    dpi = 300)

    return(clustwise_annobar)
    }else{
        ## return individual plots as a list
        sam_split <- split(sam,  f = factor(sam$clust,
                                           levels = seq_along(clusts)))
        annobar_list <- lapply(sam_split, function(anno_df){
            ##
            anno_df$Feature <- factor(anno_df$Feature,
                                      levels = levels(sam$Feature))
            anno_pl <- get_prop_anno_listplot(anno_df, txt_size = txt_size,
                                              colrs = colrs)
        })
        return(annobar_list)
    }

}
##==============================================================================

get_fixed_anno_ord <- function(){

    anno_terms_ord <- c("Promoter", "5' UTR", "3' UTR", "1st Exon",
                        "Other Exon", "1st Intron", "Intron",
                        "Downstream",
                        "Distal Intergenic")

    anno_terms_ord
}
##==============================================================================


get_named_colors <- function(anno_terms_ord, palname = "Set1"){
    use_colors <- get_ncolors(n=length(anno_terms_ord), palname=palname)
    if(palname=="Set1") use_colors <- use_colors[c(2:length(use_colors),1)]
    names(use_colors) <- anno_terms_ord
    use_colors
}
##==============================================================================

## anno_df holds information on the annoStats from the ChIPSeeker annotatePeak
## return object
## Expected columns are 'Frequency', 'Feature', 'clust'.
##
##
get_prop_anno_oneplot <- function(anno_df, txt_size, colrs){
    clustwise_annobar <- ggplot(anno_df,
                                aes(y = forcats::fct_rev(
                                    stats::reorder(clust, as.numeric(clust))),
                                    x = Frequency, fill = Feature)) +
        geom_bar(stat = "identity", width = 0.6,
                 position = position_fill(reverse=TRUE)) +
        ggplot2::scale_fill_manual(name = "", values = colrs) +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.text.y = element_blank(),
                       axis.ticks.length.y.left = unit(0.1, units = "cm"),
                       axis.ticks.length.x.bottom = unit(0.1, units = "cm"),
                       axis.text.x.bottom = element_text(size = txt_size),
                       axis.title.x.bottom = element_text(size = txt_size),
                       legend.position = "bottom",
                       legend.text = element_text(size = txt_size)
        ) +
        ggplot2::scale_y_discrete(expand = expansion(add = c(0.75, 0))) +
        ggplot2::scale_x_continuous(expand = expansion(mult = 0.01)) +
        ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE
        )) +
        ggplot2::ylab(NULL) +
        ggplot2::xlab("Percentage (%)")
    ##
    clustwise_annobar
}
##==============================================================================


get_prop_anno_listplot <- function(anno_df, txt_size, colrs){
    # message("text size is: ", txt_size)
    pl <- ggplot(anno_df, aes(y = clust, x = Frequency, fill = Feature)) +
        ggplot2::theme_bw() +
        ggplot2::geom_col(aes(fill = Feature), width=0.98,
                          position = position_fill(reverse=TRUE)) +
        ggplot2::scale_fill_manual(values= colrs) +
        ggplot2::theme(
            axis.title.y = element_text(size = txt_size),
            axis.title = element_text(size = txt_size),
            axis.text = element_text(size = txt_size),
            axis.text.y = element_text(hjust = 1.0, size = txt_size),
            # axis.text.x = element_blank(),
            # axis.title.x = element_blank(),
            # axis.ticks.x = element_blank(),
            panel.grid = element_blank(),
            legend.text = element_text(size = txt_size),
            legend.title = element_text(size = txt_size),
            legend.position = "bottom",
            legend.direction = "horizontal") +
        ggplot2::guides(fill = guide_legend(ncol = 7, byrow = FALSE
        )) +
        ggplot2::labs(x = "Percentage (%)", y = "Cluster") +
        NULL
    return(pl)
}
##==============================================================================

#' @title Handle writing of tag clusters to disk as BED files.
#'
#' @description Handle writing of tag clusters to disk as BED files. It can also
#' return corresponding promoter sequences as a DNAStringSet object or write
#' them to disk as FASTA files.
#'
#' @param sname The sample name.
#'
#' @param cage_obj The CAGEexp object from CAGEr.
#'
#' @param fl_size_up,fl_size_down Numeric. The size of the flanks in the
#' upstream and downstream directions.
#'
#' @param bsgenome The BSgenome file that will be used to obtain sequences of
#' the organism from.
#'
#' @param dir_path The path to the directory where files will be written. By
#' default, all BED files are written within a subdirectory named "BED", and
#' all FASTA files are written within a subdirectory named "FASTA", both
#' created at the `dir_path` location.
#'
#' @param fname_prefix,fname_suffix Specify any prefix or suffix string to be
#' used in the filename. This can be the organism name etc. Specify without
#' any separator. By default, an underscore is used as a separator in the
#' filename.
#'
#' @param write_to_disk Logical. Specify TRUE to write files to disk.
#' More specifically, BED files are written to disk only when this is set to
#' TRUE. For promoter sequences, FASTA files are written to disk if this arg is
#' set to TRUE, otherwise not. and a DNAStringSet object is returned if
#' `ret_seqs` is set to TRUE.
#'
#' @param ret_seqs Logical. Specify TRUE if promoter sequences are to be
#' returned as a DNAStringSet object.
#'
#' @return If `ret_seqs` is TRUE, a DNAStringSet object is returned.
#' Depending on `write_to_disk`, files are written to disk at the specified
#' location.
#'
#' @details
#' You can use the fname_prefix and fname_suffix arguments to specify strings
#' to be used as prefix and suffix for the files. For example, the organism
#' name can be used as a prefix for the filename. Similarly, for suffix.
#'
#' @importFrom GenomicRanges promoters trim
#' @importFrom IRanges IRanges
#' @importFrom CAGEr tagClustersGR
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqinfo seqinfo<-
#' @importFrom S4Vectors mcols
#'
#' @export
## Writes TCs as BED files via Granges and also as fasta files with a fixed
## flank size of +/- 500 bp
##
##
handle_tc_from_cage <- function(sname, cage_obj,
                                fl_size_up = 500, fl_size_down = 500,
                                bsgenome, dir_path = NULL,
                                fname_prefix = NULL, fname_suffix = NULL,
                                write_to_disk = TRUE, ret_seqs = TRUE){
    cli::cli_h1(paste0("Handling TCs from CAGE (BED and FASTA)"))
    cli::cli_h2(paste0("Sample: ", sname))

    bed_fname <- file.path(dir_path, "BED",
                            paste0(fname_prefix, "_TC_sample_",
                            sname, "_", fname_suffix, ".bed"))
    fasta_fname <- file.path(dir_path, "FASTA",
                            paste0(fname_prefix, "_sample_",
                            sname, "_promoters_",
                            fname_suffix, ".fa"))

    ##
    this_gr <- CAGEr::tagClustersGR(cage_obj, sample = sname,
                                    returnInterquantileWidth = TRUE,
                                    qLow = 0.1, qUp = 0.9)
    ## At present CAGEr tagClusters loose this information, so add them
    seqlevels(this_gr) <- GenomeInfoDb::seqlevels(CAGEr::CTSStagCountSE(cage_obj))
    seqinfo(this_gr)   <- GenomeInfoDb::seqinfo(CAGEr::CTSStagCountSE(cage_obj))

    if(any("score" == names(S4Vectors::mcols(this_gr)))){
        this_gr$tpm <- as.numeric(this_gr$score)
    }else{
        stop("Expecting an mcol named 'score', did not find.")
    }
    this_gr_df <- as.data.frame(this_gr)

    ## Storing dominant CTSS information into GRanges object...
    this_gr_domCtss <- GenomicRanges::GRanges(
                            seqnames = this_gr_df[, "seqnames"],
                            ranges = IRanges::IRanges(
                               start = this_gr_df[, "dominant_ctss"],
                               width = 1),
                            strand = this_gr_df[,"strand"])



    this_prom <- GenomicRanges::promoters(this_gr_domCtss, upstream = fl_size_up,
                    downstream = fl_size_down+1)
    this_prom <- GenomicRanges::trim(this_prom, use.names = TRUE)
    cli::cli_alert_success("Making promoters")

    ## Exclude indices where the sequences have been trimmed
    omit_ids <- which(Biostrings::width(this_prom) != ((fl_size_up  + fl_size_down)+1))

    message("Omitted IDs: ", length(omit_ids))

    if(length(omit_ids) > 0){
        this_gr <- this_gr[-omit_ids, ]
        this_prom <- this_prom[-omit_ids]
        this_gr_domCtss <- this_gr_domCtss[-omit_ids]
    }

    if(write_to_disk) write_tc_bed(this_gr, bed_fname)
    cli::cli_alert_success("BED file written")
    prom_seqs <- write_tc_fasta(this_prom, bsgenome = bsgenome,
                                fasta_fname = fasta_fname,
                                ret_val = ret_seqs,
                                also_write = write_to_disk)
    cli::cli_alert_success("FASTA file written")
    if(ret_seqs) return(prom_seqs)
}
##==============================================================================

##
#' @title Write Tag Clusters to a BED file
#'
#' @param gr The GenomicRanges object
#' @param bed_fname The BED filename that will be written to disk (with
#' full path)
#'
#' @details Writes the provided `gr` to the file named `bed_fname`
#'
#' @importFrom utils write.table
#'
#'
write_tc_bed <- function(gr, bed_fname){
    ##
    ## Write regions to BED file
    message("Writing tag cluster coordinates in BED file at: ", bed_fname)
    ##
    gr_df <- as.data.frame(gr)
    utils::write.table(
        gr_df,
        # [ , c("seqnames", "start", "end", "interquantile_width",
        #         "tpm.dominant_ctss", "strand")],
        file = bed_fname,
        sep = "\t",
        col.names = TRUE,
        quote = FALSE,
        row.names = FALSE)
}
##==============================================================================

# write_cage_tc_to_disk <- function(cage_obj){
#     message("Saving CAGE TC for ", sname)
#     saveRDS(myCAGEobject_thisSample, file = cage_tc_granges_fname)
# }

write_tc_fasta <- function(prom, fl_size_up = 500, fl_size_down = 500,
                        bsgenome, fasta_fname, ret_val = TRUE,
                        also_write = TRUE){

    cli::cli_alert_info("Fetching sequences...")

    fasta_names <- paste0("domCTSS=", seqnames(prom), ":",
                          start(prom), ";strand=",
                          strand(prom), ";",
                          "up=", fl_size_up, ";",
                          "down=", fl_size_down)

    prom_seqs <- BSgenome::getSeq(bsgenome,
                                    names= seqnames(prom),
                                    start = start(prom),
                                    end = end(prom),
                                    strand = strand(prom))
    names(prom_seqs) <- fasta_names
    ##
    if(also_write){
        cli::cli_alert_info(paste0("Writing FASTA file at: ", fasta_fname))
        Biostrings::writeXStringSet(prom_seqs,
                                filepath = fasta_fname,
                                format = "FASTA")
    }
    if(ret_val) return(prom_seqs)

}
##==============================================================================

#' @title Plot per cluster sequence logos
#'
#' @param sname The sample name
#' @param seqs The raw sequences as a DNAStringSet. These are also available
#' as part of the seqArchR result object.
#' @param clusts List of sequence ids in each cluster.
#' @param pos_lab The position labels
#' @param bits_yax The yaxis limits. See argument in seqArchR
#' @param strand_sep Logical. Whether sequences are to be separated by strand.
#' @param one_plot Logical. Whether all sequence logos should be combined
#' into one grid (with ncol = 1)?
#'
#' @param info_df The information data.frame
#'
#' @param txt_size The font size for text
#'
#' @param save_png Logical. Set TRUE if you would like to save the
#' architectures sequence logos as PNG files
#'
#' @param dir_path The path to the directory where plot will be saved.
#'
#'
#' @details Plots the sequence logos of all clusters
#'
#' @return If `one_plot` is TRUE, one plot as a grid of all sequence logos is
#' returned.
#'
#' If `one_plot` is FALSE, a set of ggplot2-based sequence logos for all
#' clusters is saved to disk with the default filename
#' `Architectures_0-max.pdf`. This is a multi-page PDF document with the
#' sequence logo for each cluster on a separate page. Also, the list of plots
#' is returned.
#'
#' @importFrom seqArchR plot_ggseqlogo_of_seqs plot_arch_for_clusters
#' collate_seqArchR_result get_seqs_clust_list seqs_str collate_clusters
#' viz_seqs_acgt_mat
#'
#' @export
#'
per_cluster_seqlogos <- function(sname, seqs = NULL, clusts,
                                 pos_lab = -45:45, bits_yax = "max",
                                 strand_sep = FALSE, one_plot = TRUE,
                                 info_df = NULL, txt_size = 12,
                                 save_png = FALSE, dir_path){
    if(strand_sep){
        if(is.null(info_df)) stop("`info_df` cannot be NULL when ",
                                  "`strand_sep` is TRUE")
        strand_sep_pl <- strand_sep_seqlogos(sname = sname, seqs = seqs,
                            clusts = clusts, info_df = info_df,
                            pos_lab = pos_lab, bits_yax = bits_yax,
                            dir_path = dir_path, txt_size = txt_size,
                            save_png = save_png)
        return(strand_sep_pl)
    }
    cli::cli_h1(paste0("All clusters' sequence logos"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##
    message("Generating architectures for clusters of sequences...")
    fname <- file.path(result_dir_path, paste0("Architectures_0-max.pdf"))
    ##
    T_F_titles <- TRUE
    if(one_plot) T_F_titles <- FALSE
    ## arch_list plots directly using seqArchR function
    arch_pl <- seqArchR::plot_arch_for_clusters(seqs = seqs,
                                clust_list = clusts,
                                pos_lab = pos_lab,
                                set_titles = T_F_titles,
                                method = "bits",
                                pdf_name = NULL, show = FALSE)

    arch_list <- lapply(seq_along(arch_pl), function(y) {
        pl <- arch_pl[[y]] + ggplot2::theme(
                    axis.text = ggplot2::element_text(size = txt_size),
                    axis.text.y = ggplot2::element_text(vjust = 0.5),
                    axis.title.y = ggplot2::element_text(size = txt_size),
                    axis.ticks.length = ggplot2::unit(0.04, "cm"),
                    plot.margin = ggplot2::unit(c(0.2,0.1,-0.4,0.1), "cm")) +
        ggplot2::scale_y_continuous(
                        sec.axis = ggplot2::dup_axis(name = paste0("C", y),
                        labels = NULL))
        pl
    })
    ## Making a plot_grid and plotting is tedious; the PDF height has to
    ## be set by trial and error to a very large value. Instead, we decide
    ## to plot this a pair on one page
    if(one_plot){
        comb_one_pl <- cowplot::plot_grid(plotlist = arch_list, ncol = 1)
        return(comb_one_pl)
    }else{
        grDevices::pdf(file = fname, width = 20, height = 2, onefile = TRUE)
        lapply(arch_list, print)
        grDevices::dev.off()

        ## save PNGs
        if(save_png)
            save_PNGs(dir_path = result_dir_path, plot_list = arch_list,
                  txt_size = txt_size)
        ##
        return(arch_list)
    }
    ##
}
##==============================================================================

save_PNGs <- function(dir_path, plot_list, strand_sep = FALSE, txt_size = 10){

    use_res_path <- file.path(dir_path, "arch_png")
    if(strand_sep) use_res_path <- file.path(dir_path, "strand_sep_arch_png")

    use_ht <- 3
    if(strand_sep) use_ht <- 6

    stopifnot(.check_and_create_dir(use_res_path))
    for(p in 1:length(plot_list)){
        fname <- file.path(use_res_path,
                           paste0("Architecture_clust", p, "_0-max.png"))
        pl <- plot_list[[p]] + ggplot2::theme(axis.text = ggplot2::element_text(size = txt_size),
                         axis.title.y = ggplot2::element_text(size = txt_size),
                         axis.text.y = ggplot2::element_text(size = txt_size))
        cowplot::ggsave2(fname, plot = pl, device = "svg",
                         width = 20, height = use_ht, units = "cm", dpi = 300)
    }
}
##==============================================================================





strand_sep_seqlogos <- function(sname, seqs, clusts, info_df, pos_lab,
                                bits_yax, dir_path, txt_size = 12,
                                save_png = FALSE){
    cli::cli_h1(paste0("All clusters' strand-separated sequence logos"))
    cli::cli_h2(paste0("Sample: ", sname))
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##

    ## STRAND SEPARTATED SEQUENCE LOGOS (AUTO LIMITS)
    fname <- file.path(result_dir_path,
                       paste0("Architectures_0-max_strand_separated.pdf"))


    per_cl_idx_by_strand <- lapply(clusts, function(x){
        idx_p_n <- lapply(c("+", "-"), function(y){
            .get_strand_specific_indices(df = info_df,
                                         seq_ids_in_cl = x,
                                         strand_val = y)
        })
        names(idx_p_n) <- c('+', '-')
        idx_p_n
    })


    per_cl_title_by_strand <- lapply(seq_along(clusts), function(x){

        title_p_n <- lapply(c("+", "-"), function(y){
            .get_strand_plot_title(x, nclust = length(clusts),
                                this_n = length(per_cl_idx_by_strand[[x]][[y]]),
                                tot_n = length(clusts[[x]]),
                                strand_val = y)
        })
        names(title_p_n) <- c('+', '-')
        title_p_n
    })

    plots_p_n <- lapply(seq_along(clusts), function(x){
        if(length(per_cl_idx_by_strand[[x]]$`+`) < 1){
            print("plus is zero")
            pl_p <- ggplot2::ggplot() + ggplot2::theme_void()
        }else{
            pl_p <- seqArchR::plot_ggseqlogo_of_seqs(
                seqs = seqs[per_cl_idx_by_strand[[x]]$`+`],
                pos_lab = pos_lab,
                xt_freq = 5,
                title = per_cl_title_by_strand[[x]]$`+`,
                bits_yax = bits_yax)
            pl_p <- pl_p + theme(
                axis.text = element_text(size = txt_size),
                axis.text.y = element_text(vjust = 0.5),
                axis.title.y = element_text(size = txt_size),
                axis.ticks.length = unit(0.04, "cm"),
                plot.margin = unit(c(0.2,0.1,-0.4,0.1), "cm")) +
                ggplot2::scale_y_continuous(
                    sec.axis = dup_axis(name = paste0("C", x),
                                        labels = NULL))
        }
        if(length(per_cl_idx_by_strand[[x]]$`-`) < 1){
            print("negative is zero")
            pl_n <- ggplot() + theme_void()
        }else{
            pl_n <- seqArchR::plot_ggseqlogo_of_seqs(
                seqs = seqs[per_cl_idx_by_strand[[x]]$`-`],
                pos_lab = pos_lab,
                xt_freq = 5,
                title = per_cl_title_by_strand[[x]]$`-`,
                bits_yax = bits_yax)
            pl_n <- pl_n + theme(
                axis.text = element_text(size = txt_size),
                axis.text.y = element_text(vjust = 0.5),
                axis.title.y = element_text(size = txt_size),
                axis.ticks.length = unit(0.04, "cm"),
                plot.margin = unit(c(0.2,0.1,-0.4,0.1), "cm")) +
                ggplot2::scale_y_continuous(
                    sec.axis = dup_axis(name = paste0("C", x),
                                        labels = NULL))
        }
        ##
        pl_p
        ##
        comb_pl_p_n <- cowplot::plot_grid(pl_p, pl_n, nrow = 2, align = "v")
        comb_pl_p_n
    })
    ## Making a plot_grid and plotting is tedious; the PDF height has to be set
    ## by trial and error to a very large value. Instead, we decide to plot this
    ## a pair on one page
    grDevices::pdf(file = fname, width = 20, height = 4.5, onefile = TRUE)
    lapply(plots_p_n, print)
    grDevices::dev.off()

    ## save PNGs
    if(save_png)
        save_PNGs(dir_path = result_dir_path, plot_list = plots_p_n,
              strand_sep = TRUE, txt_size = 10)


    return(plots_p_n)
}
##==============================================================================



#' @title Plot heatmaps of motifs occuring in seqArchR clusters
#'
#' @param sname The sample name
#'
#' @param seqs The sequences as a DNAStringSet object
#'
#' @param flanks Flank size. The same flank is used upstream and downstream.
#' A vector of values is also accepted when more than oone flanks should be
#' visualized.
#'
#' @param clusts List of sequence Ids in each cluster.
#'
#' @param use_colors Specify colors to use
#'
#' @param motifs A vector of motifs to be visualized in the sequence. This can
#' be any words formed by the
#' \href{https://www.bioinformatics.org/sms2/iupac.html}{IUPAC code}.
#' For example, TATAA, CG, WW, SS etc.
#'
#' @param dir_path The path to the directory
#'
#' @param fheight,fwidth,funits Height and width of the PDF file, and the units
#'  in which they are specified.
#'
#'
#' @return Nothing. PNG images are written to disk using the provided filenames.
#'
#' @importFrom Biostrings width
#'
#' @export
#'
plot_motif_heatmaps <- function(sname, seqs, flanks = c(50), clusts,
                            use_colors = NULL, motifs, dir_path,
                            fheight = 500, fwidth = 500, funits = "px"){
    cli::cli_h1(paste0("Motif heatmaps"))
    cli::cli_h2(paste0("Sample: ", sname))
    ##
    nClust <- length(clusts)
    if(is.null(use_colors)){
        message("Using default colors")
        nClust_colors <- get_ncolors(n = nClust, palname = "Set1")
    }
    ##
    all_motifs_str <- paste0(motifs, collapse="_")

    ###
    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
    ##
    ## Iterate over different flanks
    sam <- lapply(flanks, function(x){
        cli::cli_h3(paste0("Using flank: ", x))
        fl_up <- fl_down <- x

        maxL <- Biostrings::width(seqs[1])
        midP <- base::ceiling(maxL/2)
        seqs <- Biostrings::subseq(seqs, start = midP - fl_up,
                                   end = midP + fl_down)
        ##
        seq_order <- unlist(clusts)
        #### Using heatmaps pkg
        patt_hm_list500 <- lapply(motifs,
                          function(x){
                              hm <- heatmaps::PatternHeatmap(
                                  seq = seqs[seq_order],
                                  pattern = x,
                                  coords = c(-1*fl_up, fl_down),
                                  label = x)
                              shm <- heatmaps::smoothHeatmap(hm, sigma=c(2,2))
                          }
        )
        ##
        fname_suffix <- paste0(paste(fl_up, "up", fl_down, "down",
                                     "motifHeatmaps", all_motifs_str,
                                     sep = "_"), ".png")
        fname <- file.path(result_dir_path, fname_suffix)

        clust_lens <- unlist(lapply(clusts, length))
        ##
        grDevices::png(fname, height = fheight, width = fwidth, units = funits)
        pl_hms <- heatmaps::plotHeatmapList(patt_hm_list500,
                                        box.width = 1.3,
                                        cex.label = 1.1,
                                        cex.axis = 0.7, scale.lwd = 0.5,
                                        tcl=-0.3, padj=-1.5,
                                        partition = clust_lens,
                                        partition.legend = TRUE,
                                        partition.lines = TRUE,
                                        partition.col = nClust_colors,
                                        legend = TRUE,
                                        legend.width = 0.3, cex.legend = 0.8,
                                        legend.pos = "r")
        grDevices::dev.off()
        cli::cli_alert_success(paste0("Written to: ", fname))
    })
    ##
}
## =============================================================================


## Use this function to get n colors (in sequence) from the specified palette
## -- Can specify n greater than that existing in a palette, in which case
## the colors are recycled and returned
## -- Can also specify a list of colors of choice < n, which are recycled and
## returned. Useful when random colors from a palette are required
get_ncolors <- function(n, palname = "Set1", clrs = NULL){

    ## Recycle color vector from RColorBrewer
    use_colors <- clrs
    if(is.null(clrs)){
        use_colors <- suppressWarnings(
            RColorBrewer::brewer.pal(n = n, name = palname))
    }

    nColor <- length(use_colors)
    if(n <= nColor){
        n_colors <- use_colors[seq_len(n)]
        return(n_colors)
    }

    rep_times <- base::ceiling((n-nColor)/nColor)
    if(n %% nColor == 0) rep_times <- rep_times + 1
    additional <- ((n-nColor) %% nColor)
    col_idx <- c(rep(seq_len(nColor), rep_times), seq_len(additional))
    n_colors <- use_colors[col_idx]
    n_colors
}
## =============================================================================

#' @title Visualize as barplots how promoters in each cluster are distributed
#' on different chromosomes and strands
#'
#' @param sname The sample name
#' @param clusts List of sequence ids in each cluster.
#' @param info_df The information data.frame
#' @param dir_path Specify the path to the directory on disk where plots
#' will be saved
#' @param colrs Specify colors used for two strands
#'
#' @importFrom ggplot2 theme_classic scale_fill_manual element_line
#'
#' @export
#'
per_cluster_strand_dist <- function(sname, clusts, info_df, dir_path,
                                    colrs = "Paired"){
    cli::cli_h1(paste0("All clusters' strand distributions"))
    cli::cli_h2(paste0("Sample: ", sname))

    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)


    fname <- file.path(result_dir_path,
                paste0("Per_cluster_strand_distributions.pdf"))

    chr_set <- unique(info_df$chr)
    nChr <- length(chr_set)
    empty_df <- data.frame("Chromosomes" = rep(chr_set, each=2),
                           "Strand" = rep(c("-", "+"), times = nChr),
                           "Frequency" = rep(0, times = nChr))
    ## Using 'times' instead of 'each' for Strand column and
    ## 'each' instead of 'times in Chromosomes column
    ## in empty_df makes
    ## it easier when iterating to pick each chromosome for pie_list

    ## Used later when + and - frequencies are zero
    # dummy_df <- empty_df
    # dummy_df[, "Frequency"] <- 5
    # ##
    chrStr <- paste0(empty_df$Chromosomes, empty_df$Strand)

    plot_list <- lapply(seq_along(clusts), function(x){
        df_strand <- data.frame("chr" = info_df$chr[ clusts[[x]] ],
                                "strand" = info_df$strand[ clusts[[x]] ] )
        df_strand_tab <- table(df_strand)
        df_strand_df <- data.frame(df_strand_tab)
        RchrStr <- paste0(df_strand_df$chr, df_strand_df$strand)

        filled_df <- empty_df

        filled_df$chrStr <- chrStr
        fill_idx <- unlist(lapply(RchrStr, function(a){
            which(filled_df$chrStr == a)
        }))

        filled_df$Frequency[fill_idx] <- df_strand_df$Freq

        ## Barplots
        pl <- ggplot(filled_df) +
            ggplot2::geom_col(aes(x=Chromosomes, y = Frequency, fill = Strand,
                                  width = 0.5),
                              position = ggplot2::position_dodge(width = 0.5)) +
            scale_fill_manual(values=colrs) +
            xlab("") +
            ggtitle(label = .get_strand_plot_title(this_id = x,
                                nclust = length(clusts),
                                tot_n = length(clusts[[x]]),
                                strand_val = NULL)) +
            theme_classic() +
            theme(legend.text = element_text(size = 12, face = "plain",
                                             color = "black"),
                  legend.position = "right",
                  panel.grid.major.y = element_line()
                  )
        return(pl)

    })
    grDevices::pdf(file = fname, width = 20, height = 2, onefile = TRUE)
    lapply(plot_list, print)
    grDevices::dev.off()
    plot_list
}
## =============================================================================

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
#' #' TO-DO What if two elements from independent clusters are to
#' be put into a totally new cluster, will using 0 work in this case?
#'
#' @param final Logical, set to TRUE or FALSE
#'
#' @param dir_path The path to the directory where files will be written
#'
#' @details
#'
#'
#'
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
#' The return value of the function is the clustering with re-assignments
#' executed.
#'
#' @importFrom stats cutree median
curate_clusters <- function(sname, use_aggl = "ward.D", use_dist = "euclid",
                            seqArchR_result, iter, pos_lab = NULL,
                            regularize = TRUE, topn = 50, use_cutk = 2,
                            need_change = NULL, change_to = NULL,
                            final = FALSE, dir_path){
    cli::cli_h1(paste0("seqArchR result clusters curation"))
    cli::cli_h2(paste0("Sample: ", sname))

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


#' @title Order clusters by median or mean interquantile widths
#'
#' @param sname The sample name
#' @param clusts List of sequence ids in each cluster.
#' @param info_df The data.frame with all tag clusters information. The
#' following columns are expected in the data.frame:"chr", "start", "end",
#' "width", "strand", "score", "nr_ctss", "dominant_ctss", "domTPM",
#' "IQW", "tpm" and two additional columns based on qLow and qUp used.#'
#'
#' @param order_by_median Logical. Whether to order to clusters by their
#' median (when TRUE) or mean (when FALSE) interquantile widths.
#'
#' @importFrom stats median
#'
#' @export
#'
order_clusters_iqw <- function(sname, clusts, info_df,
                               order_by_median = TRUE){
    cli::cli_h1(paste0("Order clusters by IQW"))
    cli::cli_h2(paste0("Sample: ", sname))
    cluster_medians_IQW <- unlist(lapply(clusts, function(x){
        stats::median(info_df$IQW[x])
    }))
    cluster_means_IQW <- unlist(lapply(clusts, function(x){
        base::mean(info_df$IQW[x])
    }))

    if(order_by_median){
        ascending_order_IQW <- sort(cluster_medians_IQW, decreasing = FALSE,
                                    index.return = TRUE)
    }else{
        ascending_order_IQW <- sort(cluster_means_IQW, decreasing = FALSE,
                                    index.return = TRUE)
    }
    ##
    clusts_list_ordered <- lapply(ascending_order_IQW$ix, function(x){
                   clusts[[x]]
               })

    clusts_list_ordered
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
        ## clusters, i.e., new clusters. This are marked by numeral 0 in
        ## change_to. Instead of leting them create confusions, we first
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
                                       rel_widths = c(0.40,1,1), rel_heights = c(1,1,1),
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
#'
#' @export
#'
#'
seqs_acgt_image <- function(sname, seqs, seqs_ord, pos_lab, xt_freq, yt_freq,
                            f_height, f_width, dir_path){

    result_dir_path <- .handle_per_sample_result_dir(sname, dir_path)
        fname <- file.path(result_dir_path,
                           paste0(sname, "_ClusteringImage.png"))
        seqArchR::viz_seqs_acgt_mat(as.character(seqs[unlist(seqs_ord)]),
                                    pos_lab = pos_lab,
                                    xt_freq = 5,
                                    yt_freq = 500,
                                    f_height = 1200, f_width = 600,
                                    save_fname = fname,
                                    file_type = "png")
}
## =============================================================================



## 1. Get PFM matrix with custom background frequencies for the DNA_BASES
## 2. Use TFBSTools::toICM to convert to information content matrix
## 3. Visualize with ggseqlogo

# enrichment_seqlogos <- function(seqs, clusts){
#
#     nonuni_bg_seqlogos <- lapply(seq_along(clusts), function(x){
#         print(paste0("Clust ", x))
#         seqs_bg <- Biostrings::letterFrequency(seqs[clusts[[x]]],
#                             letters = Biostrings::DNA_BASES, as.prob = TRUE)
#         seqs_bg <- colSums(seqs_bg)/nrow(seqs_bg)
#         print(seqs_bg)
#         pfm_mat <- Biostrings::consensusMatrix(seqs[clusts[[x]]],
#                                         baseOnly = TRUE)[1:4,]
#
#         seqs_pfm <- TFBSTools::PFMatrix(bg = seqs_bg,
#                                         profileMatrix = pfm_mat)
#
#         seqs_pwm <- TFBSTools::toPWM(seqs_pfm, bg = seqs_bg, type = "prob")
#
#         pl <- ggseqlogo::ggseqlogo(data = seqs_pwm@profileMatrix)
#         pl
#     })
#     return(nonuni_bg_seqlogos)
# }

