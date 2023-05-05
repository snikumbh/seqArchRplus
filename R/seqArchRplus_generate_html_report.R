

#' @title Generate HTML report with scrollable combined panel plots
#'
#' @description This function generates an HTML report with large scrollable
#' combined panels for multiple samples that eases comparison of changes
#' between samples
#'
#' @param snames Sample names to be included
#'
#' @param file_type "PDF" (default) or "SVG" (to be supported in the future).
#' The type of files to look for in the sample-specific results folder
#'
#' @param img_ht,img_wd The height and width (in pixels) of the images in the
#' output HTML report. Default values are '1200px' (height) and
#' '1600px' (width)
#'
#' @param page_wd The width of the body in the HTML. Default is '1800px'.
#'
#' @param render_silently Logical. TRUE or FALSE
#'
#' @param dir_path Specify the `/path/to/directory` where sample-specific
#' results folders are located. This is a required argument and cannot be NULL.
#' A directory named `combined_results` is created at the given location, and
#' the HTML report is written into it
#'
#' @details This functionality requires suggested libraries \code{
#' `slickR`} and \code{`pdftools`} installed.
#' The function assumes requires that the candidate figure files have combined_panel
#' Note that the combined plot panels are arranged horizontally and
#' therefore are best viewed in wide desktop monitors.
#'
#' @return Nothing. Report is written to disk at the provided \code{`dir_path`
#' } using the filename \code{
#' 'Combined_panels_report_samples_<samples_names>.html'}.
#'
#' @export
#'
#' @examples
#'
#' ## Need these packages to run these examples
#' if(require("slickR", 'pdftools")){
#'
#' ## Make IQW-TPM plots
#'
#' bed_fname <- system.file("extdata", "example_info_df.bed.gz",
#'          package = "seqArchRplus", mustWork = TRUE)
#'
#' info_df <- read.delim(file = bed_fname,
#'          sep = "\t", header = TRUE,
#'          col.names = c("chr", "start", "end", "width",
#'                  "dominant_ctss", "domTPM",
#'                  "strand", "score", "nr_ctss",
#'                  "q_0.1", "q_0.9", "IQW", "tpm"))
#'
#' use_clusts <- readRDS(system.file("extdata", "example_clust_info.rds",
#'          package = "seqArchRplus", mustWork = TRUE))
#'
#' use_dir <- tempdir()
#'
#' iqw_tpm_pl <- iqw_tpm_plots(sname = "sample1",
#'                             dir_path = use_dir,
#'                             info_df = info_df,
#'                             iqw = TRUE,
#'                             tpm = TRUE,
#'                             cons = FALSE,
#'                             clusts = use_clusts,
#'                             txt_size = 14)
#'
#' ## Make sequence logos
#' library(Biostrings)
#' raw_seqs <- Biostrings::readDNAStringSet(
#'                           filepath = system.file("extdata",
#'                             "example_promoters45.fa.gz",
#'                             package = "seqArchRplus",
#'                             mustWork = TRUE)
#'                         )
#'
#'
#' seqlogo_oneplot_pl <- per_cluster_seqlogos(sname = "sample1",
#'                                    seqs = raw_seqs,
#'                                    clusts = use_clusts,
#'                                    pos_lab = -45:45,
#'                                    bits_yax = "auto",
#'                                    strand_sep = FALSE,
#'                                    one_plot = TRUE,
#'                                    dir_path = use_dir,
#'                                    txt_size = 14)
#'
#' ## Need the TxDb object to run these examples
#' if(require("TxDb.Dmelanogaster.UCSC.dm6.ensGene")){
#'     annotations_pl <- per_cluster_annotations(sname = "sample1",
#'                          clusts = NULL,
#'                          tc_gr = bed_fname,
#'                          txdb_obj = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
#'                          one_plot = FALSE,
#'                          dir_path = use_dir,
#'                          tss_region = c(-500,100))
#' }else{
#'     annotations_pl <- NULL
#' }
#'
#' # Combine them together
#' if(!is.null(annotations_pl)){
#'     panel_pl <- form_combined_panel(iqw_tpm_pl = iqw_tpm_pl,
#'                     seqlogos_pl = seqlogos_oneplot_pl,
#'                     annot_pl = annotations_oneplot_pl)
#'  }else{
#'     panel_pl <- form_combined_panel(iqw_tpm_pl = iqw_tpm_pl,
#'                     seqlogos_pl = seqlogos_oneplot_pl)
#'
#'  }
#'
#'  cowplot::save_plot(filename = file.path(use_dir,
#'                                         paste0("sample1_combined_panel.pdf"),
#'                    plot = panel_pl)
#'
#' # Call function to generate HTML report
#' generate_html_report(snames = c("sample1", "sample1"),
#'                 dir_path = use_dir)
#' }
#'
#' @author Sarvesh Nikumbh
generate_html_report <- function(snames, file_type = "PDF",
                                img_ht = '1200px', img_wd = '1600px',
                                page_wd = '1800px',
                                render_silently = TRUE, dir_path) {

    if(!requireNamespace(c("slickR", "pdftools"), quietly = TRUE)){
        stop("Please install R package(s) 'slickR' and 'pdftools' to use ",
            "this functionality")
    }else{

        dir_path <- normalizePath(dir_path)
        this_dir_path <- file.path(dir_path, paste0(snames, "_results"))

        ## Fail if the sample-specific results folders are not available.
        if(any(!dir.exists(this_dir_path))){
            stop("Did not find sample-specific results folder(s):",
                which(!dir.exists(this_dir_path)), call. = FALSE)
        }


        result_dir_path <- file.path(dir_path, paste0("combined_results"))

        if(!dir.exists(result_dir_path)){
            message("Creating folder: ", result_dir_path)
            dir.create(result_dir_path)
        }else{
            message("Results directory exists: ", result_dir_path)
        }
        ##

        ## This will change to incorporate file formats in the future.
        ## Current plan is to support SVG next, so that high-quality figures are
        ## produced.
        extn <- ifelse(file_type == "PDF", ".pdf", ".pdf")

        pdf_files <- list.files(this_dir_path,
                        pattern = paste0("*", "combined_panel", extn),
                        full.names = TRUE)

        if(length(pdf_files) < length(snames))
            warning("Some samples are missing comnbined panel plot PDFs",
                immediate. = TRUE, call. = FALSE)

        use_template <- system.file("extdata",
            "combined-panel-report-template.Rmd",
            package = "seqArchRplus")

        dest_template <- file.path(result_dir_path, basename(use_template))

        file.copy(from = use_template,
            to = dest_template,
            overwrite = TRUE)

        message("Generating report...", appendLF = FALSE)

        html_name <- file.path(result_dir_path,
            paste0("Combined_panels_report_samples_",
            paste(snames, collapse = "_"), ".html"))


        rmarkdown::render(dest_template, output_dir = result_dir_path,
            quiet = render_silently,
            output_file = html_name,
            params = list(this_dir_path = this_dir_path,
                            snames = snames,
                            img_ht = img_ht, img_wd = img_wd,
                            page_wd = page_wd,
                            lsnames = length(snames),
                            snames_str =  paste(snames, collapse = ", ")
                        )
            )

        file.remove(dest_template)
        message("Done", appendLF = TRUE)

        message("Output HTML file: ", html_name)
    }
}
## =============================================================================
