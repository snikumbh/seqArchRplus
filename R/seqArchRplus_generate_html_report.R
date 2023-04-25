

#' @title Generate HTML report with scrollable combined panel plots
#'
#' @description This function generates an HTML report with large scrollable
#' combined panels for multiple samples that eases comparison of changes
#' between samples
#'
#' @param snames Sample names to be included
#'
#' @param file_type "PDF" (default) or "SVG". The type of files to look for
#' in the sample-specific results folder
#'
#' @param img_ht,img_wd The height and width (in pixels) of the images in the
#' output HTML report. Default values are '1200px' (height) and
#' '1600px' (width)
#'
#' @param page_wd The width of the body in the HTML. Default is '1800px'.
#'
#' @param render_silently Logical. TRUE or FALSE
#'
#' @param dir_path Specify the /path/to/directory where sample-specific results
#' folders are located. This is a required argument and cannot be NULL. A
#' directory named `combined_results` is created at the given location, and the
#' HTML report is written into it
#'
#' @details This functionality requires suggested libraries \code{
#' `slickR`} and \code{`pdftools`} installed. Note that the combined
#' plot panels are arranged horizontally and therefore are best viewed in wide
#' desktop monitors.
#'
#' @return Nothing. Report is written to disk at the provided \code{`dir_path`
#' } using the filename \code{
#' 'Combined_panels_report_samples_<samples_names>.html'}.
#'
#' @export
#'
#' @examples
#'
#' use_dir <- tempdir()
#'
#'
#' generate_html_report(snames = c("sample1", "sample2"),
#'                 dir_path = use_dir)
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

        result_dir_path <- file.path(dir_path, paste0("combined_results"))

        if(!dir.exists(result_dir_path)){
            message("Creating folder: ", result_dir_path)
            dir.create(result_dir_path)
        }else{
            message("Results directory exists: ", result_dir_path)
        }
        ##

        pdf_files <- list.files(this_dir_path, pattern = "*combined_panel.pdf",
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
