## Write CAGE tagclusters to disk
##

#' @title Handle writing of tag clusters (TCs) from CAGE experiment to disk as
#'   BED files.
#'
#' @description Handle writing of tag clusters obtained from a CAGE experiment
#'   to disk as BED files. It can also return corresponding promoter sequences
#'   as a DNAStringSet object or write them to disk as FASTA files.
#'
#'   For information on tag clusters obtained by clustering of CAGE-derived TSS,
#'   we refer the reader to the CAGEr vignette, Section 3.7
#'   \url{https://www.bioconductor.org/packages/release/bioc/vignettes/CAGEr/inst/doc/CAGEexp.html#ctss-clustering}
#'
#'
#'
#' @param sname The sample name
#'
#' @param tc_gr Tag clusters as \code{\link[GenomicRanges]{GRanges}}. If
#'   `cager_obj` is not provided (is NULL), this argument is required. Also note
#'   that the columns should match the columns (with mcols) as returned by
#'   \code{\link[CAGEr]{tagClusters}}. If `tc_gr` is provided, `cager_obj` is
#'   ignored
#'
#' @param cager_obj A CAGEexp object obtained from the CAGEr package, if and
#'   when CAGEr was used to process the raw CAGE data
#'
#' @param qLow,qUp The interquantile boundaries to be considered for obtaining
#'   tag clusters from the CAGEexp object. See \code{\link[CAGEr]{tagClusters}}
#'
#' @param fl_size_up,fl_size_down Numeric. The size of the flanks in the
#'   upstream and downstream directions.
#'
#' @param bsgenome The BSgenome file that will be used to obtain sequences of
#'   the organism from.
#'
#' @param dir_path The path to the directory where files will be written. By
#'   default, all BED files are written within a subdirectory named "BED", and
#'   all FASTA files are written within a subdirectory named "FASTA", both
#'   created at the `dir_path` location.
#'
#' @param fname_prefix,fname_suffix Specify any prefix or suffix string to be
#'   used in the filename. This can be the organism name etc. Specify without
#'   any separator. By default, an underscore is used as a separator in the
#'   filename.
#'
#' @param write_to_disk Logical. Specify TRUE to write files to disk. More
#'   specifically, BED files are written to disk only when this is set to TRUE.
#'   For promoter sequences, FASTA files are written to disk if this arg is set
#'   to TRUE, otherwise not. and a DNAStringSet object is returned if `ret_seqs`
#'   is set to TRUE.
#'
#' @param ret_seqs Logical. Specify TRUE if promoter sequences are to be
#'   returned as a DNAStringSet object.
#'
#' @return If `ret_seqs` is TRUE, a DNAStringSet object is returned. Depending
#'   on `write_to_disk`, files are written to disk at the specified location.
#'
#' @details You can use the fname_prefix and fname_suffix arguments to specify
#' strings to be used as prefix and suffix for the files. For example, the
#' organism name can be used as a prefix for the filename. Similarly, for
#' suffix.
#'
#' @return If `ret_seq = TRUE`, the promoter sequences are returned as a
#' \code{\link[Biostrings]{DNAStringSet}} object.
#'
#' @importFrom GenomicRanges promoters trim
#' @importFrom IRanges IRanges
#' @importFrom GenomeInfoDb seqlevels seqlevels<- seqinfo seqinfo<-
#' @importFrom S4Vectors mcols
#' @importFrom utils read.delim
#'
#' @export
#'
#' @examples
#'
#' if(require("BSgenome.Hvulgare2.NCBI.MorexV3")){
#'
#'   seqs <- seqArchRplus::handle_tc_from_cage(sname = sn,
#'                                  cager_obj = ce,
#'                                  fl_size_up = 500,
#'                                  fl_size_down = 500,
#'                                  dir_path = NULL,
#'                                  fname_prefix = "pre",
#'                                  fname_suffix = "suff",
#'                                  write_to_disk = FALSE,
#'                                  bsgenome = BSgenome.Hvulgare2.NCBI.MorexV3,
#'                                  ret_seqs = TRUE)
#'
#' }
#'
#'
## Writes TCs as BED files via GRanges and also as FASTA files with a fixed
## flank size of +/- 500 bp (or provided flank size)
##
##
##
handle_tc_from_cage <- function(sname, tc_gr, cager_obj, qLow = 0.1, qUp = 0.9,
                                fl_size_up = 500, fl_size_down = 500,
                                bsgenome, dir_path = NULL,
                                fname_prefix = NULL, fname_suffix = NULL,
                                write_to_disk = TRUE, ret_seqs = TRUE) {
    cli::cli_h1(paste0("Handling TCs from CAGE (BED and FASTA)"))
    cli::cli_h2(paste0("Sample: ", sname))

    bed_fname <- file.path(
        dir_path, "BED",
        paste0(
            fname_prefix, "_TC_sample_",
            sname, "_", fname_suffix, ".bed"
        )
    )
    fasta_fname <- file.path(
        dir_path, "FASTA",
        paste0(
            fname_prefix, "_sample_",
            sname, "_promoters_",
            fname_suffix, ".fa"
        )
    )
    ##
    this_gr <- .handle_tc_cager(tc_gr, cager_obj, sname, qLow, qUp)
    ##
    # this_gr <- CAGEr::tagClustersGR(cager_obj, sample = sname,
    #                                 returnInterquantileWidth = TRUE,
    #                                 qLow = 0.1, qUp = 0.9)
    ## At present CAGEr tagClusters loose this information, so add them
    # seqlevels(this_gr) <- GenomeInfoDb::seqlevels(
    #     CAGEr::CTSStagCountSE(cager_obj))
    # seqinfo(this_gr)   <- GenomeInfoDb::seqinfo(
    #     CAGEr::CTSStagCountSE(cager_obj))

    if (any("score" == names(S4Vectors::mcols(this_gr)))) {
        this_gr$tpm <- as.numeric(this_gr$score)
    } else {
        stop("Expecting an mcol named 'score', did not find.")
    }
    this_gr_df <- as.data.frame(this_gr)

    ## Storing dominant CTSS information into GRanges object...
    this_gr_domCtss <- GenomicRanges::GRanges(
        seqnames = this_gr_df[, "seqnames"],
        ranges = IRanges::IRanges(
            start = this_gr_df[, "dominant_ctss"],
            width = 1
        ),
        strand = this_gr_df[, "strand"]
    )



    this_prom <- GenomicRanges::promoters(this_gr_domCtss,
        upstream = fl_size_up,
        downstream = fl_size_down + 1
    )
    this_prom <- GenomicRanges::trim(this_prom, use.names = TRUE)
    cli::cli_alert_success("Making promoters")

    ## Exclude indices where the sequences have been trimmed
    omit_ids <- which(Biostrings::width(this_prom) != ((fl_size_up +
        fl_size_down) + 1))

    message("Omitted IDs: ", length(omit_ids))

    if (length(omit_ids) > 0) {
        this_gr <- this_gr[-omit_ids, ]
        this_prom <- this_prom[-omit_ids]
        this_gr_domCtss <- this_gr_domCtss[-omit_ids]
    }

    if (write_to_disk) .write_tc_bed(this_gr, bed_fname)
    cli::cli_alert_success("BED file written")
    prom_seqs <- .write_tc_fasta(this_prom,
        bsgenome = bsgenome,
        fasta_fname = fasta_fname,
        ret_val = ret_seqs,
        also_write = write_to_disk
    )
    cli::cli_alert_success("FASTA file written")
    if (ret_seqs) {
        return(prom_seqs)
    }
}
## =============================================================================

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
#' @keywords internal
#' @noRd
.write_tc_bed <- function(gr, bed_fname) {
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
        row.names = FALSE
    )
}
## =============================================================================


.write_tc_fasta <- function(prom, fl_size_up = 500, fl_size_down = 500,
                            bsgenome, fasta_fname, ret_val = TRUE,
                            also_write = TRUE) {
    cli::cli_alert_info("Fetching sequences...")

    fasta_names <- paste0(
        "domCTSS=", seqnames(prom), ":",
        start(prom), ";strand=",
        strand(prom), ";",
        "up=", fl_size_up, ";",
        "down=", fl_size_down
    )

    prom_seqs <- BSgenome::getSeq(bsgenome,
        names = seqnames(prom),
        start = start(prom),
        end = end(prom),
        strand = strand(prom)
    )
    names(prom_seqs) <- fasta_names
    ##
    if (also_write) {
        cli::cli_alert_info(paste0("Writing FASTA file at: ", fasta_fname))
        Biostrings::writeXStringSet(prom_seqs,
            filepath = fasta_fname,
            format = "FASTA"
        )
    }
    if (ret_val) {
        return(prom_seqs)
    }
}
## =============================================================================
