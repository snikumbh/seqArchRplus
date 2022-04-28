sklearn <- NULL
.onLoad <- function(libname, pkgname) {

    #
    # CMD check NOTE avoidance See: https://stackoverflow.com/a/12429344
    # Example:
    # https://github.com/HughParsonage/grattan/blob/master/R/zzz.R
    if (getRversion() >= "2.15.1") {
        utils::globalVariables(
            c(
                "Feature",
                "Frequency",
                "Chromosomes",
                "Strand",
                "clust",
                "clust_ID",
                "start",
                "seqnames",
                "strand",
                "end",
                "IQW",
                "domTPM"
            )
        )
    }
}
