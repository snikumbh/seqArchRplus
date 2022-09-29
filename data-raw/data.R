## This script produces the reduced version of the data to be used as a minimal
## example with the seqArchRplus package
##
## Data:
## The required data (from which a reduced version is created) is stored in the
## data-raw folder in the seqArchRplus package source directory, only in a
## separate dedicated branch named "example-data-raw" available on
## seqArchRplus' github link
##
read_path <- file.path(".")
write_path <- file.path("../inst/extdata")

prom <- Biostrings::readDNAStringSet(
    filepath = file.path(read_path,
    "dm6_samarth_schor_et_al_TC_sample_RAL28_10_to_12_sample12_minTPM1_flank_up500_flank_down500.fa"))
prom_small <- Biostrings::subseq(prom, start = 501-45, end = 501+45)
prom_small_200 <- Biostrings::subseq(prom, start = 501-200, end = 501+200)

## Pick 8, 9, 20, 25, 28 from 10-12h Schor et al result
seqArchR_result <- readRDS(file.path(read_path, "schor2017_archRresult.rds"))
seqArchR_clusts <- seqArchR::get_seqs_clust_list(
                            seqArchR_result$seqsClustLabels[[5]])
clust_idx_iter5 <- c(8,9,20,25,28)
pick_clusts <- lapply(clust_idx_iter5, function(x) seqArchR_clusts[[x]])
pick_idx <- unlist(pick_clusts)
print(length(pick_idx))

tc_gr <- readRDS(file.path(read_path,
    "dm6_samarth_schor_et_al_TC_sample_RAL28_10_to_12_sample12_minTPM1.rds"))

pick_idx <- unlist(pick_clusts)


pick_prom_small_200 <- prom_small_200[pick_idx]
pick_prom_small_45 <- prom_small[pick_idx]


## Make info_df from tagClusters granges object
info_df <- GenomicRanges::as.data.frame(tc_gr)
info_df$tpm <- info_df$score
colnames(info_df) <- c("chr", "start", "end", "width",
    "strand",	"score", "nr_ctss",
    "dominant_ctss", "domTPM",
    "q_0.1", "q_0.9",
    "IQW", "tpm")

info_df <- info_df[pick_idx,]
tc_gr <- tc_gr[pick_idx,]

## info_df rownames are sequence IDs from the original set of sequences
## (the parent fasta from which only subsets have been selected)
## While pick_clusts holds these sequence IDs, sam_upd reinitializes the clustIDs
## such that the row numbers in info_df are used and not the rownames
##
sam_lens <- cumsum(lengths(pick_clusts))
sam_upd <- lapply(seq_along(pick_clusts), function(x) {
    if(x == 1) {seq(length(pick_clusts[[x]]))}
    else {sam_lens[x-1] + seq(length(pick_clusts[[x]]))}
}
)

## Uncomment to write to disk
saveRDS(sam_upd, file = "../inst/extdata/example_clust_info.rds")


## Write fasta files to disk
## ## ## Uncomment to write to disk
Biostrings::writeXStringSet(pick_prom_small_200,
    filepath = file.path(write_path, "example_promoters200.fa.gz"),
    compress = TRUE)

## ## Uncomment to write to disk
Biostrings::writeXStringSet(pick_prom_small_45,
    filepath = file.path(write_path, "example_promoters45.fa.gz"),
    compress = TRUE)

## Write gzipped info_df to disk
## ## Uncomment to write to disk
gzfile <- gzfile(file.path(write_path, "example_info_df.bed.gz"), "w")
write.table(info_df, file = gzfile, sep = "\t")
close(gzfile)

## Write tc_gr to disk (as retuened by CAGEr, but just the selected pick_idx)
## Uncomment to write to disk
saveRDS(tc_gr, file = file.path(write_path, "example_tc_gr.rds"))


## Modify seqArchR_result to retain only required data from iter 5 and remove
## everything else

new_result <- seqArchR_result
new_result$clustBasisVectors[[1]]$basisVectors <- NULL
new_result$clustBasisVectors[[2]]$basisVectors <- NULL
new_result$clustBasisVectors[[3]]$basisVectors <- NULL
new_result$clustBasisVectors[[4]]$basisVectors <- NULL
new_result$clustBasisVectors[[5]]$basisVectors <-
    as.matrix(new_result$clustBasisVectors[[5]]$basisVectors[, clust_idx_iter5])

## Handle sequence cluster labels

old_labels <- unique(new_result$seqsClustLabels[[5]][pick_idx])
new_seqClustLabels <- lapply(seq_along(old_labels), function(x){
    rep(x, length(which(new_result$seqsClustLabels[[5]] == old_labels[x])))
})
new_result$seqsClustLabels[[5]] <- unlist(new_seqClustLabels)

## Handle raw sequences
new_result$rawSeqs <- seqArchR_result$rawSeqs[pick_idx]

## Uncomment to write to disk
saveRDS(new_result, file =  file.path(write_path, "seqArchR_result.rds"))
