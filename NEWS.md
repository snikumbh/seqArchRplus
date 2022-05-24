# seqArchRplus 0.99.0.4

## New features
* (User-facing) In `write_seqArchR_cluster_track_bed()`: 
  * new argument `use_q_bound` to choose if you wish to use quantiles as 
  tag cluster boundaries
  * new argument `use_as_names` to specify any column in `info_df` to be used 
  as the name column in the output track BED file of clusters
  * dominant_ctss information for each tag cluster is presented as thickStart 
  thickEnd for ease of visualising. See documentation for more details.
* (User-facing) `per_cluster_annotations()`:
  * tc_gr can now accept a bedfile to read records as a GRanges object
  * Details added in documentation for ways to selectively pick clusters to 
  annotate



# seqArchRplus 0.99.0.3

## Bux-fixes
* Fixed bugs in `curate_clusters()` function, and touch up its documentation


# seqArchRplus 0.99.0.2

## New features
* (User-facing) Parallelization support to speed up annotating genomic regions
