

# seqArchRplus 1.1.2

## New
* (User-facing) Ability to generate HTML reports to view large combined panels 
  for multiple processed samples as scrollable carousels

# seqArchRplus 0.99.10

## New
* (User-facing) Alternative function `plot_motif_heatmaps2()` combines the 
  heatmap scale with the heatmaps. Also, can choose to output PNG or JPEG 
  files


# seqArchRplus 0.99.0

## New
* Package ready for Bioconductor submission

# seqArchRplus 0.99.0.10

## New
* (User-facing) An alternative function to plot motif heatmaps using the R pkg 
  seqPattern. Function name: `plot_motif_heatmaps2()`


# seqArchRplus 0.99.0.7

## New
* (User-facing) Now possible to perform per cluster GO term enrichment analysis 
  via `per_cluster_go_term_enrichments()` when a orgDb package is available for
  the organism


# seqArchRplus 0.99.0.6

## New
* (User-facing) Examples added majority of functions using example/dummy data

# seqArchRplus 0.99.0.5

## Bug-fixes
* (User-facing) In `per_cluster_annotations()`:
  * Default for `clusts` set to NULL. This enables cleanly specifying just the 
  BED file as input to `tc_gr`
  * Help pages now have examples except for the `curate_clusters` function. 
  Elaborate documentation for this function is available as part of the 
  vignette

# seqArchRplus 0.99.0.4

## New features
* (User-facing) In `write_seqArchR_cluster_track_bed()`: 
  * new argument `use_q_bound` to choose if you wish to use quantiles as 
  tag cluster boundaries
  * new argument `use_as_names` to specify any column in `info_df` to be used 
  as the name column in the output track BED file of clusters
  * dominant_ctss information for each tag cluster is presented as thickStart 
  thickEnd for ease of visualising. See documentation for more details
* (User-facing) In `per_cluster_annotations()`:
  * tc_gr can now accept a bedfile to read records as a GRanges object
  * Details added in documentation for ways to selectively pick clusters to 
  annotate



# seqArchRplus 0.99.0.3

## Bux-fixes
* Fixed bugs in `curate_clusters()` function, and touch up its documentation


# seqArchRplus 0.99.0.2

## New features
* (User-facing) Parallelization support to speed up annotating genomic regions
