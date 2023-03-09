<!-- badges: start -->
<!--  [![Bioc release status](http://www.bioconductor.org/shields/build/release/bioc/seqArchRplus.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/seqArchRplus) -->
  [![Bioc devel status](http://www.bioconductor.org/shields/build/devel/bioc/seqArchRplus.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/seqArchRplus)
  [![Bioc downloads rank](https://bioconductor.org/shields/downloads/devel/seqArchRplus.svg)](http://bioconductor.org/packages/stats/bioc/seqArchRplus/)
  [![Bioc support](https://bioconductor.org/shields/posts/seqArchRplus.svg)](https://support.bioconductor.org/tag/seqArchRplus)
  [![Bioc history](https://bioconductor.org/shields/years-in-bioc/seqArchRplus.svg)](https://bioconductor.org/packages/release/bioc/html/seqArchRplus.html#since)
  [![Bioc last commit](https://bioconductor.org/shields/lastcommit/devel/bioc/seqArchRplus.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/seqArchRplus/)
  [![Bioc dependencies](https://bioconductor.org/shields/dependencies/devel/seqArchRplus.svg)](https://bioconductor.org/packages/devel/bioc/html/seqArchRplus.html#since)
  <!-- badges: end -->

# seqArchRplus

## Downstream analyses of promoter sequence architectures and HTML report generation

seqArchRplus facilitates downstream analyses of promoter sequence 
architectures/clusters identified by seqArchR (or any other tool/method). 
With additional available information such as the TPM values and interquantile 
widths (IQWs) of the CAGE tag clusters, seqArchRplus can order the input 
promoter clusters by their shape (based on IQWs information), and write the 
cluster information as browser/IGV track files. 
Provided visualizations are of two kind: per sample/stage and per cluster 
visualizations. 
Those of the first kind include: plot panels for each sample showing per 
cluster shape, TPM and other score distributions, sequence logos, and peak 
annotations. 
The second include per cluster chromosome-wise and strand 
distributions, motif occurrence heatmaps and GO (Gene Ontology) term 
enrichments. 
Additionally, seqArchRplus can also generate HTML reports for easy viewing and 
comparison of promoter architectures between samples/stages (future).


## Installation 

The latest stable version of seqArchRplus can be installed from Bioconductor 
as shown below.

```
if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("seqArchRplus")
```

In case of any errors or help required, please consider looking up: 
https://github.com/snikumbh/seqArchRplus and file a new issue.
