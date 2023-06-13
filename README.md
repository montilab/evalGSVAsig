
# evalGSVAsig

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of evalGSVAsig is to ...

## Installation

You can install the development version of evalGSVAsig like so:

``` r
devtools:
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(evalGSVAsig)
#example gene list
sig <- c("GENE1", "GENE2", "GENE3")
signature_list <- c(list(sig))
names(signature_list) <- c("signature1")

output <- GSVAsignatureRanking(eset, signature_list)
#print df of genes ordered by correlation to GSVA scores
output[[1]]

#show heatmap of ALL gene expression in relation to GSVA score
output[[2]]

#show heatmap of SIGNATURE gene expression in relation to GSVA score
output[[3]]

#get expression with GSVA scores saved in pData
output[[4]]
```

