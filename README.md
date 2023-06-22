
# evalGSVAsig

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

The goal of evalGSVAsig is to understand which genes in a signature contribute most to a GSVA score.

## Installation

You can install the development version of evalGSVAsig like so:

``` r
devtools::install_github("lkroeh/evalGSVAsig")
```

## Example

This is a basic example which shows you how to run the function:

``` r
library(evalGSVAsig)

#example gene list
sig <- c("GENE1", "GENE2", "GENE3")
signature_list <- c(list(sig))
names(signature_list) <- c("signature1")

#run function
output <- GSVAsignatureRanking(eset, signature_list)

#view output
#print df of genes ordered by correlation to GSVA scores
output[[1]]

#show heatmap of ALL gene expression in relation to GSVA score
output[[2]]

#show heatmap of SIGNATURE gene expression in relation to GSVA score
output[[3]]

#get expression with GSVA scores saved in pData
output[[4]]


#with our sample data
signaturelist <- data(signatures.rda)
eset <- data(eset.rda)

output <- GSVAsignatureRanking(eset, signature_list)

head(output[[1]])
head(output[[2]])
output[[3]]
output[[4]]
```

