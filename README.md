
<!-- README.md is generated from README.Rmd. Please edit that file -->

# evalGSVAsig

<!-- badges: start -->
<!-- badges: end -->

The goal of evalGSVAsig is to …

## Installation

You can install the development version of evalGSVAsig from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("lkroeh/evalGSVAsig")
```

## Example

This is a basic example which shows the format of the data:

``` r
library(evalGSVAsig)

#example gene list
sig <- c("GENE1", "GENE2", "GENE3")
signature_list <- c(list(sig))
names(signature_list) <- c("signature1")

#run function
output <- evalGSVAsig::GSVAsignatureRanking(eset, signature_list)

#view output
#print df of genes ordered by correlation to GSVA scores
output[[1]]

#show heatmap of ALL gene expression in relation to GSVA score
output[[2]]

#show heatmap of SIGNATURE gene expression in relation to GSVA score
output[[3]]

#get expression with GSVA scores saved in pData
output[[4]]
```

With sample data:

``` r
#with our sample data
data(signatures)
data(eset)

output <- evalGSVAsig::GSVAsignatureRanking(eset, signatures)
#> Warning: replacing previous import 'Biobase::combine' by 'dplyr::combine' when
#> loading 'evalGSVAsig'
#> Warning: replacing previous import 'dplyr::filter' by 'stats::filter' when
#> loading 'evalGSVAsig'
#> Warning: replacing previous import 'dplyr::lag' by 'stats::lag' when loading
#> 'evalGSVAsig'
#> Warning in .filterFeatures(expr, method): 13 genes with constant expression
#> values throuhgout the samples.
#> Warning in .filterFeatures(expr, method): Since argument method!="ssgsea", genes
#> with constant expression values are discarded.
#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero

#> Warning in stats::cor(Biobase::exprs(eset)[i, ], eset$signature_gsvascore, : the
#> standard deviation is zero
```

View tables:

``` r
#This table contains only signature genes
head(output[[1]])
#>     correlation   gene rank
#> 470   0.7223057 WFDC12    1
#> 75    0.7145882 ASPRV1    2
#> 481   0.7058520  LCE3E    3
#> 458   0.6744871   DSC1    4
#> 234   0.6293875   DSG1    5
#> 454   0.6065752   ARG1    6
#This table contains all genes; right now there is a bug so it doesn't work.
head(output[[2]])
#>      [,1] [,2] [,3]
#> [1,]   NA   NA   NA
#> [2,]   NA   NA   NA
#> [3,]   NA   NA   NA
#> [4,]   NA   NA   NA
#> [5,]   NA   NA   NA
#> [6,]   NA   NA   NA
```

View heatmap that plots all signature and non-signature genes:

<img src="man/figures/README-heatmaps-1.png" width="100%" />

View heatmap that plots only signature genes:

<img src="man/figures/README-heatmap2-1.png" width="100%" />
