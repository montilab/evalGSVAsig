#' This function generates an annotated heatmap that shows the contributions of genes to the given signature score
#'
#' @param eset an expression set
#' @param signature a list of signatures (at least 1)
#' @param sig_score an aggregate signature score
#' @param col_ha a ComplexHeatmap::heatmapAnnotation object with columns' (i.e., samples') annotation
#' @param method how to compute the aggregate score (only GSVA at the moment)
#' @param ... additional parameters to pass to the method
#'
#' @return a list of a dataframe and two heatmaps
#'
#' @import Biobase
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation rowAnnotation anno_barplot
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom stats cor
#' @importFrom methods is
#' @importFrom grid gpar
#'
#' @export
omics_signature_heatmap <- function(
    eset,
    signature,
    sig_score = NULL,
    cor_method = c("pearson", "spearman", "kendall"),
    col_ha = NULL,
    min_sigsize = 3,
    method = c("GSVA"),
    name = "expression",
    gsea = FALSE,
    ...
) {
  ## input checks
  cor_method = match.arg(cor_method)
  stopifnot( methods::is(eset, "SummarizedExperiment") || methods::is(eset, "ExpressionSet") )
  stopifnot( is.null(sig_score) || length(sig_score)==ncol(eset) )
  stopifnot( is.null(sig_score) || isTRUE(all.equal(names(sig_score), sampleNames(eset))) )
  stopifnot( is.null(col_ha) || isTRUE(all(rownames(col_ha) %in% sampleNames(eset))))

  if ( methods::is(eset, "SummarizedExperiment") ) {
    eset <- Biobase::ExpressionSet(
      assayData = SummarizedExperiment::assay(eset),
      phenoData = Biobase::AnnotatedDataFrame(colData(eset)),
      featureData = Biobase::AnnotatedDataFrame(rowData(eset))
    )
  }
  ## compute score if not provided
  if ( is.null(sig_score) ) {
    sig_score <- omics_signature_score( eset = eset, signature = signature, method = method)
  }
  ## add signature score to eset metadata
  eset$sig_score <- sig_score
  stopifnot( all(!is.na(eset$sig_score)) )

  ## add correlation (and p-value) of each gene with signature score to fData
  COR <- psych::corr.test(eset$sig_score, t(exprs(eset)), method = cor_method)
  stopifnot(nrow(COR$r)==1)
  stopifnot(nrow(COR$p)==1)
  fData(eset)$score_cor <- drop(COR$r)
  fData(eset)$pval_cor <- drop(COR$p)
  fData(eset)$insig <- factor(
    ifelse(featureNames(eset) %in% signature[[1]], 'signature', 'background'),
    levels = c("signature","background")
  )
  ## PLOTS

  ## 1) with all genes
  eset_srt <- eset[
    order(Biobase::fData(eset)$score_cor, decreasing = TRUE), # high to low correlation
    order(eset$sig_score, decreasing = TRUE)                  # high to low sig score
  ]
  ks_out <-   .kstest(
    n.x = nrow(eset),
    y = rank(-fData(eset)$score_cor)[fData(eset)$insig == "signature"],
    weights = if (gsea) fData(eset_srt)$score_cor,
    plotting = TRUE
  )
  ## from idx to names
  fData(eset_srt)$leading_edge <- NA
  if (is.null(ks_out$leading_edge)) {
    ks_out$hits <- NA
  } else if (!is.null(ks_out$leading_edge) && ks_out$leading_edge == 0) {
    ks_out$hits <- NA
  } else {
    ks_out$hits <- Biobase::fData(eset_srt) |>
      dplyr::slice_head(n = ks_out$leading_edge) |>
      dplyr::filter(insig == "signature") |>
      tibble::rownames_to_column(var = "featureID") |>
      dplyr::pull(featureID)
    ## add information about hits in leading edge
    fData(eset_srt)$leading_edge <-
      factor(ifelse(featureNames(eset_srt) %in% ks_out$hits, "yes", "no"),
             levels = c("yes", "no"))
  }
  if ( is.null(col_ha)) {
    ## the only column annotation will be the sig_score barplot
    col_ha <- ComplexHeatmap::HeatmapAnnotation(
      sig_score = ComplexHeatmap::anno_barplot(eset_srt$sig_score))
  } else {
    ## augment input column annotation with sig_score barplot
    col_ha <- c(
      ComplexHeatmap::HeatmapAnnotation(
        sig_score = ComplexHeatmap::anno_barplot(eset_srt$sig_score)),
      ## next command is not a robust solution, but couldn't find a better way
      col_ha[match(Biobase::sampleNames(eset_srt),Biobase::sampleNames(eset)),])
  }
  row_ha <- ComplexHeatmap::rowAnnotation(
    genes = Biobase::fData(eset_srt)$insig,
    leadedge = ifelse(featureNames(eset_srt) %in% ks_out$hits, "yes", "no"),
    correlation = ComplexHeatmap::anno_barplot(Biobase::fData(eset_srt)$score_cor),
    col = list(genes = c("background" = "brown", "signature" = "lightgreen"),
               leadedge = c(yes = "black", no = "white")),
    show_annotation_name = FALSE
  )
  full_heatmap <- suppressMessages(ComplexHeatmap::Heatmap(
    matrix = t(scale(t(exprs(eset_srt)))),
    top_annotation = col_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = "Genes",
    show_row_names = FALSE,
    column_names_gp = grid::gpar(fontsize = 8),
    ... )) +
    row_ha

  ## 2) with only signature genes
  eset_flt <- eset_srt[Biobase::featureNames(eset_srt) %in% signature[[1]],]
  stopifnot( nrow(eset_flt) > max(min_sigsize, length(signature) * .25) )

  sig_heatmap <- suppressMessages(ComplexHeatmap::Heatmap(
    matrix = t(scale(t(exprs(eset_flt)))),
    top_annotation = col_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = "Signature Genes",
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    show_row_names = TRUE,
    row_names_side = "left",
    ... )) +
    ComplexHeatmap::rowAnnotation(
      correlation = ComplexHeatmap::anno_barplot(Biobase::fData(eset_flt)$score_cor))

  return(list(
    score_cor = Biobase::fData(eset_srt) |>
      dplyr::select(score_cor, pval_cor, insig, leading_edge),
    sig_score = Biobase::pData(eset_srt) |>
      dplyr::select(sig_score),
    heatmap_all_genes = full_heatmap,
    heatmap_sig_genes = sig_heatmap,
    ks = ks_out
  ))
}
