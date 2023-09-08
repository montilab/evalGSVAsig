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
    col_ha = NULL,
    method = c("GSVA"),
    ...
)
{
  ## input checks
  stopifnot( methods::is(eset, "SummarizedExperiment") || methods::is(eset, "ExpressionSet") )
  stopifnot( is.null(sig_score) || length(sig_score)==ncol(eset) )
  stopifnot( isTRUE(all.equal(names(sig_score), Biobase::sampleNames(eset))) )
  stopifnot( is.null(col_ha) || isTRUE(all(rownames(col_ha) %in% Biobase::sampleNames(eset))))

  if ( methods::is(eset, "SummarizedExperiment") ) {
    eset <- Biobase::ExpressionSet(
      assayData = SummarizedExperiment::assay(eset),
      phenoData = Biobase::AnnotatedDataFrame(colData(eset)),
      featureData = Biobase::AnnotatedDataFrame(rowData(eset))
    )
  }
  ## compute score if not provided
  if ( is.null(sig_score) ) {
    sig_score <- omics_signature_score( eset = eset, signature = signature, method = method, ... )
  }
  ## add signature score to eset metadata
  eset$sig_score <- sig_score

  ## add correlation (and rank) of each gene with signature score to fData
  Biobase::fData(eset)$score_cor <-
    stats::cor(eset$sig_score, t(exprs(eset)))[1,]

  ## PLOTS

  ## 1) with all genes
  eset_srt <- eset[
    order(Biobase::fData(eset)$score_cor, decreasing = TRUE), # high to low correlation
    order(eset$sig_score, decreasing = TRUE)         # high to low sig score
  ]
  Biobase::fData(eset_srt)$insig <- factor(
    ifelse(Biobase::featureNames(eset_srt) %in% signature[[1]], 'signature', 'background')
  )
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
    correlation = ComplexHeatmap::anno_barplot(Biobase::fData(eset_srt)$score_cor),
    col = list(genes = c("background" = "brown", "signature" = "lightgreen")),
    show_annotation_name = FALSE
  )
  full_heatmap <- ComplexHeatmap::Heatmap(
    matrix = t(scale(t(exprs(eset_srt)))),
    name = "expression",
    top_annotation = col_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = "Genes",
    show_row_names = FALSE,
    column_names_gp = grid::gpar(fontsize = 8)) +
    row_ha

  ## 2) with only signature genes
  eset_flt <- eset_srt[Biobase::featureNames(eset_srt) %in% signature[[1]],]
  stopifnot( nrow(eset_flt) > max(5, length(signature) * .25) )

  sig_heatmap <- ComplexHeatmap::Heatmap(
    matrix = t(scale(t(exprs(eset_flt)))),
    name = "expression",
    top_annotation = col_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = "Signature Genes",
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 8),
    show_row_names = TRUE,
    row_names_side = "left") +
    ComplexHeatmap::rowAnnotation(
      correlation = ComplexHeatmap::anno_barplot(Biobase::fData(eset_flt)$score_cor))

  return(list(
    score_cor = Biobase::fData(eset_srt) |> dplyr::select(score_cor),
    sig_score = Biobase::pData(eset_srt) |> dplyr::select(sig_score),
    heatmap_all_genes = full_heatmap,
    heatmap_sig_genes = sig_heatmap
  ))
}
