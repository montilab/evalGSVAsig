#' This function generates an annotated heatmap that shows the contributions of genes to the given signature score
#'
#' @param eset an expression set
#' @param signature a list of signatures (at least 1)
#' @param sig_score an aggregate signature score
#' @param col_ann a ComplexHeatmap::heatmapAnnotation object with columns' (i.e., samples') annotation
#' @param method how to compute the aggregate score (only GSVA at the moment)
#' @param ... additional parameters to pass to the method
#'
#' @return a list of a dataframe and two heatmaps
#' @export
#' @import ComplexHeatmap stats Biobase dplyr grid
#'
#' @examples
omics_signature_heatmap <- function(
    eset,
    signature,
    sig_score = NULL,
    col_ann = NULL,
    method = c("GSVA"),
    ...
)
{
  ## input checks
  stopifnot( is(eset, "ExpressionSet") )
  stopifnot( is.null(sig_score) || length(sig_score)==ncol(eset) )
  stopifnot( isTRUE(all.equal(names(sig_score), sampleNames(eset))) )
  stopifnot( is.null(heatmap_ann) || isTRUE(all(rownames(heatmap_ann) %in% sampleNames(eset))))

  ## compute score if not provided
  if ( is.null(sig_score) ) {
    sig_score <- compute_signature_score( eset = eset, signature = signature, method = method, ... )
  }
  ## add signature score to eset metadata
  eset$sig_score <- sig_score

  ## make dataframe of correlation of each gene with signature
  fData(eset)$score_cor <-
    cor(eset$sig_score, t(exprs(eset)))
  fData(eset)$cor_rank <-
    rank(-fData(eset)$score_cor,  ties.method = "first")

  df <- data.frame(
    gene = featureNames(eset),
    correlation = cor(eset$sig_score, t(exprs(eset)))
  ) |>
    dplyr::arrange(desc(correlation)) |>
    dplyr::mutate(rank = seq(1,length(correlation)))

  ## extract only genes within signature
  df_flt <- df |>
    dplyr::filter(gene %in% signature[[1]])
  stopifnot( nrow(d_flt) > max(5,length(signature) * 0.25) )

  ## PLOTS

  ## 1) with all genes
  eset_srt <- eset[order(fData(eset)$rank),
                   order(eset$sig_score, decreasing = TRUE)]

  fData(eset_srt)$insig <-
    ifelse(featureNames(eset_srt) %in% signature[[1]], 'SignatureGene', 'Background')

  if ( is.null(col_ann)) {
    col_ann <- HeatmapAnnotation(sig_score = anno_barplot(eset_srt$sig_score))
  } else {
    col_ann <- HeatmapAnnotation(sig_score = anno_barplot(eset_srt$sig_score),
                                 Annotation = col_ann)
  }
  row_ann <- HeatmapAnnotation(
    genes = fData(eset_srt)$insig,
    col = list(genes = c("Background" = "brown", "SignatureGene" = "lightgreen")),
    show_annotation_name = FALSE, which = 'row'
  )
  full_heatmap <- ComplexHeatmap::Heatmap(
    matrix = t(scale(t(exprs(eset_srt)))),
    name = "expression",
    top_annotation = column_ha,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = "Genes") +
    row_ann +
    rowAnnotation(correlation = anno_barplot(fData(eset_srt)$sig_cor))

  ## 2) with only signature genes
  eset_flt <- eset_srt[featureNames(eset_srt) %in% signature[[1]],]
  stopifnot( nrow(eset_flt) > max(5, length(signature) * .25) )

  sig_heatmap <- ComplexHeatmap::Heatmap(
    matrix = t(scale(t(exprs(eset_flt)))),
    name = "expression",
    top_annotation = col_ann,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    row_title = "Signature Genes",
    row_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = grid::gpar(fontsize = 4),
    show_row_names = TRUE,
    row_names_side = "left"
  ) +
    rowAnnotation(correlation = anno_barplot(df3$correlation))

  return(list(
    heatmap_all_genes = full_heatmap,
    heatmap_sig_genes = sig_heatmap
  ))
}
