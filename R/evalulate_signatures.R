#' This function evaluates the contributions of genes in GSVA scores
#'
#' @param eset an expression set
#' @param signature a list of signatures (at least 1)
#' @param metacol name of metadata column in eset to annotate top of heatmap
#'
#' @return a list of a dataframe and two heatmaps
#' @export
#' @import GSVA ComplexHeatmap stats Biobase dplyr grid
#'
#' @examples
GSVAsignatureRanking <- function(
    eset,
    signature,
    metacol
)   {
  #set up empty list in which to store things function will return
  returnobject <- list()

  #run gsva with eset and signature
  esetgsva <- GSVA::gsva(eset, signature, verbose = FALSE)

  #add gsva results to eset metadata
  eset$signature_gsvascore <- t(Biobase::exprs(esetgsva[1,]))

  #make dataframe of correlation of each gene with signature
  df <- data.frame()
  count <- 1
  for (i in 1:nrow(Biobase::exprs(eset))) {
    a <- stats::cor(Biobase::exprs(eset)[i,], eset$signature_gsvascore, method = "pearson")
    df[i,1] <- a[1,1]
    df[i,2] <- rownames(eset)[[i]]
    count = count + 1
  }
  names(df) <- c("correlation", "gene")

  #sort df
  df1 <- df[order(-df$correlation),]
  #add rank numbers
  df1$rank <- seq.int(nrow(df1))
  returnobject[[1]] <- df1

  #extract only genes within signature
  df_filtered <- df1 %>%
    na.omit(filter(df1$gene %in% signature[[1]]))

  #save df of signature genes, their correlation with GSVA score and rankings
  returnobject[[2]] <- df_filtered

  #takethisout
  #test <- ksGenescore(nrow(df1), df_filtered$rank,do.pval = FALSE, do.plot= T)

  #plot

  #with all genes
  hmexprs <- Biobase::exprs(eset)
  rownames(df1) <- df1$gene

  #order columns
  sampledf <- as.data.frame(eset$signature_gsvascore)
  sampledf$samples <- rownames(sampledf)
  if (is.na(metacol)) {
    sampledf <- sampledf[order(-sampledf[,1]),]
    names(sampledf) <- c("score", "samples")
  } else {
    sampledf$ano <- eset[[metacol]]
    sampledf <- sampledf[order(-sampledf[,1]),]
    names(sampledf) <- c("score", "samples", "ano")
  }

  #make matrix
  test_ordered <- hmexprs[df1$gene, sampledf$samples]
  test_ordered2 <- test_ordered
  colnames(test_ordered2) = NULL
  mat_scaled = t(scale(t(test_ordered2)))
  df3 <- df1[order(-df1$correlation),]
  df3$insig <-  with(df3, ifelse(df3$gene %in% signature[[1]], 'SignatureGene', 'Background'))

  #add top annotation
  if (is.na(metacol)) {
    column_ha = HeatmapAnnotation(GsvaScore = anno_barplot(sampledf$score))
  } else {
    column_ha = HeatmapAnnotation(GsvaScore = anno_barplot(sampledf$score), Annotation = sampledf$ano)
  }

  row_ha = HeatmapAnnotation(genes = df3$insig,
                         col = list(genes = c("Background" = "brown", "SignatureGene" = "lightgreen")),
                         show_annotation_name = F, which = 'row')

  ht_list = ComplexHeatmap::Heatmap(mat_scaled, name = "expression",
                                    top_annotation = column_ha,
                                    show_column_dend = FALSE,cluster_rows = FALSE, cluster_columns = FALSE,
                                    cluster_column_slices = FALSE, show_row_dend = FALSE, row_title = "Genes") +
    row_ha +
    rowAnnotation(correlation = anno_barplot(df3$correlation))


  #store heatmap
  returnobject[[3]] <- ht_list

  #plot

  ##with only signature genes
  hmexprs <- Biobase::exprs(eset)[intersect(rownames(Biobase::exprs(eset)), (signature[[1]])),]
  rownames(df1) <- df1$gene

  #order rows
  df2 <- df1[intersect(signature[[1]], (rownames(eset))),]
  #sort
  df2 <- df2[order(df2$rank),]

  #order columns
  sampledf <- as.data.frame(eset$signature_gsvascore)
  sampledf$samples <- rownames(sampledf)
  if (is.na(metacol)) {
    sampledf <- sampledf[order(-sampledf[,1]),]
    names(sampledf) <- c("score", "samples")
  } else {
    sampledf$ano <- eset[[metacol]]
    sampledf <- sampledf[order(-sampledf[,1]),]
    names(sampledf) <- c("score", "samples", "ano")
  }

  #make matrix
  test_ordered <- hmexprs[df2$gene, sampledf$samples]
  test_ordered2 <- test_ordered
  colnames(test_ordered2) = NULL
  mat_scaled = t(scale(t(test_ordered2)))
  df3 <- df2[order(-df2$correlation),]
  df3$insig <-  with(df3, ifelse(df3$gene %in% signature[[1]], 'SignatureGene', 'BackgroundGene'))

  if (is.na(metacol)) {
    column_ha = HeatmapAnnotation(GsvaScore = anno_barplot(sampledf$score))
  } else {
    column_ha = HeatmapAnnotation(GsvaScore = anno_barplot(sampledf$score), Annotation = sampledf$ano)
  }

  ht_list1 = ComplexHeatmap::Heatmap(mat_scaled, name = "expression",
                                     top_annotation = column_ha,
                                     show_column_dend = FALSE,cluster_rows = FALSE,
                                     cluster_columns = FALSE, cluster_column_slices = FALSE,
                                     show_row_dend = FALSE, row_title = "Signature Genes",
                                     row_title_gp=grid::gpar(fontsize=12,fontface="bold"),
                                     row_names_gp = grid::gpar(fontsize = 4),
                                     show_row_names=TRUE,
                                     row_names_side="left") +
    rowAnnotation(correlation = anno_barplot(df3$correlation))

  #store heatmap
  returnobject[[4]] <- ht_list1

  #store eset with gsva scores
  returnobject[[5]] <- eset

  #label items in returnobject
  names(returnobject) <- c("ordered_dfall", "ordered_dffiltered", "heatmap_AllGenes", "heatmap_SignatureGenes", "eset_wGSVA")

  #return everything
  return(returnobject)
}
