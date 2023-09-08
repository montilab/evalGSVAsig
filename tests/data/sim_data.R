library(Biobase)
library(ggplot2)
nsample <- 20
ngene <- 100
nsig <- 10

set.seed(123)

sig_score <- rnorm(nsample)
rnd_pheno <- sample(c("head","tail"), nsample, replace = TRUE)
dat_matrix <- rbind(
  t(replicate(nsig, rnorm(nsample,mean = sig_score))),
  matrix(rnorm(n = nsample * (ngene-nsig)), ncol = nsample)
)
rownames(dat_matrix) <- sprintf("gene_%03d",seq(1,ngene))
names(sig_score) <- colnames(dat_matrix) <- sprintf("sample_%02d", seq(1,nsample))

eset <- Biobase::ExpressionSet(
  assayData = dat_matrix,
  phenoData = AnnotatedDataFrame(data.frame(sig_score = sig_score, pheno = rnd_pheno))
)
signature <- featureNames(eset)[c(seq(1,nsig), sample(seq(nsig+1,ngene), size = 3))]

## checking
if (FALSE) {
  df <- data.frame(cor = cor(sig_score,t(dat_matrix))[1,],
                   sig = rownames(dat_matrix) %in% signature)
  boxplot(df$cor ~ df$sig, xlab = "in signature", ylab = "cor( gene, score )")
}
saveRDS(eset, file.path(".", "tests", "data", "sim_data_eset.rds"))
saveRDS(signature, file.path(".", "tests", "data", "sim_data_signature.rds"))

## test the functions on the simulated dataset
if ( FALSE )
{
  devtools::load_all(".")

  eset <- readRDS(file.path(".", "tests", "data", "sim_data_eset.rds"))
  signature <- readRDS(file.path(".", "tests", "data", "sim_data_signature.rds"))
  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = pData(eset) |> dplyr::select(pheno),
    col = list(pheno = c(head = "pink", tail = "green"))
  )
  output <- omics_signature_heatmap(
    eset = eset,
    signature = list(signature),
    sig_score = pData(eset) |> dplyr::select(sig_score) |> tibble::rownames_to_column() |> tibble::deframe(),
    col_ha = col_ha
  )
<<<<<<< HEAD
  print(head(output$score_cor))
  print(output$sig_score)
=======
  print(output$sig_score)
  print(output$score_cor[1:20,])
  print(output$heatmap_all_genes)
  print(output$heatmap_sig_genes)

  ks_out <- .kstest(
    n.x = nrow(output$score_cor),
    y = rank(-output$score_cor$score_cor)[output$score_cor$insig == "signature"],
    plotting = TRUE
  )
  ## show w/ package data
  data(eset)
  data(signatures)
  output <- omics_signature_heatmap(
    eset = eset,
    signature = signatures,
    col_ha = ComplexHeatmap::columnAnnotation(
      df = pData(eset) |> dplyr::select(hpv_status),
      col = list(hpv_status = c(neg = "gray", pos = "pink")))
  )
  print(output$sig_score[1:20,, drop = FALSE])
  print(output$score_cor[1:20,])
>>>>>>> 2600d3946e751a76bea5c3d9f2823eddbb6f6ab0
  print(output$heatmap_all_genes)
  print(output$heatmap_sig_genes)
}
