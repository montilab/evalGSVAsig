library(Biobase)
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
signature <- featureNames(eset)[seq(1,nsig)]

## checking
if (FALSE) {
  df <- data.frame(cor = cor(sig_score,t(dat_matrix))[1,],
                   sig = rownames(dat_matrix) %in% signature)
  boxplot(df$cor ~ df$sig, xlab = "in signature", ylab = "cor( gene, score )")
}
saveRDS(eset, "sim_data_eset.rds")
saveRDS(eset, "sim_data_signature.rds")

## test the functions on the simulated dataset
if ( FALSE )
{
  devtools::load_all(".")

  col_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = pData(eset) |> dplyr::select(pheno),
    col = list(pheno = c(head = "pink", tail = "green"))
  )
  heatmaps <- omics_signature_heatmap(
    eset = eset,
    signature = list(signature),
    sig_score = pData(eset) |> dplyr::select(sig_score) |> tibble::rownames_to_column() |> tibble::deframe(),
    col_ha = col_ha
  )
  print(heatmaps$heatmap_all_genes)
  print(heatmaps$heatmap_sig_genes)
}
