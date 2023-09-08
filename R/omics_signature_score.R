#' This function computes an aggregate score for a set of genes
#'
#' @param eset an expression set
#' @param signature a list of signatures (at least 1)
#' @param method how to compute the aggregate score (only GSVA at the moment)
#' @param ... additional parameters to pass to the method
#'
#' @return a named numeric vector of scores (w/ the names corresponding to colnames(eset))
#' @importFrom GSVA gsva
#' @importFrom Biobase ExpressionSet exprs featureNames
#' @importFrom tibble deframe
#' @importFrom methods is
#'
#' @export
omics_signature_score <- function(
    eset,
    signature,
    method = c("GSVA"),
    ...)
{
  ## input checks
  stopifnot( methods::is(eset,"ExpressionSet") )
  stopifnot( isTRUE(all(signature %in% Biobase::featureNames(eset))) )
  method <- match.arg(method)

  sig_score <- {
    if (method == "GSVA") {
      tmp <- GSVA::gsva(eset, signature, verbose = FALSE, ...)
      stopifnot(nrow(tmp)==1) # working w/ single signature only at the moment
      t(exprs(tmp)) |>
        data.frame(check.names = FALSE) |>
        tibble::deframe()
    } else {
      stop( "unrecognized method:", method)
    }
  }
  return(sig_score)
}
