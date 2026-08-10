# Issue: weighted .kstest() silently returns NaN on all-zero weights

## Summary
`.kstest()`'s weighted branch (used when `signature_projection_contributors(..., gsea = TRUE)`) could silently produce `NaN` for `score`/`leading_edge`/the enrichment plot, with no error or warning, whenever all signature genes had zero weight.

## Root cause
```r
Phit <- rep(0, n.x)
Phit[y] <- weights
Phit <- cumsum(Phit)
Phit <- Phit / Phit[n.x]  #this making error?
```
The author's own inline comment ("this making error?") flagged the suspicion. If `sum(weights) == 0` (e.g. all signature genes have zero correlation to the score, or a `weights.pwr` that zeroes them out), `Phit[n.x]` is `0`, and the division silently produces a vector of `NaN` that propagates into `z`, `score`, and `leading_edge` without any signal to the caller.

## Fix implemented
Added an explicit zero-weight guard immediately after computing `weights`, returning the function's existing degenerate-input result (`score = 0, pval = 1`, empty plot) with an explicit `warning()` instead of silently continuing into `NaN`:

```r
weights <- abs(weights[y])^weights.pwr
if ( sum(weights) == 0 ) {
  warning("all weights are zero; cannot compute weighted enrichment score")
  return(err)
}
```

## Files changed
- `R/ks_enrichment.R`

## Verification
- Reproduced the silent `NaN` with an all-zero weights vector before the fix.
- Confirmed the fix returns the degenerate result with a caught, informative warning instead.
- Full `testthat` suite (50/50) passing.

## Notes
Commit: `5c0b013` (same commit as the extract_eigengenes fix).
