# Issue: README lacked section-level context; gene correlation label reverted

## Summary
Two related documentation-facing changes made in the same pass.

## Requested behavior
1. `README.Rmd` moved straight from section headings into code chunks with no explanation of what each step/plot showed. Add brief explanatory text ahead of each major section.
2. `spc_heatmap_all()`'s left annotation had been renamed from `correlation` to `gene correlation` earlier in the day (see the now-closed issue #19 and its fix in commit `0532cab`, itself a re-fix of an earlier regression). The maintainer decided, after seeing it rendered, that they preferred the original `correlation` label after all.

## Fix implemented
1. Added one to two sentences ahead of each of: the toy-data intro, the toy tables, the toy full-genes heatmap, the toy signature-only heatmap, the toy KS plot, the `omics_signature_score()` methods note, the precomputed-score example, the full-dataset intro, and each of that section's tables/heatmaps/KS plot. Re-rendered `README.Rmd` -> `README.md` (via `rmarkdown::render(..., output_format = "github_document")`) and regenerated `man/figures/README-*.png` to match.
2. Reverted `spc_heatmap_all()`'s left annotation name from `` `gene correlation` `` back to `correlation` in `R/spc_heatmaps.R`. Re-rendered the README again; the affected figures came back byte-identical to their pre-`gene correlation` state, confirming deterministic rendering and no other drift.

## Files changed
- `README.Rmd`
- `README.md` (generated)
- `man/figures/README-*.png` (generated, only the ones containing the left annotation changed)
- `R/spc_heatmaps.R`

## Verification
- Full `testthat` suite (50/50) passing.
- Visually confirmed the rendered annotation label and the new explanatory text in the regenerated `README.md`/figures.

## Notes
Commit: `90ddd51` (bundles both the README text and the label revert, per explicit request not to split this one further).

**Decision record:** issue #19 ("left annotation label should display 'gene correlation'") was fixed and closed earlier the same day, but the maintainer subsequently decided to keep `correlation` instead. Issue #19 was deliberately **left closed** — the revert is a maintainer preference, not evidence the original fix was wrong. If this comes up again, check this file and `2026-07-14-summarizedexperiment-coercion-fix.md`'s sibling docs for the full history before re-"fixing" it.
