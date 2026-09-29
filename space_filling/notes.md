# Notes: Regenerate the space-filling figure

- Requested script: `code/space_filling/sf_prop_vs.R`.
- Existing output: `code/space_filling/visualization/combined_sf.png`.
- The main script computes Q and entropy before plotting, so it must not be sourced in full.
- Entropy candidates: `designs/ent.RData` (modified September 26, 2026) and `visualization/ent.RData` (modified August 28, 2026).
- Root `task_plan.md` and `notes.md` document an earlier Metal optimization task and will remain untouched.
- All five updated entropy matrices and all five saved Q matrices are 30 by 19, positive, and finite.
- The five `designs/prerun <method>.RData` files contain runtime lists that reshape to 30 by 19 matrices.
- Reuse the original code starting at `############################# computational time`; this only reshapes saved runtimes and renders the figure.
- Existing export settings: three panels, a shared legend, 9 by 3 inches at 600 dpi, followed by whitespace trimming.

## Result and verification
- Regenerated `visualization/combined_sf.png` (5323 by 1759 pixels, 477708 bytes).
- Each panel has 95 summary rows (5 methods by 19 factor counts); all quartiles are finite and ordered.
- Entropy panel medians match the log ratios computed from the updated `designs/ent.RData`.
- Visually inspected the regenerated image: all three panels, axes, titles, and shared legend render correctly.
- Output checksum matches the inspected render and differs from the previous figure.
- Checksums unchanged for the main script, both entropy files, saved Q file, and all five saved design/runtime files (9 files).
- No design generation, Q evaluation, or entropy computation was rerun.
- Previous image and verification record retained in `/tmp/sf_saved_figure_16b25106c69bc/` for this session.
- Changed existing file: `visualization/combined_sf.png` only. Added task tracking: this `notes.md` and adjacent `task_plan.md`.

## Repeat regeneration
- Reloaded `designs/ent.RData`, updated September 26, 2026 at 00:49:12 (16247 bytes), with saved Q and runtime data.
- Regenerated `visualization/combined_sf.png`: 5323 by 1759 pixels, 476420 bytes; original layout preserved.
- Verified finite, ordered quartiles for all 95 rows in each panel and checked entropy medians against current saved inputs.
- Visually inspected all three panels; installed output matches inspected render and differs from previous figure.
- All 9 source/data file checksums remain unchanged. No design generation or criterion computation was rerun.
- Prior image and verification record: `/tmp/sf_saved_repeat_16ca96612cf55/`.
- Updated project files: the figure, this task note, and `task_plan.md`.

## Why the revised criterion separates the design families
- Question: explain the large MmLHD/MaxPro scores and small MOFAT/SOFT scores for `1 / sum_{i<j} prod_k beta/(beta + delta_k^2)`.
- Inspected `generate_designs.R`; all methods use 126 points at p=20, with coordinates on the unit-cube scale. All 30 replicates per method were checked.
- With beta=1/144, the saved p=20 median scores are: MmLHD 2.6879630e13; MaxPro 5.5497284e14; MOFAT 0.24534050; SOFT-Stand 0.27214703; SOFT-Strict 0.32830976.
- Recomputed pairwise sums match the saved scores to at most 2.31e-14 relative error, excluding stale data or numerical underflow as the cause of this gap in the checked runs.
- Every p=20 MOFAT/SOFT replicate contains exactly 120 pairs differing in one coordinate. Their shares of the total sum are MOFAT 79.569892%, SOFT-Stand 80.941723%, SOFT-Strict 97.645590%.
- MOFAT and SOFT-Stand additionally have 1140 pairs differing in exactly two coordinates; these together with the one-coordinate pairs account for essentially the entire sum. SOFT-Strict has 114 two-coordinate pairs; one- and two-coordinate pairs account for 99.944721%.
- All pairs in the saved MmLHD and MaxPro designs differ in all 20 coordinates. Median kernel sums are approximately 3.72e-14 and 1.80e-15, respectively.
- MOFAT one-coordinate step is 0.5, giving a pair contribution a=1/37. Six stars contribute `6 * (20*a + choose(20,2)*a^2) = 4.0759679`, predicting score 0.2453405. Contributions between stars are negligible at p=20.
- For MOFAT, the one-coordinate pairs alone imply an upper bound on its score of 37/120=0.3083333.
- Thus the high-dimensional product makes unrestricted pairs nearly uncorrelated while the OFAT construction retains fixed within-structure contributions. This is a geometric consequence of the criterion, not evidence of a comparable ratio in predictive accuracy or information.
- For fixed p and beta the minimized sum is a regularized version of the MaxPro pairwise inverse-product objective (constant beta^p omitted). It is a pairwise similarity sum, not the GP log-determinant entropy criterion.
- Source: Joseph, Gul and Ba (2015), equations (5) and Section 5, https://www.asc.ohio-state.edu/statistics/comp_exp/jour.club/max_projection_designs_roshan2015.pdf .
- This explanation updated only task tracking notes. Scripts, saved datasets, and figure are unchanged.
