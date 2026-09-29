# Task Plan: Regenerate the space-filling figure from saved data

## Goal
Regenerate the existing combined figure using the updated saved entropy data, without rerunning designs or criterion calculations.

## Actions
- [x] Inspect the plotting block and confirm the saved inputs, including which entropy file is current.
- [x] Load saved Q, entropy, and timing data and execute only the existing plotting code.
- [x] Inspect the regenerated image and verify that source scripts and saved data were not modified.

## Decisions
- Keep the root planning files for the previous task intact; track this task here.
- Preserve the existing figure style and output location.
- Use `designs/ent.RData`, the September 26 update; the visualization copy is still from August 28.
- Execute only the timing reshape and visualization block in a temporary R session; leave source scripts unchanged.

## Errors Encountered
- None.

## Status
Repeat regeneration complete. Reloaded `designs/ent.RData` (updated September 26 at 00:49:12), regenerated `visualization/combined_sf.png`, inspected the result, and verified unchanged source/data checksums.

## Follow-up: Explain the criterion gap
- [x] Inspect saved design sizes, coordinate scales, and entropy medians.
- [x] Decompose pairwise kernel sums by the number of differing coordinates and verify representative saved values.
- [x] Explain the geometric mechanism, metric direction, and limits using the measured values.

## Follow-up execution note
- `designs/designs.R` does not exist; locate the actual generation script before inspecting design construction.
- Found and inspected `generate_designs.R`. Evaluated saved designs only; no design generation or source/data/figure changes.
- OSTI copy of the MaxPro paper timed out; the Ohio State-hosted original paper was accessible.
- Follow-up complete: checked all 30 replicates for all 5 methods at p=20; saved criterion values agree within 2.4e-14 relative error.
