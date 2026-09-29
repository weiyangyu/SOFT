# SOFT

R code for comparing Space-Filling One-Factor-at-a-Time (SOFT) designs and their augmented versions with MOFAT, MaxPro, and maximin Latin hypercube designs.

| Code | Purpose |
| --- | --- |
| `functions.R` | Benchmark functions and inert-factor effects. |
| `screening/` | Factor screening on Morris, g-function, and Levy examples. |
| `modeling/` | Gaussian and MIM surrogate modeling on eight benchmark functions. |
| `space_filling/` | Space-filling quality and design-generation runtime comparisons. |
| `trade_off/conclusion.R` | Summary figure comparing space-fillingness and screening capability. |

`modeling/lib.R` provides shared modeling routines and compiles `modeling/MIM_GP.cpp` through Rcpp. These helpers are loaded automatically by the experiment scripts.

## Setup

Tested with R 4.5.2 on macOS. A C++ compiler is required; installing `magick` from source also requires ImageMagick.

From the repository root, run in R:

```r
install.packages(c(
  "dplyr", "future", "future.apply", "ggplot2", "gridExtra", "gtools",
  "magick", "MaxPro", "nloptr", "png", "progressr", "Rcpp",
  "RcppArmadillo", "RcppHungarian", "reshape2", "rkriging", "scales",
  "sensitivity", "SFDesign", "SLHD", "spacefillr", "tidyr", "withr"
), repos = "https://cloud.r-project.org")
install.packages("MOFAT_2.0.tar.gz", repos = NULL, type = "source")
```

## Run experiments

Run the following in a terminal from the repository root. Each script starts in a fresh R session; paths are resolved relative to the repository.

```sh
for script in \
  screening/generate_designs.R \
  screening/screening_vs.R \
  modeling/generate_sfd.R \
  modeling/generate_ofat.R \
  modeling/generate_aug_ofat.R \
  modeling/ofat_mim.R \
  modeling/ofat_gauss.R \
  modeling/sfd_mim.R \
  modeling/sfd_gauss.R \
  modeling/aug_ofat_mim.R \
  modeling/aug_ofat_gauss.R \
  modeling/mim_mse_vs.R \
  modeling/aug_mim_mse_vs.R \
  modeling/aug_gauss_mse_vs.R \
  space_filling/generate_designs.R \
  space_filling/sf_prop_vs.R \
  space_filling/aug_sf_prop_vs.R \
  trade_off/conclusion.R
do
  Rscript --vanilla "$script" || break
done
```

Experiments use 30 design replicates. Modeling uses parallel workers, leaving two available CPU cores free when possible. Run generators before analyses, and analyses before figures; the trade-off figure requires both space-filling analyses.

## Outputs

Rerunning overwrites generated outputs in the existing folders below. A complete run produces 122 `.RData` files and 32 `.png` figures.

| Folder | Contents |
| --- | --- |
| `screening/designs/`, `modeling/designs/`, `space_filling/designs/` | Generated design collections. |
| `screening/results/` | Screening indices and figures. |
| `modeling/mim/`, `modeling/gauss/` | Surrogate-model results. |
| `modeling/visual_mim/`, `modeling/visual_gauss/` | Modeling comparison figures. |
| `space_filling/visualization/` | Space-filling summaries and figures. |
| `trade_off/` | `Q_vs_perturbation.png`. |

In modeling result filenames, `_0` denotes the original benchmark and `_1` its 20-factor extension with inert-factor effects; `aug_` denotes augmented designs.
