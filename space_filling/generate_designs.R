# Comparison of the Space-filling properties for MmLHD, MaxPro, MOFAT, SOFT (standard), SOFT (strict)
# p = 2,...,20
# l = 6

rm(list = ls())
# Resolve paths from the SOFT repository root or any subdirectory.
soft_root <- normalizePath(".")
while (!file.exists(file.path(soft_root, "functions.R"))) {
  if (identical(soft_root, dirname(soft_root))) {
    stop("Run this script from the SOFT repository or one of its subdirectories.")
  }
  soft_root <- dirname(soft_root)
}

setwd(file.path(soft_root, "space_filling", "designs"))
l = 6

# MmLHD
MmLHD_info <- lapply(2:20, function(i) {
  print(i)
  n <- l * (i + 1)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- SLHD::maximinSLHD(1, n, i)$StandDesign
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
MmLHD_D <- lapply(MmLHD_info, function(i) {
  lapply(i, `[[`, "design")
})
MmLHD_time <- lapply(MmLHD_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(MmLHD_D, MmLHD_time, file = "prerun MmLHD.RData")


# MaxPro
maxpro_info <- lapply(2:20, function(i) {
  print(i)
  n <- l * (i + 1)
  lapply(1:30, function(j) {
    timing <- system.time({
      InitialDesign <- MaxPro::MaxProLHD(n = n, p = i)$Design 
      D <- MaxPro::MaxPro(InitialDesign)$Design
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
maxpro_D <- lapply(maxpro_info, function(i) {
  lapply(i, `[[`, "design")
})
maxpro_time <- lapply(maxpro_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(maxpro_D, maxpro_time, file = "prerun maxpro.RData")


# MOFAT
mofat_info <- lapply(2:20, function(i) {
  print(i)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- MOFAT::mofat(i, l, "projection", init = "MmLHD")
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
mofat_D <- lapply(mofat_info, function(i) {
  lapply(i, `[[`, "design")
})
mofat_time <- lapply(mofat_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(mofat_D, mofat_time, file = "prerun mofat.RData")


# soft - standard
stand_info <- lapply(2:20, function(i) {
  print(i)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- MOFAT::soft(i, l)
      D <- D$design
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
stand_D <- lapply(stand_info, function(i) {
  lapply(i, `[[`, "design")
})
stand_time <- lapply(stand_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(stand_D, stand_time, file = "prerun stand.RData")


# soft - strict
strict_info <- lapply(2:20, function(i) {
  print(i)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- MOFAT::soft(i, l, structure = "strict")
      D <- D$design
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
strict_D <- lapply(strict_info, function(i) {
  lapply(i, `[[`, "design")
})
strict_time <- lapply(strict_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(strict_D, strict_time, file = "prerun strict.RData")


######################################################### augmentation
######################################################### 
l_aug <- 4

# Augmented MOFAT
aug_mofat_info <- lapply(2:20, function(i) {
  print(i)
  n_aug <- (l-l_aug) * (i+1)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- MOFAT::mofat(i, l_aug, "projection", init = "MmLHD")
      D <- MOFAT::a_soft(i, l_aug, n_aug = n_aug, ofat_design = D)
    })
    list(design = D$design, elapsed_time = unname(timing["elapsed"]))
  })
})
aug_mofat_D <- lapply(aug_mofat_info, function(i) {
  lapply(i, `[[`, "design")
})
aug_mofat_time <- lapply(aug_mofat_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(aug_mofat_D, aug_mofat_time, file = "prerun aug mofat.RData")


# soft - standard
aug_stand_info <- lapply(2:20, function(i) {
  print(i)
  n_aug <- (l-l_aug) * (i+1)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- MOFAT::a_soft(i,l_aug,n_aug = n_aug)
      D <- D$design
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
aug_stand_D <- lapply(aug_stand_info, function(i) {
  lapply(i, `[[`, "design")
})
aug_stand_time <- lapply(aug_stand_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(aug_stand_D, aug_stand_time, file = "prerun aug stand.RData")


# soft - strict
aug_strict_info <- lapply(2:20, function(i) {
  print(i)
  n_aug <- (l-l_aug) * (i+1)
  lapply(1:30, function(j) {
    timing <- system.time({
      D <- MOFAT::a_soft(i, l_aug, structure = "strict", n_aug = n_aug)
      D <- D$design
    })
    list(design = D, elapsed_time = unname(timing["elapsed"]))
  })
})
aug_strict_D <- lapply(aug_strict_info, function(i) {
  lapply(i, `[[`, "design")
})
aug_strict_time <- lapply(aug_strict_info, function(i) {
  unlist(lapply(i, `[[`, "elapsed_time"))
})
save(aug_strict_D, aug_strict_time, file = "prerun aug strict.RData")
