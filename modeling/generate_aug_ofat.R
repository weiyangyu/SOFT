################################### OFAT designs with augmented points
rm(list = ls())
# Resolve paths from the SOFT repository root or any subdirectory.
soft_root <- normalizePath(".")
while (!file.exists(file.path(soft_root, "functions.R"))) {
  if (identical(soft_root, dirname(soft_root))) {
    stop("Run this script from the SOFT repository or one of its subdirectories.")
  }
  soft_root <- dirname(soft_root)
}

setwd(file.path(soft_root, "modeling", "designs"))
p = c(6,7,8,10,20)
l = 4

# MOFAT
aug_mofat_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  n_aug = l*(p[i]+1)
  lapply(1:30, function(j) {
    D <- MOFAT::mofat(p[i], l, "projection", init = "MmLHD")
    D = MOFAT::a_soft(p[i], l, n_aug=n_aug, ofat_design = D)
    return(D$design)
  })
})
save(aug_mofat_D, file = "prerun aug mofat.RData")


# SOFT - standard
aug_stand_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  n_aug = l*(p[i]+1)
  lapply(1:30, function(j) {
    D <- MOFAT::a_soft(p[i], l, n_aug = n_aug)
    return(D$design)
  })
})
save(aug_stand_D, file = "prerun aug stand.RData")


# SOFT - strict
aug_strict_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  n_aug = l*(p[i]+1)
  lapply(1:30, function(j) {
    D <- MOFAT::a_soft(p[i], l, "strict", n_aug = n_aug)
    return(D$design)
  })
})
save(aug_strict_D, file = "prerun aug strict.RData")
