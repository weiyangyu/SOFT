################################### OFAT designs with augmented points
rm(list = ls())
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/modeling/designs")
p = c(6,7,8,10,20)
l = 4

# MOFAT
aug_mofat_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  n_aug = l*(p[i]+1)
  lapply(1:30, function(j) {
    D <- MOFAT::mofat(p[i], l, "projection", init = "MmLHD")
    D = MOFAT::soft_augment(p[i], l, n_aug=n_aug, existing_design = D)
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
    D <- MOFAT::soft_augment(p[i], l, n_aug = n_aug)
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
    D <- MOFAT::soft_augment(p[i], l, "strict", n_aug = n_aug)
    return(D$design)
  })
})
save(aug_strict_D, file = "prerun aug strict.RData")
