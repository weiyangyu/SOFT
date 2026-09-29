################################### OFAT designs
rm(list = ls())
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/modeling/designs")
p = c(6,7,8,10,20)
l = 8

# MOFAT
mofat_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  lapply(1:30, function(j) {
    D <- MOFAT::mofat(p[i], l, "projection", init = "MmLHD")
    return(D)
  })
})
save(mofat_D, file = "prerun mofat.RData")


# SOFT - standard
stand_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  lapply(1:30, function(j) {
    D <- MOFAT::soft(p[i], l)
    return(D$design)
  })
})
save(stand_D, file = "prerun stand.RData")


# SOFT - strict
strict_D <- lapply(seq_along(p), function(i) {
  print(i)
  set.seed(i)
  lapply(1:30, function(j) {
    D <- MOFAT::soft(p[i], l, "strict")
    return(D$design)
  })
})
save(strict_D, file = "prerun strict.RData")
