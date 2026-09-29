# Generate designs for factor screening
# p = 20
# l = 8

rm(list = ls())
# Resolve paths from the SOFT repository root or any subdirectory.
soft_root <- normalizePath(".")
while (!file.exists(file.path(soft_root, "functions.R"))) {
  if (identical(soft_root, dirname(soft_root))) {
    stop("Run this script from the SOFT repository or one of its subdirectories.")
  }
  soft_root <- dirname(soft_root)
}

setwd(file.path(soft_root, "screening", "designs"))
p = 20
l = 8
n = l*(p+1)

# MmLHD
MmLHD_D <- lapply(1:30, function(i) {
  print(i)
  D <- SLHD::maximinSLHD(1, n, p)$StandDesign
  D
})
save(MmLHD_D, file = "prerun MmLHD.RData")

# maxpro
maxpro_D <- lapply(1:30, function(i) {
  print(i)
  InitialDesign <- MaxPro::MaxProLHD(n = n, p = p)$Design 
  D <- MaxPro::MaxPro(InitialDesign)$Design
  D
})
save(maxpro_D, file = "prerun maxpro.RData")

# MOFAT
mofat_D <- lapply(1:30, function(i) {
  print(i)
  D <- MOFAT::mofat(p, l, "projection", init = "MmLHD")
  D
})
save(mofat_D, file = "prerun mofat.RData")

# soft - standard
stand_D <- lapply(1:30, function(i) {
  print(i)
  D <- MOFAT::soft(p, l)
  D <- D$design
  D
})
save(stand_D, file = "prerun stand.RData")

# soft - strict
strict_D <- lapply(1:30, function(i) {
  print(i)
  D <- MOFAT::soft(p, l, structure = "strict")
  D <- D$design
  D
})
save(strict_D, file = "prerun strict.RData")
