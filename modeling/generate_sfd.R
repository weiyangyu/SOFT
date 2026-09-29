################################### Space-Filling designs
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
l = 8

# MmLHD
MmLHD_D <- lapply(1:length(p), function(i) {
  print(i)
  set.seed(i)
  n = l * (p[i]+1)
  lapply(1:30, function(j) {
    SLHD::maximinSLHD(1, n, p[i])$StandDesign
  })
})
save(MmLHD_D, file = "prerun MmLHD.RData")

# MaxPro
maxpro_D <- lapply(1:length(p), function(i) {
  print(i)
  set.seed(i)
  n = l * (p[i]+1)
  lapply(1:30, function(j) {
    InitialDesign <- MaxPro::MaxProLHD(n = n, p = p[i])$Design 
    D <- MaxPro::MaxPro(InitialDesign)$Design
    return(D)
  })
})
save(maxpro_D, file = "prerun maxpro.RData")