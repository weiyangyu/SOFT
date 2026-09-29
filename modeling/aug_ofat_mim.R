rm(list = ls())
# Resolve paths from the SOFT repository root or any subdirectory.
soft_root <- normalizePath(".")
while (!file.exists(file.path(soft_root, "functions.R"))) {
  if (identical(soft_root, dirname(soft_root))) {
    stop("Run this script from the SOFT repository or one of its subdirectories.")
  }
  soft_root <- dirname(soft_root)
}

source(file.path(soft_root, "modeling", "lib.R"))
source(file.path(soft_root, "functions.R"))

####### load prerun design
setwd(file.path(soft_root, "modeling", "designs"))
# load("prerun aug mofat.RData")
load("prerun aug stand.RData")
load("prerun aug strict.RData")
setwd(file.path(soft_root, "modeling", "mim"))

############################################################ p = 6
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,6)
# ofat_D <- list(aug_mofat_D[[1]], aug_stand_D[[1]], aug_strict_D[[1]])
ofat_D <- list(aug_stand_D[[1]], aug_strict_D[[1]])

################################################## OTL-Circuit
otl_e <- function(X) {
  val = otlcircuit(X)
  return(val)
}
aug_otl_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = otl_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_otl_ofat_mim_0, file = "aug_otl_ofat_mim_0.RData")


############################################################ p = 7
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,7)
# ofat_D <- list(aug_mofat_D[[2]], aug_stand_D[[2]], aug_strict_D[[2]])
ofat_D <- list(aug_stand_D[[2]], aug_strict_D[[2]])

################################################## piston
piston_e <- function(X) {
  val = piston(X)
  return(val)
}
aug_piston_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = piston_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_piston_ofat_mim_0, file = "aug_piston_ofat_mim_0.RData")


############################################################ p = 8
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,8)
# ofat_D <- list(aug_mofat_D[[3]], aug_stand_D[[3]], aug_strict_D[[3]])
ofat_D <- list(aug_stand_D[[3]], aug_strict_D[[3]])

################################################## g-function
gfun_e <- function(X) {
  val = g_fun(X)
  return(val)
}
aug_gfun_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = gfun_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_gfun_ofat_mim_0, file = "aug_gfun_ofat_mim_0.RData")


################################################## levy
levy_e <- function(X) {
  val = levy(X)
  return(val)
}
aug_levy_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = levy_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_levy_ofat_mim_0, file = "aug_levy_ofat_mim_0.RData")


################################################## ackley
ackley_e <- function(X) {
  val = ackley(X)
  return(val)
}
aug_ackley_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = ackley_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_ackley_ofat_mim_0, file = "aug_ackley_ofat_mim_0.RData")


################################################## borehole
borehole_e <- function(X) {
  val = borehole(X)
  return(val)
}
aug_borehole_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = borehole_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_borehole_ofat_mim_0, file = "aug_borehole_ofat_mim_0.RData")


################################################## robot
robot_e <- function(X) {
  val = robot(X)
  return(val)
}
aug_robot_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = robot_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_robot_ofat_mim_0, file = "aug_robot_ofat_mim_0.RData")


################################################## dette
dette_e <- function(X) {
  val = Dette(X)
  return(val)
}
aug_dette_ofat_mim_0 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = dette_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_dette_ofat_mim_0, file = "aug_dette_ofat_mim_0.RData")


############################################################
############################################################
############################################################ p = 20
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,20)
# ofat_D <- list(aug_mofat_D[[length(aug_mofat_D)]], aug_stand_D[[length(aug_stand_D)]], aug_strict_D[[length(aug_strict_D)]])
ofat_D <- list(aug_stand_D[[length(aug_stand_D)]], aug_strict_D[[length(aug_strict_D)]])

################################################## OTL-Circuit
X = spacefillr::generate_sobol_set(10000,6)
y = otlcircuit(X)
sd = sqrt(var(y) / 100)
otl_e <- function(X) {
  val = otlcircuit(X)
  e = apply(X, 1, inert, p = 6, k = 14, seed = 0, sd = sd)
  return(val + e)
}
aug_otl_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = otl_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_otl_ofat_mim_1, file = "aug_otl_ofat_mim_1.RData")


################################################## piston
X = spacefillr::generate_sobol_set(10000,7)
y = piston(X)
sd = sqrt(var(y) / 100)
piston_e <- function(X) {
  val = piston(X)
  e = apply(X, 1, inert, p = 7, k = 13, seed = 0, sd = sd)
  return(val + e)
}
aug_piston_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = piston_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_piston_ofat_mim_1, file = "aug_piston_ofat_mim_1.RData")


################################################## g-function
X = spacefillr::generate_sobol_set(10000,8)
y = g_fun(X)
sd = sqrt(var(y) / 100)
gfun_e <- function(X) {
  val = g_fun(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
aug_gfun_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = gfun_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_gfun_ofat_mim_1, file = "aug_gfun_ofat_mim_1.RData")


################################################## levy
X = spacefillr::generate_sobol_set(10000,8)
y = levy(X)
sd = sqrt(var(y) / 100)
levy_e <- function(X) {
  val = levy(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
aug_levy_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = levy_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_levy_ofat_mim_1, file = "aug_levy_ofat_mim_1.RData")


################################################## ackley
X = spacefillr::generate_sobol_set(10000,8)
y = ackley(X)
sd = sqrt(var(y) / 100)
ackley_e <- function(X) {
  val = ackley(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
aug_ackley_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = ackley_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_ackley_ofat_mim_1, file = "aug_ackley_ofat_mim_1.RData")


################################################## borehole
X = spacefillr::generate_sobol_set(10000,8)
y = borehole(X)
sd = sqrt(var(y) / 100)
borehole_e <- function(X) {
  val = borehole(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
aug_borehole_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = borehole_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_borehole_ofat_mim_1, file = "aug_borehole_ofat_mim_1.RData")


################################################## robot
X = spacefillr::generate_sobol_set(10000,8)
y = robot(X)
sd = sqrt(var(y) / 100)
robot_e <- function(X) {
  val = robot(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
aug_robot_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = robot_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_robot_ofat_mim_1, file = "aug_robot_ofat_mim_1.RData")


################################################## dette
X = spacefillr::generate_sobol_set(10000,8)
y = Dette(X)
sd = sqrt(var(y) / 100)
dette_e <- function(X) {
  val = Dette(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
aug_dette_ofat_mim_1 <- lapply(1:length(ofat_D), function(i) {
  surrogate_parallel(X_ls = ofat_D[[i]], 
  X_type = "OFAT", l = 4,fun = dette_e, kernel = "MIM", 
  test_X = test_X)
})
save(aug_dette_ofat_mim_1, file = "aug_dette_ofat_mim_1.RData")
