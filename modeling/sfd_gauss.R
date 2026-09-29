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
load("prerun MmLHD.RData")
load("prerun maxpro.RData")
setwd(file.path(soft_root, "modeling", "gauss"))

############################################################ p = 6
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,6)
sf_D <- list(MmLHD_D[[1]], maxpro_D[[1]])

################################################## OTL-Circuit
otl_e <- function(X) {
  val = otlcircuit(X)
  return(val)
}
otl_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = otl_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(otl_sfd_gauss_0, file = "otl_sfd_gauss_0.RData")


############################################################ p = 7
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,7)
sf_D <- list(MmLHD_D[[2]], maxpro_D[[2]])

################################################## piston
piston_e <- function(X) {
  val = piston(X)
  return(val)
}
piston_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = piston_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(piston_sfd_gauss_0, file = "piston_sfd_gauss_0.RData")


############################################################ p = 8
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,8)
sf_D <- list(MmLHD_D[[3]], maxpro_D[[3]])

################################################## g-function
gfun_e <- function(X) {
  val = g_fun(X)
  return(val)
}
gfun_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = gfun_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(gfun_sfd_gauss_0, file = "gfun_sfd_gauss_0.RData")


################################################## levy
levy_e <- function(X) {
  val = levy(X)
  return(val)
}
levy_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = levy_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(levy_sfd_gauss_0, file = "levy_sfd_gauss_0.RData")


################################################## ackley
ackley_e <- function(X) {
  val = ackley(X)
  return(val)
}
ackley_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = ackley_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(ackley_sfd_gauss_0, file = "ackley_sfd_gauss_0.RData")


################################################## borehole
borehole_e <- function(X) {
  val = borehole(X)
  return(val)
}
borehole_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = borehole_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(borehole_sfd_gauss_0, file = "borehole_sfd_gauss_0.RData")


################################################## robot
robot_e <- function(X) {
  val = robot(X)
  return(val)
}
robot_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = robot_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(robot_sfd_gauss_0, file = "robot_sfd_gauss_0.RData")


################################################## dette
dette_e <- function(X) {
  val = Dette(X)
  return(val)
}
dette_sfd_gauss_0 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = dette_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(dette_sfd_gauss_0, file = "dette_sfd_gauss_0.RData")


############################################################
############################################################
############################################################ p = 20
set.seed(0)
test_X <- MOFAT::qmc_generate(10000,20)
sf_D <- list(MmLHD_D[[length(MmLHD_D)]], maxpro_D[[length(maxpro_D)]])

################################################## OTL-Circuit
X = spacefillr::generate_sobol_set(10000,6)
y = otlcircuit(X)
sd = sqrt(var(y) / 100)
otl_e <- function(X) {
  val = otlcircuit(X)
  e = apply(X, 1, inert, p = 6, k = 14, seed = 0, sd = sd)
  return(val + e)
}
otl_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = otl_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(otl_sfd_gauss_1, file = "otl_sfd_gauss_1.RData")


################################################## piston
X = spacefillr::generate_sobol_set(10000,7)
y = piston(X)
sd = sqrt(var(y) / 100)
piston_e <- function(X) {
  val = piston(X)
  e = apply(X, 1, inert, p = 7, k = 13, seed = 0, sd = sd)
  return(val + e)
}
piston_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = piston_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(piston_sfd_gauss_1, file = "piston_sfd_gauss_1.RData")


################################################## g-function
X = spacefillr::generate_sobol_set(10000,8)
y = g_fun(X)
sd = sqrt(var(y) / 100)
gfun_e <- function(X) {
  val = g_fun(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
gfun_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = gfun_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(gfun_sfd_gauss_1, file = "gfun_sfd_gauss_1.RData")


################################################## levy
X = spacefillr::generate_sobol_set(10000,8)
y = levy(X)
sd = sqrt(var(y) / 100)
levy_e <- function(X) {
  val = levy(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
levy_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = levy_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(levy_sfd_gauss_1, file = "levy_sfd_gauss_1.RData")


################################################## ackley
X = spacefillr::generate_sobol_set(10000,8)
y = ackley(X)
sd = sqrt(var(y) / 100)
ackley_e <- function(X) {
  val = ackley(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
ackley_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = ackley_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(ackley_sfd_gauss_1, file = "ackley_sfd_gauss_1.RData")


################################################## borehole
X = spacefillr::generate_sobol_set(10000,8)
y = borehole(X)
sd = sqrt(var(y) / 100)
borehole_e <- function(X) {
  val = borehole(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
borehole_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = borehole_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(borehole_sfd_gauss_1, file = "borehole_sfd_gauss_1.RData")


################################################## robot
X = spacefillr::generate_sobol_set(10000,8)
y = robot(X)
sd = sqrt(var(y) / 100)
robot_e <- function(X) {
  val = robot(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
robot_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = robot_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(robot_sfd_gauss_1, file = "robot_sfd_gauss_1.RData")


################################################## dette
X = spacefillr::generate_sobol_set(10000,8)
y = Dette(X)
sd = sqrt(var(y) / 100)
dette_e <- function(X) {
  val = Dette(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}
dette_sfd_gauss_1 <- lapply(1:2, function(i) {
  surrogate_parallel(X_ls = sf_D[[i]], 
  X_type = "space-filling", l = 0,fun = dette_e, kernel = "Gaussian", 
  test_X = test_X)
})
save(dette_sfd_gauss_1, file = "dette_sfd_gauss_1.RData")
