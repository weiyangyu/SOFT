rm(list = ls())
source("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/functions.R")

####### load prerun design
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/screening/designs")
load("prerun MmLHD.RData")
load("prerun maxpro.RData")
load("prerun mofat.RData")
load("prerun stand.RData")
load("prerun strict.RData")
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/screening/results")

####### For factor screening
screen <- function(X, X_type, fun, seed) {
  set.seed(seed)
  p <- ncol(X)
  l <- nrow(X) / (p+1)
  rdx <- sample(p,p) # randomly arrange columns to factors
  X_rdx <- X[,rdx]
  y <- fun(X_rdx)

  if(X_type == "space-filling") {
    fit = rkriging::Fit.Kriging(X_rdx, y, kernel.parameters=list(type="Gaussian"))
    m = 10000
    A=matrix(runif(m*p),nrow=m)
    B=matrix(runif(m*p),nrow=m)
    colnames(A) <- paste0("X", 1:p)
    colnames(B) <- paste0("X", 1:p)
    pred_fun <- function(new_X){rkriging::Predict.Kriging(fit, new_X)$mean}
    a.sen=sensitivity::soboljansen(model = pred_fun, X1 = A, X2 = B)
    t=a.sen$T[,1]
  } else if(X_type == "OFAT") {
    t = MOFAT::measure(X_rdx, y)$t
  }
  return(t/sum(t)) # normalized total Sobol' index
}

####### reference normalized total Sobol' index
true_t <- function(fun, p) {
  m = 10000
  A=matrix(runif(m*p),nrow=m)
  B=matrix(runif(m*p),nrow=m)
  colnames(A) <- paste0("X", 1:p)
  colnames(B) <- paste0("X", 1:p)
  a.sen=sensitivity::soboljansen(model = fun, X1 = A, X2 = B)
  t=a.sen$T[,1]
  return(t/sum(t))
}

####### Visualization
plot_t_boxplots <- function(boxplot_matrices, true_t,
    plot_titles = c("MmLHD","MaxPro","MOFAT","SOFT-Stand","SOFT-Strict"),
    main_title = "Morris Function", reference_value = 1e-3) { 
  number_of_boxes <- vapply(
    boxplot_matrices,
    ncol,
    integer(1)
  )

  # Save the current graphics settings
  old_par <- par(no.readonly = TRUE)

  on.exit({
    layout(1)
    par(old_par)
  })

  # ---------------------------
  # Define the panel layout
  # ---------------------------
  # First row:  two centered plots
  # Second row: three plots
  layout_matrix <- matrix(c(0, 1, 1, 2, 2, 0,
                            3, 3, 4, 4, 5, 5),nrow = 2,byrow = TRUE)
  layout(layout_matrix, widths = rep(1, 6), heights = c(1, 1))

  # Smaller margins produce larger plotting regions.
  par(
    mar = c(3.2, 3.5, 2.2, 0.5),
    oma = c(0.2, 0.2, 2.8, 0.2),
    mgp = c(2.1, 0.65, 0),
    tcl = -0.3
  )

  # ---------------------------
  # Generate the five box plots
  # ---------------------------
  shared_ylim <- range(unlist(boxplot_matrices), true_t, reference_value, finite = TRUE)
  for (i in seq_along(boxplot_matrices)) {
    boxplot(
      boxplot_matrices[[i]],
      ylim = shared_ylim,
      main = plot_titles[i],
      cex.main = 1.3,
      cex.axis = 0.9,
      col = "lightblue",
      ylab = expression("Normalized " * t),
      cex.lab = 1.1
    )

    # Add the true T values
    points(
      x = seq_along(true_t),
      y = true_t,
      col = "red",
      cex = 1.5,
      lwd = 1.3,
      pch = 4
    )

    # Add the horizontal reference line
    abline(
      h = reference_value,
      col = "green",
      lwd = 1.3
    )
  }

  # Add the overall title
  mtext(
    main_title,
    side = 3,
    outer = TRUE,
    line = 0.7,
    cex = 1.1,
    font = 2
  )

  invisible(NULL)
}


###################################################### Morris function
true_t_morris <- true_t(sensitivity::morris.fun, 20)

MmLHD_t_morris <- c()
for(i in 1:30) {
  MmLHD_t_morris <- rbind(MmLHD_t_morris, screen(MmLHD_D[[i]], 
    "space-filling", sensitivity::morris.fun, i))
}
maxpro_t_morris <- c()
for(i in 1:30) {
  maxpro_t_morris <- rbind(maxpro_t_morris, screen(maxpro_D[[i]], 
    "space-filling", sensitivity::morris.fun, i))
}
mofat_t_morris <- c()
for(i in 1:30) {
  mofat_t_morris <- rbind(mofat_t_morris, screen(mofat_D[[i]],
    "OFAT", sensitivity::morris.fun, i))
}
stand_t_morris <- c()
for(i in 1:30) {
  stand_t_morris <- rbind(stand_t_morris, screen(stand_D[[i]], 
    "OFAT", sensitivity::morris.fun, i))
}
strict_t_morris <- c()
for(i in 1:30) {
  strict_t_morris <- rbind(strict_t_morris, screen(strict_D[[i]], 
    "OFAT", sensitivity::morris.fun, i))
}
save(true_t_morris, MmLHD_t_morris, maxpro_t_morris, 
  mofat_t_morris, stand_t_morris, strict_t_morris, file = "t_morris.RData")

load("t_morris.RData")
t_morris <- list(MmLHD_t_morris, maxpro_t_morris, mofat_t_morris, stand_t_morris, strict_t_morris)


png("t_morris.png", width = 4000, height = 3000, res = 600)
plot_t_boxplots(
  boxplot_matrices = t_morris,
  true_t = true_t_morris
)
dev.off()


###################################################### g-function
X = spacefillr::generate_sobol_set(10000,8, seed = 0)
y = g_fun(X)
sd = sqrt(var(y) / 100)
g_fun_e <- function(X) {
  val = g_fun(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}

true_t_gfun <- true_t(g_fun_e, 20)

MmLHD_t_gfun <- c()
for(i in 1:30) {
  MmLHD_t_gfun <- rbind(MmLHD_t_gfun, screen(MmLHD_D[[i]], 
    "space-filling", g_fun_e, i))
}
maxpro_t_gfun <- c()
for(i in 1:30) {
  maxpro_t_gfun <- rbind(maxpro_t_gfun, screen(maxpro_D[[i]], 
    "space-filling", g_fun_e, i))
}
mofat_t_gfun <- c()
for(i in 1:30) {
  mofat_t_gfun <- rbind(mofat_t_gfun, screen(mofat_D[[i]],
    "OFAT", g_fun_e, i))
}
stand_t_gfun <- c()
for(i in 1:30) {
  stand_t_gfun <- rbind(stand_t_gfun, screen(stand_D[[i]], 
    "OFAT", g_fun_e, i))
}
strict_t_gfun <- c()
for(i in 1:30) {
  strict_t_gfun <- rbind(strict_t_gfun, screen(strict_D[[i]], 
    "OFAT", g_fun_e, i))
}
save(true_t_gfun, MmLHD_t_gfun, maxpro_t_gfun, 
  mofat_t_gfun, stand_t_gfun, strict_t_gfun, file = )

load("t_gfun.RData")
t_gfun <- list(MmLHD_t_gfun, maxpro_t_gfun, mofat_t_gfun, stand_t_gfun, strict_t_gfun)

png("t_gfun.png", width = 4000, height = 3000, res = 600)
plot_t_boxplots(
  boxplot_matrices = t_gfun,
  true_t = true_t_gfun, main_title = "G-Function"
)
dev.off()



###################################################### levy
X = spacefillr::generate_sobol_set(10000,8, seed = 0)
y = levy(X)
sd = sqrt(var(y) / 100)
levy_e <- function(X) {
  val = levy(X)
  e = apply(X, 1, inert, p = 8, k = 12, seed = 0, sd = sd)
  return(val + e)
}

true_t_levy <- true_t(levy_e,20)

MmLHD_t_levy <- c()
for(i in 1:30) {
  MmLHD_t_levy <- rbind(MmLHD_t_levy, screen(MmLHD_D[[i]], "space-filling", levy_e, i))
}
maxpro_t_levy <- c()
for(i in 1:30) {
  maxpro_t_levy <- rbind(maxpro_t_levy, screen(maxpro_D[[i]], "space-filling", levy_e, i))
}
mofat_t_levy <- c()
for(i in 1:30) {
  mofat_t_levy <- rbind(mofat_t_levy, screen(mofat_D[[i]], "OFAT", levy_e, i))
}
stand_t_levy <- c()
for(i in 1:30) {
  stand_t_levy <- rbind(stand_t_levy, screen(stand_D[[i]], "OFAT", levy_e, i))
}
strict_t_levy <- c()
for(i in 1:30) {
  strict_t_levy <- rbind(strict_t_levy, screen(strict_D[[i]], "OFAT", levy_e, i))
}
save(true_t_levy, MmLHD_t_levy, maxpro_t_levy, 
  mofat_t_levy, stand_t_levy, strict_t_levy, file = "t_levy.RData")

load("t_levy.RData")
t_levy <- list(MmLHD_t_levy, maxpro_t_levy, mofat_t_levy, stand_t_levy, strict_t_levy)


png("t_levy.png", width = 4000, height = 3000, res = 600)
plot_t_boxplots(
  boxplot_matrices = t_levy,
  true_t = true_t_levy, main_title = "Levy"
)
dev.off()
