rm(list = ls())
# Resolve paths from the SOFT repository root or any subdirectory.
soft_root <- normalizePath(".")
while (!file.exists(file.path(soft_root, "functions.R"))) {
  if (identical(soft_root, dirname(soft_root))) {
    stop("Run this script from the SOFT repository or one of its subdirectories.")
  }
  soft_root <- dirname(soft_root)
}

setwd(file.path(soft_root, "modeling", "gauss"))

######################### load mse with the original dimension
load("ackley_sfd_gauss_0.RData")
load("ackley_ofat_gauss_0.RData")
load("aug_ackley_ofat_gauss_0.RData")
load("borehole_sfd_gauss_0.RData")
load("borehole_ofat_gauss_0.RData")
load("aug_borehole_ofat_gauss_0.RData")
load("dette_sfd_gauss_0.RData")
load("dette_ofat_gauss_0.RData")
load("aug_dette_ofat_gauss_0.RData")
load("gfun_sfd_gauss_0.RData")
load("gfun_ofat_gauss_0.RData")
load("aug_gfun_ofat_gauss_0.RData")
load("levy_sfd_gauss_0.RData")
load("levy_ofat_gauss_0.RData")
load("aug_levy_ofat_gauss_0.RData")
load("otl_sfd_gauss_0.RData")
load("otl_ofat_gauss_0.RData")
load("aug_otl_ofat_gauss_0.RData")
load("piston_sfd_gauss_0.RData")
load("piston_ofat_gauss_0.RData")
load("aug_piston_ofat_gauss_0.RData")
load("robot_sfd_gauss_0.RData")
load("robot_ofat_gauss_0.RData")
load("aug_robot_ofat_gauss_0.RData")
load("wing_sfd_gauss_0.RData")
load("wing_ofat_gauss_0.RData")
load("aug_wing_ofat_gauss_0.RData")


######################### load mse with the augmented dimension
load("ackley_sfd_gauss_1.RData")
load("ackley_ofat_gauss_1.RData")
load("aug_ackley_ofat_gauss_1.RData")
load("borehole_sfd_gauss_1.RData")
load("borehole_ofat_gauss_1.RData")
load("aug_borehole_ofat_gauss_1.RData")
load("dette_sfd_gauss_1.RData")
load("dette_ofat_gauss_1.RData")
load("aug_dette_ofat_gauss_1.RData")
load("gfun_sfd_gauss_1.RData")
load("gfun_ofat_gauss_1.RData")
load("aug_gfun_ofat_gauss_1.RData")
load("levy_sfd_gauss_1.RData")
load("levy_ofat_gauss_1.RData")
load("aug_levy_ofat_gauss_1.RData")
load("otl_sfd_gauss_1.RData")
load("otl_ofat_gauss_1.RData")
load("aug_otl_ofat_gauss_1.RData")
load("piston_sfd_gauss_1.RData")
load("piston_ofat_gauss_1.RData")
load("aug_piston_ofat_gauss_1.RData")
load("robot_sfd_gauss_1.RData")
load("robot_ofat_gauss_1.RData")
load("aug_robot_ofat_gauss_1.RData")
load("wing_sfd_gauss_1.RData")
load("wing_ofat_gauss_1.RData")
load("aug_wing_ofat_gauss_1.RData")



###########################################################################
# Extract the two SFD, three OFAT MSE, and the two augmented OFAT MSE samples in plotting order.
extract_mse_matrices <- function(sfd_gauss, ofat_gauss, aug_ofat_gauss) {
  if (!is.list(sfd_gauss) || length(sfd_gauss) < 2) {
    stop("`sfd_gauss` must be a list containing at least two results.")
  }
  if (!is.list(ofat_gauss) || length(ofat_gauss) < 3) {
    stop("`ofat_gauss` must be a list containing at least three results.")
  }

  if (!is.list(aug_ofat_gauss) || length(aug_ofat_gauss) != 2) {
    stop("`aug_ofat_gauss` must contain Aug.SOFT.Standard and Aug.SOFT.Strict, in that order.")
  }

  results <- c(sfd_gauss[1:2], ofat_gauss[1:3], aug_ofat_gauss[1:2])
  method_names <- c("MmLHD", "MaxPro", "MOFAT", "SOFT.Standard", "SOFT.Strict", "Aug.SOFT.Standard", "Aug.SOFT.Strict")

  mse_values <- lapply(seq_along(results), function(i) {
    mse_mat <- results[[i]]$mse_mat
    if (is.null(mse_mat) || !is.numeric(mse_mat)) {
      stop("Missing or nonnumeric `mse_mat` for ", method_names[i], ".")
    }

    values <- as.numeric(mse_mat)
    values <- values[is.finite(values)]
    if (length(values) == 0) {
      stop("`mse_mat` for ", method_names[i], " has no finite values.")
    }
    values
  })

  setNames(mse_values, method_names)
}


# Draw one two-panel log-MSE row. Each panel contains seven method boxplots.
plot_mse_row <- function(
    sfd_gauss_0,
    ofat_gauss_0,
    aug_ofat_gauss_0,
    sfd_gauss_1,
    ofat_gauss_1,
    aug_ofat_gauss_1,
    panel_titles = c("p = 8", "p = 8 + 12"),
    method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
    colors = c(
      "#5f5f5fCC",
      "#2e8b57CC",
      "#2f38ffCC",
      "#f4a261CC",
      "#cc45e0CC",
      "#b38b63CC",
      "#b980ff"
    ),
    row_title = NULL,
    ylab = expression(log(MSE)),
    ylims = NULL,
    outline = FALSE
) {
  if (length(method_labels) != 7 || length(colors) != 7) {
    stop("Provide exactly seven method labels and colors.")
  }
  if (is.null(ylims)) {
    ylims <- list(NULL, NULL)
  }
  if (!is.list(ylims) || length(ylims) != 2) {
    stop("`ylims` must be NULL or a list containing two panel-specific limits.")
  }
  valid_ylim <- function(x) {
    is.null(x) ||
      (is.numeric(x) && length(x) == 2 && all(is.finite(x)) && x[1] < x[2])
  }
  if (!all(vapply(ylims, valid_ylim, logical(1)))) {
    stop("Each non-NULL element of `ylims` must contain two increasing finite numbers.")
  }

  panel_data <- list(
    extract_mse_matrices(sfd_gauss_0, ofat_gauss_0, aug_ofat_gauss_0),
    extract_mse_matrices(sfd_gauss_1, ofat_gauss_1, aug_ofat_gauss_1)
  )
  if (any(unlist(panel_data, use.names = FALSE) <= 0)) {
    stop("All MSE values must be positive to calculate log(MSE).")
  }
  panel_data <- lapply(
    panel_data,
    function(panel) lapply(panel, log)
  )

  # Calculate each panel's limits independently. When outliers are hidden,
  # use the boxplot whiskers so hidden extremes do not compress the boxes.
  calculate_ylim <- function(samples) {
    scale_values <- if (outline) {
      unlist(samples, use.names = FALSE)
    } else {
      unlist(lapply(samples, function(x) boxplot.stats(x)$stats), use.names = FALSE)
    }
    panel_ylim <- range(scale_values, finite = TRUE)

    if (diff(panel_ylim) == 0) {
      padding <- if (panel_ylim[1] == 0) 0.5 else abs(panel_ylim[1]) * 0.04
      panel_ylim <- panel_ylim + c(-padding, padding)
    }
    panel_ylim
  }

  for (i in seq_along(panel_data)) {
    if (is.null(ylims[[i]])) {
      ylims[[i]] <- calculate_ylim(panel_data[[i]])
    }
  }

  old_par <- par(no.readonly = TRUE)
  on.exit({
    layout(1)
    par(old_par)
  })

  layout(matrix(1:2, nrow = 1))
  par(
    mar = c(4.2, 4.4, 2.2, 0.8),
    oma = c(0, 0, if (is.null(row_title)) 0 else 1.8, 0),
    mgp = c(2.5, 0.7, 0),
    tcl = -0.3
  )

  for (i in seq_along(panel_data)) {
    boxplot(
      panel_data[[i]],
      names = FALSE,
      main = panel_titles[i],
      col = colors,
      border = "black",
      ylab = ylab,
      ylim = ylims[[i]],
      outline = outline,
      xaxt = "n",
      cex.axis = 1.3,
      cex.lab = 1.5,
      cex.main = 1.6
    )
    abline(v = 3.5, col = "red", lty = 2)
    axis(1, at = seq_along(method_labels), labels = FALSE)
    mtext(
      method_labels,
      side = 1,
      at = seq_along(method_labels),
      line = 2.5,
      cex = 0.9
    )
  }

  if (!is.null(row_title)) {
    mtext(row_title, side = 3, outer = TRUE, line = 0.1, cex = 1.8, font = 2)
  }

  invisible(list(
    log_mse_0 = panel_data[[1]],
    log_mse_1 = panel_data[[2]],
    ylims = ylims
  ))
}

setwd(file.path(soft_root, "modeling", "visual_gauss"))
png("aug_ackley_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  ackley_sfd_gauss_0,
  ackley_ofat_gauss_0,
  aug_ackley_ofat_gauss_0,
  ackley_sfd_gauss_1,
  ackley_ofat_gauss_1,
  aug_ackley_ofat_gauss_1,
  row_title = "Ackley",
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_borehole_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  borehole_sfd_gauss_0,
  borehole_ofat_gauss_0,
  aug_borehole_ofat_gauss_0,
  borehole_sfd_gauss_1,
  borehole_ofat_gauss_1,
  aug_borehole_ofat_gauss_1,
  row_title = "Borehole",
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_dette_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  dette_sfd_gauss_0,
  dette_ofat_gauss_0,
  aug_dette_ofat_gauss_0,
  dette_sfd_gauss_1,
  dette_ofat_gauss_1,
  aug_dette_ofat_gauss_1,
  row_title = "Dette and Pepelyshev",
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_gfun_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  gfun_sfd_gauss_0,
  gfun_ofat_gauss_0,
  aug_gfun_ofat_gauss_0,
  gfun_sfd_gauss_1,
  gfun_ofat_gauss_1,
  aug_gfun_ofat_gauss_1,
  row_title = "G-function",
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_levy_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  levy_sfd_gauss_0,
  levy_ofat_gauss_0,
  aug_levy_ofat_gauss_0,
  levy_sfd_gauss_1,
  levy_ofat_gauss_1,
  aug_levy_ofat_gauss_1,
  row_title = "Levy",
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_otl_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  otl_sfd_gauss_0,
  otl_ofat_gauss_0,
  aug_otl_ofat_gauss_0,
  otl_sfd_gauss_1,
  otl_ofat_gauss_1,
  aug_otl_ofat_gauss_1,
  row_title = "OTL-Circuit",
  panel_titles = c("p = 6", "p = 6 + 14"),
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_piston_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  piston_sfd_gauss_0,
  piston_ofat_gauss_0,
  aug_piston_ofat_gauss_0,
  piston_sfd_gauss_1,
  piston_ofat_gauss_1,
  aug_piston_ofat_gauss_1,
  row_title = "Piston",
  panel_titles = c("p = 7", "p = 7 + 13"),
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_robot_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  robot_sfd_gauss_0,
  robot_ofat_gauss_0,
  aug_robot_ofat_gauss_0,
  robot_sfd_gauss_1,
  robot_ofat_gauss_1,
  aug_robot_ofat_gauss_1,
  row_title = "Robot",
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()

png("aug_wing_mse_gauss.png", width = 8000, height = 2700, res = 600)
plot_mse_row(
  wing_sfd_gauss_0,
  wing_ofat_gauss_0,
  aug_wing_ofat_gauss_0,
  wing_sfd_gauss_1,
  wing_ofat_gauss_1,
  aug_wing_ofat_gauss_1,
  row_title = "Wingweight",
  panel_titles = c("p = 10", "p = 10 + 10"),
  method_labels = c("MmLHD", "MaxPro", "MOFAT", "SOFT\nStand", "SOFT\nStrict", "A-SOFT\nStand", "A-SOFT\nStrict"),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
dev.off()


# library(dplyr)
# library(tidyr)
# library(ggplot2)
# library(reshape2)
# library(png)
# library(grid)
# library(gridExtra)
# img1 <- rasterGrob(readPNG("aug_levy_mse.png"), interpolate = TRUE)
# img2 <- rasterGrob(readPNG("aug_ackley_mse.png"), interpolate = TRUE)
# png("aug_combined_mse.png", width = 8000, height = 5400, res = 600)
# grid.arrange(img1, img2, nrow = 2)
# dev.off()
