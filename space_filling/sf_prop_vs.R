rm(list = ls())
# Resolve paths from the SOFT repository root or any subdirectory.
soft_root <- normalizePath(".")
while (!file.exists(file.path(soft_root, "functions.R"))) {
  if (identical(soft_root, dirname(soft_root))) {
    stop("Run this script from the SOFT repository or one of its subdirectories.")
  }
  soft_root <- dirname(soft_root)
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
setwd(file.path(soft_root, "space_filling", "designs"))

####### load prerun design
load("prerun MmLHD.RData")
load("prerun maxpro.RData")
load("prerun mofat.RData")
load("prerun stand.RData")
load("prerun strict.RData")

############################# Q
l = 6
beta = 1/(2*l)^2

MmLHD_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, MmLHD_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
MmLHD_Q <- Reduce(cbind, MmLHD_Q)


maxpro_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, maxpro_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
maxpro_Q <- Reduce(cbind, maxpro_Q)


mofat_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, mofat_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
mofat_Q <- Reduce(cbind, mofat_Q)


stand_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, stand_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
stand_Q <- Reduce(cbind, stand_Q)


strict_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, strict_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
strict_Q <- Reduce(cbind, strict_Q)

setwd(file.path(soft_root, "space_filling", "visualization"))
save(MmLHD_Q, maxpro_Q, mofat_Q, stand_Q, strict_Q, file = "Q.RData")
load("Q.RData")


############################# entropy (Phi)
entropy_fun <- function(D,beta) {
  p = ncol(D)
  gamma = p*beta
  ans <- sum(1/(1 + dist(D)^2 / gamma)) # minimize
  return(1/ans) # maximize
}

MmLHD_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(MmLHD_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
MmLHD_ent <- Reduce(cbind, MmLHD_ent)


maxpro_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(maxpro_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
maxpro_ent <- Reduce(cbind, maxpro_ent)


mofat_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(mofat_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
mofat_ent <- Reduce(cbind, mofat_ent)


stand_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(stand_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
stand_ent <- Reduce(cbind, stand_ent)


strict_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(strict_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
strict_ent <- Reduce(cbind, strict_ent)

save(MmLHD_ent, maxpro_ent, mofat_ent, stand_ent, strict_ent, file = "ent.RData")
load("ent.RData")


############################# computational time
MmLHD_time = Reduce(cbind, MmLHD_time)
maxpro_time = Reduce(cbind, maxpro_time)
mofat_time = Reduce(cbind, mofat_time)
stand_time = Reduce(cbind, stand_time)
strict_time = Reduce(cbind, strict_time)


############################# visualization
plot_SF <- function(baseline, all_list,
                    method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict"),
                    factor_seq = 2:20,
                    y_label, 
                    plot_title,
                    colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0")) {
  baseline_median = apply(baseline, 2, median)
  one_method <- function(mat, method_name) {
    LR <- log(sweep(mat, 2, baseline_median, "/"))
    data.frame(factor = factor_seq, method = method_name, 
      Q1 = apply(LR,2,quantile)[2,], Q2 = apply(LR,2,quantile)[3,],
      Q3 = apply(LR,2,quantile)[4,])
  }

  sum_df <- dplyr::bind_rows(Map(one_method, all_list, method_names)) |>
    dplyr::mutate(
      method = factor(method, levels = method_names)
    )

  names(colors) <- method_names

  p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = factor, y = Q2, color = method, group = method, linetype = method)) +
    ggplot2::geom_line(ggplot2::aes(linetype = method), linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(shape = method), size = 1.1) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Q1, ymax = Q3),
      width = 0.12, linewidth = 0.25
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(breaks = factor_seq[factor_seq %% 2 == 0]) + 
    ggplot2::labs(x = "No. of Factors (p)", y = y_label, title = plot_title) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(), legend.position = "top",
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
    )

  invisible(p)
}


plot_runtime <- function(all_list,
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict"),
  factor_seq = 2:20, y_label, plot_title,
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0")) {

  sum_df <- dplyr::bind_rows(Map(
    function(mat, method_name) {
        data.frame(factor = factor_seq,method = method_name,
          Q1 = apply(mat,2,quantile)[2,], Q2 = apply(mat,2,quantile)[3,],
          Q3 = apply(mat,2,quantile)[4,])
    }, all_list, method_names)) |>
    dplyr::mutate(
      method = factor(method, levels = method_names)
    )

  names(colors) <- method_names

  p <- ggplot2::ggplot(sum_df, ggplot2::aes(x = factor, y = Q2, color = method, group = method, linetype = method)) +
    ggplot2::geom_line(ggplot2::aes(linetype = method), linewidth = 0.35) +
    ggplot2::geom_point(ggplot2::aes(shape = method), size = 1.1) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = Q1, ymax = Q3),
      width = 0.12, linewidth = 0.25
    ) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(breaks = factor_seq[factor_seq %% 2 == 0]) + 
    ggplot2::labs(x = "No. of Factors (p)", y = y_label, title = plot_title) +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
      linetype = ggplot2::guide_legend(nrow = 3, byrow = TRUE),
      shape = ggplot2::guide_legend(nrow = 3, byrow = TRUE)
    ) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(), legend.position = "top",
      legend.text = ggplot2::element_text(size = 8),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11)
    )

  invisible(p)
}

q_plot <- plot_SF(
  baseline = MmLHD_Q,
  all_list = list(MmLHD_Q, maxpro_Q, mofat_Q, stand_Q, strict_Q),
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict"),
  y_label = expression(log(Q / Median(Q[MmLHD]))),
  plot_title = expression(bold("(a) " * Q / Median(Q[MmLHD])))
)
phi_plot <- plot_SF(
  baseline = MmLHD_ent,
  all_list = list(MmLHD_ent, maxpro_ent, mofat_ent, stand_ent, strict_ent),
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict"),
  y_label = expression(log(Median(Phi[MmLHD]) / Phi)),
  plot_title = expression(bold("(b) " * Median(Phi[MmLHD]) / Phi))
)

runtime_plot <- plot_runtime(all_list = list(MmLHD_time, maxpro_time, mofat_time, stand_time, strict_time),
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict"),
  y_label = expression(Runtime (sec)),
  plot_title = expression(bold("(c) " * Runtime))
)

# Open the output device before building grobs to avoid an implicit Rplots.pdf.
local({
  grDevices::png("combined_sf.png", width = 9, height = 3,
    units = "in", res = 600, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)

  # Use one horizontal legend for the combined figure.
  legend_plot <- q_plot + ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1),
    linetype = ggplot2::guide_legend(nrow = 1),
    shape = ggplot2::guide_legend(nrow = 1)
  ) + ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
  legend_gtable <- ggplot2::ggplotGrob(legend_plot)
  legend_index <- which(
    grepl("^guide-box", legend_gtable$layout$name) &
      !vapply(legend_gtable$grobs, inherits, logical(1), "zeroGrob")
  )
  shared_legend <- legend_gtable$grobs[[legend_index[1]]]
  panels <- lapply(list(q_plot, phi_plot, runtime_plot), function(p) {
    ggplot2::ggplotGrob(p + ggplot2::theme(legend.position = "none"))
  })
  # Align plotting regions by padding the outer margin, keeping each
  # y-axis title next to its own tick labels rather than widening that gap.
  left_widths <- lapply(panels, function(p) {
    panel_col <- p$layout$l[p$layout$name == "panel"]
    sum(p$widths[seq_len(panel_col - 1)])
  })
  max_left_width <- do.call(grid::unit.pmax, left_widths)
  panels <- Map(function(p, left_width) {
    p$widths[1] <- p$widths[1] + max_left_width - left_width
    p
  }, panels, left_widths)
  combined <- gridExtra::arrangeGrob(
    shared_legend,
    gridExtra::arrangeGrob(grobs = panels, nrow = 1),
    ncol = 1,
    heights = grid::unit.c(sum(shared_legend$heights), grid::unit(1, "null"))
  )
  grid::grid.newpage()
  grid::grid.draw(combined)
})

# Trim the left, right, and bottom while preserving the original top margin.
combined_image <- magick::image_read("combined_sf.png")
image_pixels <- magick::image_data(combined_image, channels = "rgb")
content_rows <- which(apply(image_pixels != as.raw(255), 3, any))
top_margin <- if (length(content_rows)) content_rows[1] - 1L else 0L
trimmed_combined <- magick::image_trim(combined_image)
trimmed_info <- magick::image_info(trimmed_combined)
trimmed_combined <- magick::image_extent(
  trimmed_combined,
  geometry = sprintf("%dx%d", trimmed_info$width, trimmed_info$height + top_margin),
  gravity = "south", color = "white"
)
magick::image_write(trimmed_combined, path = "combined_sf.png", density = "600x600")
