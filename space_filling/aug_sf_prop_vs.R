rm(list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/space_filling/designs")

####### load prerun design
load("prerun MmLHD.RData")
load("prerun maxpro.RData")
load("prerun mofat.RData")
load("prerun stand.RData")
load("prerun strict.RData")
# load("prerun aug mofat.RData")
load("prerun aug stand.RData")
load("prerun aug strict.RData")

############################# Q
l = 6
beta = 1/(2*l)^2

# aug_mofat_Q <- lapply(1:19, function(i) {
#   print(i)
#   set.seed(i)
#   p <- i+1
#   qmc <- MOFAT::qmc_generate(10000, p) 
#   a <- lapply(1:30, function(j) {
#     MOFAT::Q(1, beta, aug_mofat_D[[i]][[j]], qmc)$Q
#   })
#   a <- unlist(a)
#   return(a)
# })
# aug_mofat_Q <- Reduce(cbind, aug_mofat_Q)


aug_stand_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, aug_stand_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
aug_stand_Q <- Reduce(cbind, aug_stand_Q)


aug_strict_Q <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  qmc <- MOFAT::qmc_generate(10000, p)
  a <- lapply(1:30, function(j) {
    MOFAT::Q(1, beta, aug_strict_D[[i]][[j]], qmc)$Q
  })
  a <- unlist(a)
  return(a)
})
aug_strict_Q <- Reduce(cbind, aug_strict_Q)

setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/space_filling/visualization")
# save(aug_mofat_Q, aug_stand_Q, aug_strict_Q, file = "aug_Q.RData")
save(aug_stand_Q, aug_strict_Q, file = "aug_Q.RData")


############################# entropy (Phi)
entropy_fun <- function(D,beta) {
  p = ncol(D)
  gamma = p*beta
  ans <- sum(1/(1 + dist(D)^2 / gamma)) # minimize
  return(1/ans) # maximize
}

# aug_mofat_ent <- lapply(1:19, function(i) {
#   print(i)
#   set.seed(i)
#   p <- i+1
#   a <- lapply(1:30, function(j) {
#     entropy_fun(aug_mofat_D[[i]][[j]], beta)
#   })
#   a <- unlist(a)
#   return(a)
# })
# aug_mofat_ent <- Reduce(cbind, aug_mofat_ent)


aug_stand_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(aug_stand_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
aug_stand_ent <- Reduce(cbind, aug_stand_ent)


aug_strict_ent <- lapply(1:19, function(i) {
  print(i)
  set.seed(i)
  p <- i+1
  a <- lapply(1:30, function(j) {
    entropy_fun(aug_strict_D[[i]][[j]], beta)
  })
  a <- unlist(a)
  return(a)
})
aug_strict_ent <- Reduce(cbind, aug_strict_ent)

# save(aug_mofat_ent, aug_stand_ent, aug_strict_ent, file = "aug_ent.RData")
save(aug_stand_ent, aug_strict_ent, file = "aug_ent.RData")


############################# computational time
MmLHD_time = Reduce(cbind, MmLHD_time)
maxpro_time = Reduce(cbind, maxpro_time)
mofat_time = Reduce(cbind, mofat_time)
stand_time = Reduce(cbind, stand_time)
strict_time = Reduce(cbind, strict_time)

# aug_mofat_time = Reduce(cbind, aug_mofat_time)
aug_stand_time = Reduce(cbind, aug_stand_time)
aug_strict_time = Reduce(cbind, aug_strict_time)


load("Q.RData")
load("ent.RData")
load("aug_Q.RData")
load("aug_ent.RData")

############################# visualization
plot_SF <- function(baseline, all_list,
                    method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT.Stand", "SOFT.Strict"),
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
    ggplot2::scale_shape_manual(values = c(18, 8, 16, 17, 15, 1, 2)) +
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
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT.Stand", "SOFT.Strict"),
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
    ggplot2::scale_shape_manual(values = c(18, 8, 16, 17, 15, 1, 2)) +
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
  all_list = list(MmLHD_Q, maxpro_Q, mofat_Q, stand_Q, strict_Q, aug_stand_Q, aug_strict_Q),
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict", 
"A-SOFT-Stand", "A-SOFT-Strict"),
  y_label = expression(log(Q / Median(Q[MmLHD]))),
  plot_title = expression(bold("(a) " * Q / Median(Q[MmLHD]))),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
phi_plot <- plot_SF(
  baseline = MmLHD_ent,
  all_list = list(MmLHD_ent, maxpro_ent, mofat_ent, stand_ent, strict_ent, aug_stand_ent, aug_strict_ent),
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict", 
"A-SOFT-Stand", "A-SOFT-Strict"),
  y_label = expression(log(Median(Phi[MmLHD]) / Phi)),
  plot_title = expression(bold("(b) " * Median(Phi[MmLHD]) / Phi)),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)
runtime_plot <- plot_runtime(all_list = list(MmLHD_time, maxpro_time, mofat_time, stand_time, strict_time, aug_stand_time, aug_strict_time),
  method_names = c("MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict", 
"A-SOFT-Stand", "A-SOFT-Strict"),
  y_label = expression(Runtime (sec)),
  plot_title = expression(bold("(c) " * Runtime)),
  colors = c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff")
)

# Open the output device before building grobs to avoid an implicit Rplots.pdf.
local({
  grDevices::png("aug_combined_sf.png", width = 9, height = 3,
    units = "in", res = 600, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)

  # Use one shared legend, on one line for all seven designs.
  legend_plot <- q_plot + ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 1, byrow = TRUE),
    shape = ggplot2::guide_legend(nrow = 1, byrow = TRUE)
  ) + ggplot2::theme(
    legend.text = ggplot2::element_text(size = 10),
    legend.key.width = grid::unit(0.35, "cm"),
    legend.spacing.x = grid::unit(0.08, "cm")
  )
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
combined_image <- magick::image_read("aug_combined_sf.png")
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
magick::image_write(trimmed_combined, path = "aug_combined_sf.png", density = "600x600")
