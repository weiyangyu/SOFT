rm(list = ls())
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/space_filling/visualization")
load("Q.RData")
load("aug_Q.RData")
Q_vals <- c(1, median(maxpro_Q[,9])/median(MmLHD_Q[,9]), median(mofat_Q[,9])/median(MmLHD_Q[,9]), 
  median(stand_Q[,9])/median(MmLHD_Q[,9]), median(strict_Q[,9])/median(MmLHD_Q[,9]),
  median(aug_stand_Q[,9])/median(MmLHD_Q[,9]), median(aug_strict_Q[,9])/median(MmLHD_Q[,9]))

####### load prerun design
setwd("/Users/weiyang/Library/CloudStorage/OneDrive-GeorgiaInstituteofTechnology/OFAT Design/code/space_filling/designs")
load("prerun MmLHD.RData")
load("prerun maxpro.RData")
load("prerun mofat.RData")
load("prerun stand.RData")
load("prerun strict.RData")
load("prerun aug stand.RData")
load("prerun aug strict.RData")

# p = 20
# l = 6
l = 6
l2 = 4
mofat_perturb <- sum(abs(mofat_D[[9]][[1]][1:l,1] - mofat_D[[9]][[1]][(l+1):(2*l),1]))
stand_perturb <- sum(abs(stand_D[[9]][[1]][1:l,1] - stand_D[[9]][[1]][(l+1):(2*l),1]))
strict_diff_mat <- strict_D[[9]][[1]][1:l,] - strict_D[[9]][[1]][(l+1):(2*l),]
a = which(colSums(abs(strict_diff_mat) > 1e-10) > 0)
strict_perturb <- sum(abs(strict_diff_mat[,a]))

aug_stand_perturb <- sum(abs(aug_stand_D[[9]][[1]][1:l2,1] - aug_stand_D[[9]][[1]][(l2+1):(2*l2),1]))
aug_strict_diff_mat <- aug_strict_D[[9]][[1]][1:l2,] - aug_strict_D[[9]][[1]][(l2+1):(2*l2),]
a = which(colSums(abs(aug_strict_diff_mat) > 1e-10) > 0)
aug_strict_perturb <- sum(abs(aug_strict_diff_mat[,a]))

perturb_mag <- c(0, 0, mofat_perturb, stand_perturb, strict_perturb, aug_stand_perturb, aug_strict_perturb)


############################# Q versus perturbation magnitude
method_names <- c(
  "MmLHD", "MaxPro", "MOFAT", "SOFT-Stand", "SOFT-Strict",
  "A-SOFT-Stand", "A-SOFT-Strict"
)
# Match the seven-design colors and symbols in aug_sf_prop_vs.R.
method_colors <- setNames(
  c("#5f5f5f", "#2e8b57", "#2f38ff", "#f4a261", "#cc45e0", "#b38b63", "#b980ff"),
  method_names
)
method_shapes <- setNames(c(18, 8, 16, 17, 15, 1, 2), method_names)
plot_data <- data.frame(
  method = factor(method_names, levels = method_names),
  Q_vals = Q_vals,
  log_Q_ratio = log(Q_vals / Q_vals[method_names == "MOFAT"]),
  perturb_mag = perturb_mag
)

# Points use the original values; only text is shifted for readability.
# Shift MmLHD left and MaxPro right, keeping both labels at the same height.
# Explicit label centers keep connectors attached to the middle of each name.
plot_data$label_x <- plot_data$log_Q_ratio + c(-0.4, 0.4, 0.2, -2, 2, 0, 0)
plot_data$label_y <- plot_data$perturb_mag + 0.30
soft_labels <- plot_data$method %in% c("SOFT-Stand", "SOFT-Strict")
plot_data$label_y[soft_labels] <- plot_data$label_y[soft_labels] + 0.10
aug_stand_label <- plot_data$method == "A-SOFT-Stand"
plot_data$label_y[aug_stand_label] <- plot_data$perturb_mag[aug_stand_label] - 0.30
# Connect to the label edge nearest its point, including labels below points.
plot_data$connector_y <- plot_data$label_y -
  ifelse(plot_data$label_y > plot_data$perturb_mag, 0.08, -0.08)
conclusion_plot <- ggplot2::ggplot(
  plot_data, ggplot2::aes(x = log_Q_ratio, y = perturb_mag, color = method)
) +
  ggplot2::geom_segment(
    ggplot2::aes(xend = label_x, yend = connector_y),
    color = "grey65", linewidth = 0.25, show.legend = FALSE
  ) +
  ggplot2::geom_point(ggplot2::aes(shape = method), size = 5, stroke = 1.2) +
  ggplot2::geom_text(
    ggplot2::aes(x = label_x, y = label_y, label = method),
    size = 6.5, fontface = "bold", hjust = 0.5,
    show.legend = FALSE
  ) +
  ggplot2::scale_x_continuous(
    limits = c(-1.5, max(plot_data$log_Q_ratio) + 3.5), expand = c(0, 0)
  ) +
  ggplot2::scale_y_continuous(
    limits = c(-0.25, max(perturb_mag) + 0.6), expand = c(0, 0)
  ) +
  ggplot2::scale_color_manual(values = method_colors) +
  ggplot2::scale_shape_manual(values = method_shapes) +
  ggplot2::labs(
    x = "Space-fillingness",
    y = "Screening capability"
  ) +
  ggplot2::theme_classic(base_size = 10) +
  ggplot2::theme(
    legend.position = "none",
    aspect.ratio = 0.9,
    panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.6),
    axis.line = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 22),
    axis.text = ggplot2::element_blank(),
    axis.ticks = ggplot2::element_blank()
  )

output_file <- file.path(
  "..", "..", "trade_off", "Q_vs_perturbation.png"
)
ggplot2::ggsave(output_file, conclusion_plot,
  width = 6, height = 6, units = "in", dpi = 600, bg = "white")
magick::image_write(
  magick::image_trim(magick::image_read(output_file)), output_file
)
