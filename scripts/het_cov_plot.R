library(tidyverse)
library(cowplot)

het_cov_plot <- plot_grid(snp_window_plot, genome_cov_plot, ncol = 1, align = "v")
