library(tidyverse)
library(cowplot)

snps <-
  read_tsv("reads/fAloSap1_snps.bed",
           col_names = c("Chromosome", "Position0", "End", "Depth"))

cov_nogaps <- read_tsv(
  "reads/fAloSap1_pacbio_dedup_ucsc.coverage",
  col_names = c("Chromosome", "Position", "Depth")
) %>%
  filter(Chromosome %in% paste0("chr", c(1:24))) %>%
  mutate(Chromosome = as.integer(str_sub(Chromosome, 4, str_length(Chromosome))))

mean_cov  <- cov_nogaps %>%
  pull(Depth) %>% mean()

chroms <- paste0("chr", c(1:24))

cov <-
  read_tsv(
    "reads/fAloSap1_pacbio_coverage.ucsc.bedGraph",
    col_names = c("Chromosome", "Start", "End", "Depth")
  ) %>%
  filter(Chromosome %in% paste0("chr", c(1:24))) %>%
  mutate(Chromosome = as.integer(str_sub(Chromosome, 4, str_length(Chromosome))))

snps_filtered <- snps %>%
  filter(Chromosome %in% paste0("chr", c(1:24)),
         Depth > mean_cov / 3,
         Depth < mean_cov * 2) %>%
  mutate(Chromosome = as.integer(str_sub(Chromosome, 4, str_length(Chromosome))))

chrom_size <-
  cov %>%
  group_by(Chromosome) %>%
  summarize(chrom_length = max(End)) %>%
  select(Chromosome, chrom_length)

no_cov <- cov %>% filter(Depth == 0)

cov_binned <- cov_nogaps %>%
  mutate(ints = cut_width(Position, width = 5e4, boundary = 0)) %>%
  group_by(Chromosome, ints) %>%
  summarize(mean_depth = mean(Depth),
            int_start = min(Position)) %>%
  mutate(chrom_odd = ifelse(Chromosome %% 2 == 1, "odd", "even"))

genome_cov_plot <- cov_binned %>%
  ggplot() +
  geom_point(aes(x = int_start, y = mean_depth),
             size = 0.5) +
  ylim(c(0, 500)) +
  geom_rect(
    data = no_cov,
    aes(xmin = Start, xmax = End),
    ymin = 0,
    ymax = 1000,
    fill = "red",
    alpha = 1
  ) +
  facet_wrap(~ `Chromosome`,
             ncol = length(chroms),
             scales = "free_x") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.spacing.x = unit(0, "lines"),
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )



# Break SNPs into 1mb windows
snp_1mb_window <- snps_filtered %>%
  mutate(ints = cut_width(End, width = 1e6, boundary = 0)) %>%
  group_by(Chromosome, ints) %>%
  summarize(snp_per_kb = n() / (1e6 / 1e3),
            int_start = min(End)) %>%
  mutate(chrom_odd = ifelse(Chromosome %% 2 == 1, "odd", "even"))


## Snps per kb, 1mb windows
snp_window_plot <- snp_1mb_window %>%
  ggplot() +
  geom_col(
    aes(x = int_start, y = snp_per_kb, fill = chrom_odd),
    position = "jitter",
    width = 1000000
  ) +
  facet_wrap(
    ~ `Chromosome`,
    ncol = 24,
    strip.position = "top",
    scales = "free_x"
  ) +
  ggtitle("Heterozygosity by Chromosome (1mb windows) \n fAloSap1.pri") +
  xlab("") +
  ylab("Heterozygosity per kb") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.spacing.x = unit(0, "lines")
  ) +
  scale_fill_manual(values = c("dodgerblue4", "cornflowerblue"))





het_cov_plot <- plot_grid(snp_window_plot,
                          genome_cov_plot,
                          ncol = 1,
                          align = "v")

het_cov_plot

het_plot_density <- snps %>% ggplot() +
  geom_density(aes(x = `Position`, fill = as.character(Chromosome)), alpha = 0.5) +
  facet_wrap(~ `Chromosome`, ncol = 3, scales = "free_x")  +
  ggtitle("Heterozygosity by Chromosome \n fAloSap1.pri") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    plot.background = element_rect(fill = "white")
  )


snp_50kb_window <- snps_filtered %>%
  mutate(ints = cut_width(End, width = 5e4, boundary = 0)) %>%
  group_by(Chromosome, ints) %>%
  summarize(snp_per_kb = n() / (5e4 / 1e3),
            int_start = min(End)) %>%
  mutate(chrom_odd = ifelse(Chromosome %% 2 == 1, "odd", "even"))


snp_window_plot_50kb <- snp_50kb_window %>%
  ggplot() +
  geom_point(aes(x = int_start, y = snp_per_kb, color = chrom_odd),
             size = 0.5) +
  geom_hline(aes(yintercept = nrow(snps_filtered) / ((nrow(cov_nogaps) - sum(no_cov$End - no_cov$Start)) / 1000)),
             linetype = "dashed",
             color = "red") +
  facet_wrap(
    ~ `Chromosome`,
    ncol = 24,
    strip.position = "top",
    scales = "free_x"
  ) +
  ggtitle("Heterozygosity by Chromosome (50kb windows) \n fAloSap1.pri") +
  xlab("") +
  ylab("Heterozygosity per kb") +
  theme(
    legend.position = "none",
    axis.text.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.x = element_blank(),
    plot.background = element_rect(fill = "white"),
    panel.spacing.x = unit(0, "lines")
  ) +
  scale_color_manual(values = c("dodgerblue4", "cornflowerblue"))

het_cov_plot_50kb <- plot_grid(snp_window_plot_50kb,
                               genome_cov_plot,
                               ncol = 1,
                               align = "v")

ggsave("plots/het_cov_plot_50kb.png", 
       het_cov_plot_50kb,
       units = "in",
       dpi = "retina",
       width = 20,
       height = 10)



ggsave(
  plot = het_plot_density,
  filename = "plots/het_plot_density.png",
  dpi = "retina",
  width = 12,
  height = 15
)

ggsave(
  plot = het_plot_point,
  filename = "plots/het_plot_point.png",
  dpi = "retina",
  width = 20,
  height = 4
)

ggsave()

write_tsv(snps, "snp_density.txt")
