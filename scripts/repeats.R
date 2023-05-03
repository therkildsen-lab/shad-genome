library(tidyverse)
library(viridis)

shad_repeats <-
  read_tsv(
    "repeats/genome_mask/fAloSap1_repeats.out",
    col_names = c("chrom", "start", "end", "type")
  ) %>% mutate(type = gsub("\\/.*","", type))

chroms <- paste0("chr", c(1:24))

chrom_sizes <-
  read_tsv("assembly/sizes.genome.ucsc",
           col_names = c("chrom", "start", "length")) %>%
  select(c(chrom, length)) %>%
  filter(chrom %in% chroms) %>%
  mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom))))

repeat_plot <- shad_repeats %>% filter(chrom %in% chroms) %>%
  mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom)))) %>%
  ggplot()  +
  geom_rect(
    data = chrom_sizes,
    aes(ymin = 0, ymax = length),
    xmin = 0,
    xmax = 1,
    col = NA,
    fill = "white"
  ) +
  geom_rect(
    aes(ymin = start, ymax = end, fill = type),
    xmin = 0,
    xmax = 1,
    alpha = 1
  ) +
  geom_rect(
    data = chrom_sizes,
    aes(ymin = 0, ymax = length),
    xmin = 0,
    xmax = 1,
    col = "black" ,
    fill = NA
  ) +
  scale_fill_viridis(discrete = TRUE, direction = -1) +
  theme_minimal() +
  theme(
        strip.text.y.left = element_text(angle = 0),
        plot.background = element_rect(fill = "gray80")) +
  coord_flip() +
  facet_wrap( ~ `chrom`, nrow = length(chroms), strip.position = "left")


repeat_plot 
