library(tidyverse)
library(viridis)
library(reshape)
library(hrbrthemes)
library(gridExtra)
library(ggplot2)

shad_repeats <- read_tsv("repeats/genome_mask/fAloSap1_repeats.out", col_names = c("chrom",
    "start", "end", "type")) %>%
    mutate(type = gsub("\\/.*", "", type))

chroms <- paste0("chr", c(1:24))

chrom_sizes <- read_tsv("assembly/sizes.genome.ucsc", col_names = c("chrom", "start",
    "length")) %>%
    select(c(chrom, length)) %>%
    filter(chrom %in% chroms) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom))))

repeat_plot <- shad_repeats %>%
    filter(chrom %in% chroms) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom)))) %>%
    ggplot() + geom_rect(data = chrom_sizes, aes(ymin = 0, ymax = length), xmin = 0,
    xmax = 1, col = NA, fill = "gray90") + geom_rect(aes(ymin = start, ymax = end,
    fill = type), xmin = 0, xmax = 1, alpha = 1) + geom_rect(data = chrom_sizes,
    aes(ymin = 0, ymax = length), xmin = 0, xmax = 1, col = "black", fill = NA) +
    scale_fill_viridis(discrete = TRUE, direction = -1) + theme_minimal() + theme(strip.text.y.left = element_text(angle = 0),
    plot.background = element_rect(fill = "white")) + coord_flip() + facet_wrap(~chrom,
    nrow = length(chroms), strip.position = "left")


repeat_plot



## Code from Kristina Galagova; See: https://github.com/oushujun/EDTA/issues/92

KimuraDistance <- read.csv("repeats/genome_mask_2/GCF_018492685.fa.distance", sep = " ")

# add here the genome size in bp
genomes_size = 903581644

kd_melt = melt(KimuraDistance, id = "Div")
kd_melt$norm = kd_melt$value/genomes_size * 100

ggplot(kd_melt, aes(fill = variable, y = norm, x = Div)) + geom_bar(position = "stack",
    stat = "identity", color = "black") + scale_fill_viridis(discrete = T) + theme_classic() +
    xlab("Kimura substitution level") + ylab("Percent of the genome") + labs(fill = "") +
    coord_cartesian(xlim = c(0, 55)) + theme(axis.text = element_text(size = 11),
    axis.title = element_text(size = 12))


