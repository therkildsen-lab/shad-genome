Repetitive Regions
================

# Annotating Repetitive Regions

Pacbio HiFi reads used to generate the reference sequence for
fAloSap1.pri were mapped back onto the reference genome with [minimap2
2.24](https://github.com/lh3/minimap2) with the `map-hifi` setting and
sorted using `samtools sort`. To annotate the repetitive regions in the
genome, [RepeatModeler
(v2.0.1)](http://www.repeatmasker.org/RepeatModeler/) was used to build
a repeat library for the American shad genome employing the
[Repbase](https://www.girinst.org/repbase/) database using the latest
version 28.04 (04/27/2023). The repeat consensus database was then
filtered for known protein sequences with uniref90 (known transposases
provided with Repeatmasker were pre-removed from the uniref90 database).
The remaining repeat database was then classified with Repbase, and used
to generate a repeat annotation gff file using RepeatMasker (Methods
adapted from [Stanhope et
al. 2022](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9826928/)) and
provided by the Cornell BioHPC User Guide for
[RepeatModeler](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=259#c)
and
[RepeatMasker](https://biohpc.cornell.edu/lab/userguide.aspx?a=software&i=62#c).

The repeat sequences in the genome assembly was soft masked with
RepeatMasker (version 4.1.0), using the repeat library created by
combining the Repbase database (latest version 10/26/2018) with de novo
TEs identified by RepeatModeler.

#### Summary of repetitive regions in the assembly

#### Genome-wide distribution of repetitive regions

``` r
setwd("/workdir/azwad/shad-genome")

shad_repeats <-
    read_tsv(
        "repeats/genome_mask/fAloSap1_repeats.out",
        col_names = c("chrom", "start", "end", "type")
    ) %>% mutate(type = gsub("\\/.*", "", type), type = ifelse(type == "SINE?", "SINE", type))

chroms <- paste0("chr", c(1:24))

chrom_sizes <-
    read_tsv("assembly/sizes.genome.ucsc",
        col_names = c("chrom", "start", "length")
    ) %>%
    select(c(chrom, length)) %>%
    filter(chrom %in% chroms) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom))))

repeat_plot <- shad_repeats %>%
    filter(chrom %in% chroms) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom)))) %>%
    dplyr::rename(Class = type) %>%
    ggplot() +
    geom_rect(
        data = chrom_sizes,
        aes(ymin = 0, ymax = length),
        xmin = 0,
        xmax = 1,
        col = NA,
        fill = "gray90"
    ) +
    geom_rect(
        aes(ymin = start, ymax = end, fill = Class),
        xmin = 0,
        xmax = 1,
        alpha = 1
    ) +
    geom_rect(
        data = chrom_sizes,
        aes(ymin = 0, ymax = length),
        xmin = 0,
        xmax = 1,
        col = "black",
        fill = NA
    ) +
    ylab("Position") +
    scale_fill_viridis(discrete = TRUE) +
    theme_minimal() +
    theme(
        strip.text.y.left = element_text(angle = 0),
        plot.background = element_rect(fill = "white", colour = "#ffffff00"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 14),
        strip.text = element_text(size = 16)
    ) +
    coord_flip() +
    facet_wrap(~`chrom`, nrow = length(chroms), strip.position = "left")


repeat_plot
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->
<br>

#### Separating repeats by type to identify any spacial patterns of distribution on the genome

``` r
shad_repeats %>%
    filter(chrom %in% chroms) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom)))) %>%
    left_join(chrom_sizes) %>%
    mutate(relative_pos = ((end + start) / 2) / length) %>%
    ggplot() +
    geom_density(aes(x = relative_pos, fill = type), alpha = 0.5) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap(~type, scales = "free_y") +
    ylab("Density") +
    xlab("Relative Genomic Position") +
    theme_classic() +
    theme(legend.position = "none")
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
shad_repeats %>%
    filter(chrom %in% chroms) %>%
    mutate(chrom = as.integer(str_sub(chrom, 4, str_length(chrom)))) %>%
    ggplot() +
    geom_histogram(aes(x = start, fill = type), alpha = 0.5) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap(~chrom, scales = "free", nrow = length(chroms), strip.position = "right") +
    ylab("Density") +
    xlab("Genomic Position") +
    theme_classic()
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->

#### Composition of Repetitive Elements

``` r
## Code from Kristina Galagova; See: https://github.com/oushujun/EDTA/issues/92
setwd("/workdir/azwad/shad-genome")
KimuraDistance <- read.csv("repeats/genome_mask_2/GCF_018492685.fa.distance", sep = " ")

# add here the genome size in bp
genomes_size <- 903581644

kd_melt <- melt(KimuraDistance, id = "Div") %>% filter(variable != "X")
kd_melt$norm <- kd_melt$value / genomes_size * 100

ggplot(kd_melt, aes(fill = variable, y = norm, x = Div)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    scale_fill_viridis(discrete = T) +
    theme_classic() +
    xlab("Kimura substitution level") +
    ylab("Percent of the genome (%)") +
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 55)) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
# Simplified by general class

kd_melt %>%
    mutate(variable = gsub("\\..*", "", variable), variable = ifelse(variable == "SINE?", "SINE", variable)) %>%
    ggplot(aes(fill = variable, y = norm, x = Div)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    scale_fill_viridis_d() +
    theme_classic() +
    xlab("Kimura substitution level") +
    ylab("Percent of the genome (%)") +
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 55)) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

``` r
### Just DNA elements
kd_melt %>%
    mutate(variable = gsub("\\..*", "", variable), variable = ifelse(variable == "SINE?", "SINE", variable)) %>%
    filter(variable == "DNA") %>%
    ggplot(aes(fill = variable, y = norm, x = Div)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    scale_fill_viridis_d() +
    theme_classic() +
    xlab("Kimura substitution level") +
    ylab("Percent of the genome (%)") +
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 55)) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-3-3.png)<!-- -->

``` r
### Just "Unknown" elements
kd_melt %>%
    mutate(variable = gsub("\\..*", "", variable), variable = ifelse(variable == "SINE?", "SINE", variable)) %>%
    filter(variable == "Unknown") %>%
    ggplot(aes(fill = variable, y = norm, x = Div)) +
    geom_bar(position = "stack", stat = "identity", color = "black") +
    scale_fill_viridis_d() +
    theme_classic() +
    xlab("Kimura substitution level") +
    ylab("Percent of the genome (%)") +
    labs(fill = "") +
    coord_cartesian(xlim = c(0, 55)) +
    theme(axis.text = element_text(size = 16), axis.title = element_text(size = 18))
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-3-4.png)<!-- -->

#### Comparing composition of repetitive elements between Atlantic herring, Allis shad, and American shad

``` r
setwd("/workdir/azwad/shad-genome")

shad <- read_tsv("repeats/genome_mask_2/GCF_018492685.fa.div")
allis <- read_tsv("other_genomes/Alosa_alosa/repeats/genome_mask_final_2/GCA_017589495.2_AALO_Geno_1.1_genomic.fna.div")
herring <- read_tsv("other_genomes/Ch_v2.0.2/repeats/genome_mask_final_2/GCF_900700415.2_Ch_v2.0.2_genomic.fna.div")

# add here the genome size in bp
genome_size_shad <- 903581644
genome_size_allis <- 854464681
genome_size_herring <- 725686888

shad_summary <- shad %>%
    mutate(Class = gsub("\\/.*", "", Class), Class = ifelse(Class == "SINE?", "SINE", Class)) %>%
    group_by(Class) %>%
    summarize(
        sumlen = sum(absLen),
        len_percent = sumlen / genome_size_shad * 100
    )

shad_summary <- shad_summary %>%
    add_row(Class = "Unmasked", sumlen = genome_size_shad - sum(shad_summary$sumlen), len_percent = 100 - sum(shad_summary$len_percent)) %>%
    mutate(species = "Alosa sapidissima")

allis_summary <- allis %>%
    mutate(Class = gsub("\\/.*", "", Class), Class = ifelse(Class == "SINE?", "SINE", Class)) %>%
    group_by(Class) %>%
    summarize(
        sumlen = sum(absLen),
        len_percent = sumlen / genome_size_allis * 100
    )

allis_summary <- allis_summary %>%
    add_row(Class = "Unmasked", sumlen = genome_size_allis - sum(allis_summary$sumlen), len_percent = 100 - sum(allis_summary$len_percent)) %>%
    mutate(species = "Alosa alosa")

herring_summary <- herring %>%
    mutate(Class = gsub("\\/.*", "", Class), Class = ifelse(Class == "SINE?", "SINE", Class)) %>%
    group_by(Class) %>%
    summarize(
        sumlen = sum(absLen),
        len_percent = sumlen / genome_size_herring * 100
    )

herring_summary <- herring_summary %>%
    add_row(Class = "Unmasked", sumlen = genome_size_herring - sum(herring_summary$sumlen), len_percent = 100 - sum(herring_summary$len_percent)) %>%
    mutate(species = "Clupea harengus")


total_summary <- rbind(shad_summary, allis_summary, herring_summary)

mycolor <- c(viridis(6, option = "D"), "grey60")

total_summary %>%
    filter(len_percent > 0.7) %>%
    mutate(species = factor(species, levels = c("Clupea harengus", "Alosa alosa", "Alosa sapidissima"))) %>%
    ggplot(aes(x = species, y = len_percent, fill = Class)) +
    geom_col(width = 0.5, position = position_fill()) +
    geom_text(aes(label = round(len_percent, digits = 1)),
        color = "white",
        position = position_fill(vjust = 0.5),
        size = 4
    ) +
    coord_flip() +
    scale_fill_manual(values = mycolor) +
    theme_classic() +
    xlab("") +
    ylab("Percentage of sequence (%)") +
    scale_y_continuous(labels = c(0, 25, 50, 75, 100)) +
    theme(
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 16, colour = "black", face = "italic"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)
    )
```

![](repetitive_regions_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->
