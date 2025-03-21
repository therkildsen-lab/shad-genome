Synteny
================

``` r
# read in tab and convert to ucsc naming
synteny_tab_og <- read_delim("last/last_fAloSap1_medaka.tab",
    delim = "\t", escape_double = FALSE,
    col_names = FALSE, trim_ws = TRUE, skip = 31
)

synteny_tab <- read_delim("last/last_fAloSap1_medaka.tab",
    delim = "\t", escape_double = FALSE,
    col_names = FALSE, trim_ws = TRUE, skip = 31
) %>% relocate(X7, X8)

write_tsv(synteny_tab, "last/synteny_tab_mod.tab", col_names = F)
```

``` bash
reads/chromToUcsc -i last/synteny_tab_mod.tab -o last/synteny_tab_mod.ucsc.tab -a reads/fAloSap1.chromAlias.tsv
```

``` r
synteny_tab_ucsc_1 <- read_tsv("last/synteny_tab_mod.ucsc.tab", col_names = F) %>%
    relocate(X4, X5) %>%
    mutate(X1 = paste0("shad_", X1))

write_tsv(synteny_tab_ucsc_1,
    "last/synteny_tab_mod.tab",
    col_names = F
)
```

``` bash
reads/chromToUcsc -i last/synteny_tab_mod.tab -o last/synteny_tab_mod.ucsc.tab -a last/medaka.chromAlias.tsv
```

``` r
synteny_tab_final <- read_tsv("/workdir/azwad/shad-genome/last/synteny_tab_mod.ucsc.tab", col_names = F) %>%
    mutate(X1 = paste0("medaka_", X1))
```

``` r
synteny_tab_filtered <-
    synteny_tab_final %>%
    dplyr::rename(
        "medaka_chr" = "X1",
        "medaka_start" = "X2",
        "shad_chr" = "X3",
        "shad_start" = "X4",
        "length" = "X6"
    ) %>%
    select(c(medaka_chr, medaka_start, shad_chr, shad_start, length)) %>%
    mutate(medaka_end = medaka_start + length, .after = medaka_start) %>%
    mutate(shad_end = shad_start + length, .after = shad_start) %>%
    mutate(
        medaka_chr = str_remove(medaka_chr, "medaka_")
    ) %>%
    filter(
        shad_chr %in% c(paste0("shad_chr", c(1:24))),
        medaka_chr %in% c(paste0("chr", c(1:24)))
    )
```

``` r
shad_ref <-
    read_tsv("/workdir/azwad/shad-genome/assembly/sizes.genome.ucsc", col_names = F) %>%
    mutate(X1 = paste0("shad_", X1)) %>%
    filter(X1 %in% c(paste0("shad_chr", c(1:24))))

shad_cyto <-
    shad_ref %>%
    dplyr::rename(
        chr = X1,
        starts = X2,
        lengths = X3
    ) %>%
    mutate(cumsum = cumsum(lengths)) %>%
    select(1, 2, 4, 3)
shad_cyto$starts <- c(0, shad_cyto$cumsum[1:23])

medaka_ref <-
    read_tsv("/workdir/azwad/shad-genome/other_genomes/medaka/ref/medaka.ucsc.genome", col_names = F) %>% filter(X1 %in% c(paste0("chr", c(1:24))))

medaka_cyto <-
    medaka_ref %>%
    dplyr::rename(
        chr = X1,
        starts = X2,
        lengths = X3
    ) %>%
    mutate(cumsum = cumsum(lengths)) %>%
    select(1, 2, 4, 3)
medaka_cyto$starts <- c(0, medaka_cyto$cumsum[1:23])


cytoband <- rbind(shad_cyto, medaka_cyto)[, -4]

medaka_ref_rev <- medaka_ref %>% map_df(rev)

combined_ref <- rbind(medaka_ref_rev, shad_ref)
```

``` r
## break into two segments
# shad <- synteny_tab_corrected_filtered %>% select(shad_chr, shad_start, shad_end)
# medaka <- synteny_tab_corrected_filtered %>% select(medaka_chr, medaka_start, medaka_end)

shad <- synteny_tab_filtered %>%
    filter(length > 200, shad_chr == "shad_chr6") %>%
    select(shad_chr, shad_start, shad_end)
medaka <- synteny_tab_filtered %>%
    filter(length > 200, shad_chr == "shad_chr6") %>%
    select(medaka_chr, medaka_start, medaka_end)
```

``` r
chromcounts <- as.data.frame(table(shad$shad_chr))
col_palette <- turbo(n = nrow(chromcounts), alpha = 0.03)
col_list <- rep(col_palette, chromcounts$Freq)

circos.clear()
col_text <- "white"
circos.par("track.height" = 0.8, gap.degree = 0.5, cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = combined_ref$X1, xlim = matrix(c(rep(0, 48), combined_ref$X3), ncol = 2))

# genomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), gsub(".*chr", "", chr),
        cex = 1, col = col_text,
        facing = "bending.inside", niceFacing = TRUE
    )
}, bg.col = c(rep("gray70", 24), rep("cornflowerblue", 24)), bg.border = F, track.height = 0.12)


brk <- c(0, 10, 20, 30, 40, 50, 60) * 10^6
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(
        h = "top", major.at = brk, labels = round(brk / 10^7, 1), labels.cex = 0.5,
        col = "grey30", labels.col = "grey30", lwd = 1, labels.facing = "clockwise"
    )
}, bg.border = F)

# circos.track(
#     ylim = c(0,10),
#     cell.padding = c(0, 0, 0, 0),
#     track.height = mm_h(1),
#     bg.lwd= 30,
#     bg.col = c(rep("gray80", 24), rep("cornflowerblue", 24)),
#     bg.border = F
# )

circos.genomicLink(shad, medaka, col = col_list, border = col_list)
```

### Synteny with Allis shad (*Alosa alosa*)

``` r
# read in tab and convert to ucsc naming
synteny_tab_og <- read_delim("/workdir/azwad/shad-genome/last/last_fAloSap1_allis.tab",
    delim = "\t", escape_double = FALSE,
    col_names = FALSE, trim_ws = TRUE, skip = 31
)

synteny_tab <- read_delim("/workdir/azwad/shad-genome/last/last_fAloSap1_allis.tab",
    delim = "\t", escape_double = FALSE,
    col_names = FALSE, trim_ws = TRUE, skip = 31
) %>% relocate(X7, X8)

write_tsv(synteny_tab, "/workdir/azwad/shad-genome/last/allis_synteny_tab_mod.tab", col_names = F)
```

``` bash
/workdir/azwad/shad-genome/reads/chromToUcsc -i /workdir/azwad/shad-genome/last/allis_synteny_tab_mod.tab -o /workdir/azwad/shad-genome/last/allis_synteny_tab_mod.ucsc.tab -a /workdir/azwad/shad-genome/reads/fAloSap1.chromAlias.tsv
```

``` r
synteny_tab_ucsc_1 <- read_tsv("/workdir/azwad/shad-genome/last/allis_synteny_tab_mod.ucsc.tab", col_names = F) %>%
    relocate(X4, X5) %>%
    mutate(X1 = paste0("shad_", X1)) %>%
    filter(!X4 == "NA") %>%
    filter(!str_detect(X4, "J"))


write_tsv(synteny_tab_ucsc_1,
    "/workdir/azwad/shad-genome/last/allis_synteny_tab_mod.tab",
    col_names = F
)
```

``` bash
/workdir/azwad/shad-genome/reads/chromToUcsc -i /workdir/azwad/shad-genome/last/allis_synteny_tab_mod.tab -o /workdir/azwad/shad-genome/last/allis_synteny_tab_mod.ucsc.tab -a /workdir/azwad/shad-genome/last/GCF_017589495.1.chromAlias.txt
```

``` r
synteny_tab_final <- read_tsv("/workdir/azwad/shad-genome/last/allis_synteny_tab_mod.ucsc.tab", col_names = F) %>%
    mutate(X1 = paste0("allis_", X1))
```

``` r
synteny_tab_filtered <-
    synteny_tab_final %>%
    dplyr::rename(
        "allis_chr" = "X1",
        "allis_start" = "X2",
        "shad_chr" = "X3",
        "shad_start" = "X4",
        "length" = "X6"
    ) %>%
    select(c(allis_chr, allis_start, shad_chr, shad_start, length)) %>%
    mutate(allis_end = allis_start + length, .after = allis_start) %>%
    mutate(shad_end = shad_start + length, .after = shad_start) %>%
    mutate(
        allis_chr = str_remove(allis_chr, "allis_")
    ) %>%
    filter(
        shad_chr %in% c(paste0("shad_chr", c(1:24))),
        allis_chr %in% c(paste0("chr", c(1:24)))
    )
```

What percentage of alignments map to one or more chroms? (How high is
synteny between chroms)

``` r
# visualize pairwise relationships
synteny_tab_filtered %>% 
  group_by(shad_chr, allis_chr) %>% 
  summarize(count = n()) %>% 
  mutate(pct = count/sum(count) * 100) %>% 
  arrange(-pct) %>% 
  ggplot(aes(x = shad_chr, y = pct, fill = allis_chr, group = pct)) +
  geom_bar(position = "stack", stat = "identity")
```

![](synteny_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

``` r
# what is the mean pct of alignments (weighted mean)
synteny_tab_filtered %>% 
  group_by(shad_chr, allis_chr) %>% 
  summarize(count = n()) %>% 
  mutate(pct = count/sum(count) * 100, total = sum(count)) %>% 
  slice_max(pct, n = 1) %>% 
  ungroup() %>% 
  mutate(weight = total/sum(total)) %>% 
  summarize(mean_weighted = sum(pct * weight))
```

    ## # A tibble: 1 × 1
    ##   mean_weighted
    ##           <dbl>
    ## 1          83.9

``` r
# what is the mean pct of alignments (weighted mean)
synteny_tab_filtered %>% 
  group_by(shad_chr, allis_chr) %>% 
  summarize(count = n()) %>% 
  mutate(pct = count/sum(count) * 100, total = sum(count)) %>% 
  slice_max(pct, n = 1) %>% 
  ungroup() %>% 
  mutate(weight = total/sum(total)) %>% 
  summarize(mean_weighted = sum(pct * weight))
```

    ## # A tibble: 1 × 1
    ##   mean_weighted
    ##           <dbl>
    ## 1          83.9

``` r
# ID rearrangements
synteny_tab_filtered %>% 
  group_by(shad_chr, allis_chr) %>% 
  summarize(count = n()) %>% 
  mutate(pct = count/sum(count) * 100, total = sum(count)) %>% 
  ungroup() %>% 
  filter(pct < 50) %>% arrange(-pct)
```

    ## # A tibble: 552 × 5
    ##    shad_chr   allis_chr count   pct total
    ##    <chr>      <chr>     <int> <dbl> <int>
    ##  1 shad_chr1  chr17       639  2.49 25667
    ##  2 shad_chr4  chr4        452  2.44 18551
    ##  3 shad_chr4  chr16       435  2.34 18551
    ##  4 shad_chr2  chr1        449  1.95 23072
    ##  5 shad_chr24 chr6        322  1.93 16662
    ##  6 shad_chr16 chr22       319  1.89 16885
    ##  7 shad_chr14 chr11       288  1.81 15892
    ##  8 shad_chr3  chr22       365  1.81 20205
    ##  9 shad_chr3  chr7        363  1.80 20205
    ## 10 shad_chr2  chr23       398  1.73 23072
    ## # ℹ 542 more rows

``` r
# what about by length?
rearrangements_len <- synteny_tab_filtered %>% 
  group_by(shad_chr, allis_chr) %>% 
  summarize(count = n(), cum_len = sum(length)) %>% 
   mutate(pct = count/sum(count) * 100, total = sum(count)) %>% 
  ungroup() %>% 
  mutate(arr = ifelse(pct > 50, "syntenic", "rearrangement")) %>% 
  group_by(shad_chr, arr) %>% 
  summarize(total_len = sum(cum_len)) %>% 
  ungroup() %>% group_by(shad_chr) %>% mutate(pct_len = total_len/sum(total_len) * 100) 

rearrangements_len %>% filter(arr == "rearrangement") %>% pull(pct_len) %>% mean()
```

    ## [1] 5.234034

``` r
shad_ref <-
    read_tsv("/workdir/azwad/shad-genome/assembly/sizes.genome.ucsc", col_names = F) %>%
    mutate(X1 = paste0("shad_", X1)) %>%
    filter(X1 %in% c(paste0("shad_chr", c(1:24))))

shad_cyto <-
    shad_ref %>%
    dplyr::rename(
        chr = X1,
        starts = X2,
        lengths = X3
    ) %>%
    mutate(cumsum = cumsum(lengths)) %>%
    select(1, 2, 4, 3)
shad_cyto$starts <- c(0, shad_cyto$cumsum[1:23])

allis_ref <-
    read_tsv("/workdir/azwad/shad-genome/other_genomes/Alosa_alosa/assembly/sizes.genome.ucsc", col_names = F) %>% filter(X1 %in% c(paste0("chr", c(1:24))))

allis_cyto <-
    allis_ref %>%
    dplyr::rename(
        chr = X1,
        starts = X2,
        lengths = X3
    ) %>%
    mutate(cumsum = cumsum(lengths)) %>%
    select(1, 2, 4, 3)
allis_cyto$starts <- c(0, allis_cyto$cumsum[1:23])


cytoband <- rbind(shad_cyto, allis_cyto)[, -4]

allis_ref_rev <- allis_ref %>% map_df(rev)

combined_ref <- rbind(allis_ref_rev, shad_ref)

combined_ref %>% write_tsv("/local/workdir/azwad/shad-genome/markdowns/synteny_files/allis_american_shad_synteny_chroms.tsv")
```

``` r
## break into two segments

shad <- synteny_tab_filtered %>%
    filter(length > 500) %>%
    select(shad_chr, shad_start, shad_end)

palette <- tibble(color = turbo(n = 24, alpha = 0.01), shad_chr = paste0("shad_chr", c(1:24)))

shad_colors <- left_join(shad, palette) %>% select(color)


allis <- synteny_tab_filtered %>%
    filter(length > 500) %>%
    select(allis_chr, allis_start, allis_end)
```

``` r
combined_ref <- read_tsv("/local/workdir/azwad/shad-genome/markdowns/synteny_files/allis_american_shad_synteny_chroms.tsv")

chromcounts <- as.data.frame(table(shad$shad_chr))

circos.clear()
col_text <- "white"
circos.par("track.height" = 0.8, gap.degree = 0.5, cell.padding = c(0, 0, 0, 0))
circos.initialize(factors = combined_ref$X1, xlim = matrix(c(rep(0, 48), combined_ref$X3), ncol = 2))

# genomes
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr <- CELL_META$sector.index
    xlim <- CELL_META$xlim
    ylim <- CELL_META$ylim
    circos.text(mean(xlim), mean(ylim), gsub(".*chr", "", chr),
        cex = 1, col = col_text,
        facing = "bending.inside", niceFacing = TRUE
    )
}, bg.col = c(rep("gray70", 24), rep("cornflowerblue", 24)), bg.border = F, track.height = 0.12)


brk <- c(0, 10, 20, 30, 40, 50, 60) * 10^6
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.axis(
        h = "top", major.at = brk, labels = round(brk / 10^7, 1), labels.cex = 0.5,
        col = "grey30", labels.col = "grey30", lwd = 1, labels.facing = "clockwise"
    )
}, bg.border = F)

# circos.track(
#     ylim = c(0,10),
#     cell.padding = c(0, 0, 0, 0),
#     track.height = mm_h(1),
#     bg.lwd= 30,
#     bg.col = c(rep("gray80", 24), rep("cornflowerblue", 24)),
#     bg.border = F
# )

circos.genomicLink(shad, allis, col = shad_colors$color)
```

![](synteny_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

### Synteny across all chroms individually

``` r
circos_single_chrom <- function(chrom) {
    shad_chrom <- paste0("shad_chr", chrom)

    shad <- synteny_tab_filtered %>%
        filter(length > 200) %>%
        filter(shad_chr == shad_chrom) %>%
        select(shad_chr, shad_start, shad_end)

    allis <- synteny_tab_filtered %>%
        filter(length > 200) %>%
        filter(shad_chr == shad_chrom) %>%
        select(allis_chr, allis_start, allis_end)

    synteny_tab_chr <- synteny_tab_filtered %>%
        filter(length > 200) %>%
        filter(shad_chr == shad_chrom)

    comb_ref_chrom <- combined_ref %>%
        filter(X1 %in% c(synteny_tab_chr$shad_chr, synteny_tab_chr$allis_chr))
    chromcounts <- as.data.frame(table(shad$shad_chr))
    col_palette <- turbo(n = nrow(chromcounts), alpha = 0.03, begin = as.integer(gsub(".*chr", "", shad_chrom)) / 24)
    col_list <- rep(col_palette, chromcounts$Freq)

    circos.clear()
    col_text <- "white"
    circos.par("track.height" = 0.8, gap.degree = 0.5, cell.padding = c(0, 0, 0, 0))
    circos.initialize(factors = comb_ref_chrom$X1, xlim = matrix(c(rep(0, nrow(comb_ref_chrom)), comb_ref_chrom$X3), ncol = 2))

    # genomes
    circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
        chr <- CELL_META$sector.index
        xlim <- CELL_META$xlim
        ylim <- CELL_META$ylim
        circos.text(mean(xlim), mean(ylim), gsub(".*chr", "", chr),
            cex = 1, col = col_text,
            facing = "bending.inside", niceFacing = TRUE
        )
    }, bg.col = c(rep("gray70", (nrow(comb_ref_chrom) - 1)), rep("cornflowerblue", 1)), bg.border = F, track.height = 0.12)

    brk <- c(0, 10, 20, 30, 40, 50, 60) * 10^6
    circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
        circos.axis(
            h = "top", major.at = brk, labels = round(brk / 10^7, 1), labels.cex = 0.5,
            col = "grey30", labels.col = "grey30", lwd = 1, labels.facing = "clockwise"
        )
    }, bg.border = F)

    circos.genomicLink(shad, allis, col = col_list, border = col_list)
}
```

``` r
for (i in 1:24) {
    circos_single_chrom(i)
}
```

#### Linear Synteny plot with syntenyPlotteR

``` r
synteny_tab_linear <-
    read_tsv("/workdir/azwad/shad-genome/last/allis_synteny_tab_mod.ucsc.tab", col_names = F) %>%
    dplyr::rename(
        "allis_chr" = "X1",
        "allis_start" = "X2",
        "shad_chr" = "X3",
        "shad_start" = "X4",
        "orientation" = "X7",
        "length" = "X6"
    ) %>%
    select(c(allis_chr, allis_start, shad_chr, shad_start, length, orientation)) %>%
    mutate(allis_end = allis_start + length, .after = allis_start) %>%
    mutate(shad_end = shad_start + length, .after = shad_start) %>%
    mutate(
        shad_chr = str_remove(shad_chr, "shad_")
    ) %>%
    filter(
        shad_chr %in% c(paste0("chr", c(1:24))),
        allis_chr %in% c(paste0("chr", c(1:24)))
    ) %>%
    filter(length > 200) %>%
    select(-length) %>%
    relocate(c(allis_chr, allis_start, allis_end), .after = shad_end) %>%
    mutate(sps2 = "American shad", sps1 = "Allis shad") %>%
    mutate(allis_chr = as.integer(gsub(".*chr", "", allis_chr))) %>%
    mutate(shad_chr = as.integer(gsub(".*chr", "", shad_chr)))


shad_ref <-
    read_tsv("/workdir/azwad/shad-genome/assembly/sizes.genome.ucsc", col_names = F) %>%
    mutate(species = "American shad") %>%
    filter(X1 %in% c(paste0("chr", c(1:24))))


allis_ref <-
    read_tsv("/workdir/azwad/shad-genome/other_genomes/Alosa_alosa/assembly/sizes.genome.ucsc", col_names = F) %>%
    mutate(species = "Allis shad") %>%
    filter(X1 %in% c(paste0("chr", c(1:24))))


combined_ref_linear <- rbind(allis_ref, shad_ref) %>%
    select(-X2) %>%
    mutate(X1 = as.integer(gsub(".*chr", "", X1)))

write_tsv(combined_ref_linear, "/workdir/azwad/shad-genome/last/allis_american_sizes.tsv", col_names = F)
write_tsv(synteny_tab_linear, "/workdir/azwad/shad-genome/last/allis_american_synteny.tsv", col_names = F)
```

``` r
draw.linear(
    "linear_synteny_allis",
    w = 25,
    h = 4,
    colours = c(turbo(n = 24)),
    "/workdir/azwad/shad-genome/last/allis_american_sizes.tsv",
    "/workdir/azwad/shad-genome/last/allis_american_synteny.tsv",
    directory = "/workdir/azwad/shad-genome/markdowns/synteny_files"
)
```
