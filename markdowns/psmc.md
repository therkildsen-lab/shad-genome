PSMC
================

### Plot PSMC results for Shad & Herring

``` r
setwd("/local/workdir/azwad/shad-genome/psmc")

# png(filename="/workdir/azwad/shad-genome/markdowns/psmc_files/figure-gfm/test.png", height=8, width = 10, units = "in", type = "cairo", res = 320)

# Set up plot
plot(c(10000, 40000000), c(0, 110),
    pch = "",
    log = "x", ylab = "", xlab = "", xaxt = "n", type = "n",
    cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25
)
axis(1,
    tcl = -0.7,
    at = c(10^(4:10)),
    labels = sapply(c(4:10), function(i) as.expression(bquote(10^.(i)))),
    cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25
)
axis(1,
    tcl = -0.4, labels = NA, cex.lab=1.25, cex.axis=1.25, cex.main=1.25, cex.sub=1.25,
    at = c(
        2 * 10^(4:10), 3 * 10^(4:10), 4 * 10^(4:10), 5 * 10^(4:10),
        6 * 10^(4:10), 7 * 10^(4:10), 8 * 10^(4:10), 9 * 10^(4:10)
    )
)
mtext("Years", side = 1, line = 2, cex=1.25)
mtext(expression(paste("Effective Population Size (x10"^"4", ")")),
    side = 2, line = 2, cex=1.25,
)



## Plot shad data

# Plot bootstrapped lines
for (i in 0:100) {
    psmc <- read_table(paste0("plot_data/shad_psmc_combined.", i, ".txt"), col_names = F)
    lines(psmc$X1, psmc$X2, type = "s", lwd = 2, col = "#FF000008")
}

# Plot initial PSMC line
psmc <- read_table("plot_data/shad_psmc_combined.0.txt", col_names = F)
lines(psmc$X1, psmc$X2, type = "s", lwd = 4, col = "#880000")



# Plot in herring data

setwd("/local/workdir/azwad/shad-genome/other_genomes/herring_2015/psmc_no_A_flag_mpileup")

# Plot bootstrapped lines
for (i in 0:96) {
    psmc <- read_table(paste0("plot_data/herring_psmc_combined.", i, ".txt"), col_names = F)
    lines(psmc$X1, psmc$X2, type = "s", lwd = 2, col = "#0050910e")
}

psmc <- read_table("plot_data/herring_psmc_combined.0.txt", col_names = F)
lines(psmc$X1, psmc$X2, type = "s", lwd = 4, col = "#002c50")


# Plot in Allis shad data

setwd("/local/workdir/azwad/shad-genome/other_genomes/Alosa_alosa/reads/psmc")

# Plot bootstrapped lines
for (i in 0:100) {
    psmc <- read_table(paste0("plot_data/allis_psmc_combined.", i, ".txt"), col_names = F) %>%
        filter(X2 < 400)
    lines(psmc$X1, psmc$X2, type = "s", lwd = 2, col = "#0088560c")
}

psmc <- read_table("plot_data/allis_psmc_combined.0.txt", col_names = F) %>%
    filter(X2 < 400)

lines(psmc$X1, psmc$X2, type = "s", lwd = 4, col = "#008856")



legend("topright",
    legend = c("Atlantic herring", "American shad", "Allis shad"),
    col = c("#002c50", "#880000", "#008856"),
    lty = 1,
    lwd = 4,
    cex = 1.3
)


rect(1e4, -1000, 115000, 1e6, col="#FF7F0E1A", border = NA)
```

![](psmc_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
# dev.off()
```
