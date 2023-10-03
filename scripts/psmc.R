library(psmcr)
library(pegas)
library(tidyverse)


setwd("/local/workdir/azwad/shad-genome")

# Create a diploid consensus sequence from VCF and Ref genome
fAloSap1_diploid <-
  VCF2DNAbin(
    "reads/gatk_processed_variants_filtered.vcf.gz",
    refgenome = "assembly/GCF_018492685.1.fa",
    individual = 1,
    quiet = FALSE
  )

saveRDS(fAloSap1_diploid, "psmc/fAloSap1_diploid_gatk.rds")

fAloSap1_diploid <- readRDS("psmc/fAloSap1_diploid_gatk.rds")


fAloSap1_diploid_binned <- seqBinning(fAloSap1_diploid, bin.size = 100)

fAloSap1_psmc <-
  psmc(
    fAloSap1_diploid_binned,
    parapattern = "4+25*2+4+6",
    maxt = 15,
    niters = 25,
    trratio = 5,
    B = 100,
    trunksize = 5e5,
    decoding = FALSE,
    quiet = FALSE,
    raw.output = FALSE,
    mc.cores = 16
  )


saveRDS(fAloSap1_psmc, "psmc/output_run_psmc_gatk.rds")

fAloSap1_psmc_gatk <- readRDS("psmc/output_run_psmc_gatk.rds")

fAloSap1_psmc <- readRDS("psmc/output_run_psmc.rds")

fAloSap1_psmc <- readRDS("psmc/output_run_psmc.rds")

plot(
  fAloSap1_psmc_gatk,
  col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
  mutation.rate = 2e-09,
  g = 3,
  bin.size = 100,
  xlim = c(0,100000)
)

glines(fAloSap1_psmc_gatk, col = "red")
