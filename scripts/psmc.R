library(psmcr)
library(pegas)
library(tidyverse)


setwd("/local/workdir/azwad/shad-genome")

# Create a diploid consensus sequence from VCF and Ref genome
fAloSap1_diploid <-
  VCF2DNAbin(
    "reads/freebayes_results.vcf.gz",
    refgenome = "assembly/GCF_018492685.1.fa",
    individual = 1,
    quiet = FALSE
  )

saveRDS(fAloSap1_diploid, "psmc/fAloSap1_diploid.rds")


fAloSap1_diploid_binned <-
  seqBinning(fAloSap1_diploid, bin.size = 100)

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
    mc.cores = 12
  )

fAloSap1_psmc %>% as.data.frame()


plot(
  fAloSap1_psmc,
  col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5),
  mutation.rate = 2e-09,
  g = 4.5,
  bin.size = 100,
  xlim = c(0,120000)
)

saveRDS(fAloSap1_psmc, "psmc/output_run_psmc.rds")

