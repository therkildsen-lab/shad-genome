library(psmcr)
library(pegas)



setwd("/local/workdir/azwad/shad-genome")

fAloSap1_diploid <- VCF2DNAbin("psmc/fAloSap1_minimap_dedup_readgroups_for_diploid_filtered.vcf.gz", refgenome = "assembly/GCF_018492685.1.fa", individual = 1, quiet = FALSE)


# Create a diploid consensus sequence from VCF and Ref genome fAloSap1_diploid
# <- VCF2DNAbin( 'reads/gatk_processed_variants_filtered.vcf.gz', refgenome =
# 'assembly/GCF_018492685.1.fa', individual = 1, quiet = FALSE )



# saveRDS(fAloSap1_diploid, 'psmc/fAloSap1_diploid_gatk.rds')

# fAloSap1_diploid <- readRDS('psmc/fAloSap1_diploid_gatk.rds')


fAloSap1_diploid_binned <- seqBinning(fAloSap1_diploid, bin.size = 100)

fAloSap1_psmc_test <- psmc(fAloSap1_diploid_binned,
    parapattern = "4+25*2+4+6",
    maxt = 15, niters = 25, trratio = 5, B = 0, trunksize = 5e5, decoding =
        FALSE, quiet = FALSE, raw.output = FALSE, mc.cores = 16
)


# saveRDS(fAloSap1_psmc, 'psmc/output_run_psmc_gatk.rds')

fAloSap1_psmc_gatk <- readRDS("/local/workdir/azwad/shad-genome/psmc/output_run_psmc_gatk.rds")

fAloSap1_psmc <- readRDS("/local/workdir/azwad/shad-genome/fAloSap1_minimap_dedup_readgroups_for_diploid.vcf_psmc.rds")

herring_psmc <- import.psmc("/local/workdir/azwad/shad-genome/other_genomes/ChASM_A3/reads/illumina/psmc/illumina_merged_filtered_diploid.psmc")



plot(fAloSap1_psmc_test,
    col = rgb(red = 0, green = 0, blue = 1, alpha = 0.5), mutation.rate = 2e-09,
    g = 4, xlim = c(10000, 1500000)
)

lines(herring_psmc,
    col = "red", mutation.rate = 2e-09,
    g = 6
)
