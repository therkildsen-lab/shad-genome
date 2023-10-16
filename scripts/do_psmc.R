library(psmcr)
library(pegas)

args <- commandArgs(trailingOnly = TRUE)
vcf <- as.character(args[1])
ref <- as.character(args[2])

# Create a diploid consensus DNAbin object (equivalent to vcf2fq in the original PSMC)
sample_diploid <- VCF2DNAbin(
    vcf,
    refgenome = ref,
    individual = 1,
    quiet = FALSE
)

# Bin the consensus sequence
sample_diploid_binned <- seqBinning(
    sample_diploid,
    bin.size = 100
)

# Run PSMC with settings recommended for humans, 100 bootstraps, 16 threads for bootstrapping
sample_psmc <- psmc(
    sample_diploid_binned,
    parapattern = "4+25*2+4+6",
    maxt = 15,
    niters = 25,
    trratio = 5,
    B = 100,
    trunksize = 5e+05,
    decoding = FALSE,
    quiet = FALSE,
    raw.output = FALSE,
    mc.cores = 16
)

# Save the PSMC object for plotting later
saveRDS(
    sample_psmc,
    paste0(tools::file_path_sans_ext(basename(vcf)), "_psmc.rds")
)
