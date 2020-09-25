#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Variant QC") %>%
    add_argument("gds_file", help="GDS file") %>%
    add_argument("--out_prefix", help="Prefix for output files", default="") %>%
    add_argument("--sample_id", help="File with vector of sample IDs") %>%
    add_argument("--pheno_file", help="Phenotype file in annotated dataframe format (used to select only females for the X chromosome)") %>%
    add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1)
argv <- parse_args(argp)

# load libraries
library(SeqArray)
library(SeqVarTools)
library(Biobase)
library(ggplot2)

# log versions and arguments for reproducibility
sessionInfo()
print(argv)

# open GDS file
gds <- seqOpen(argv$gds_file)

# create SeqVarData
if (!is.na(argv$pheno_file)) {
    pheno <- readRDS(argv$pheno_file)
    gds <- SeqVarData(gds, sampleData = pheno)
}

# select samples
if (!is.na(argv$sample_id)) {
    sample.id <- readRDS(argv$sample_id)
    seqSetFilter(gds, sample.id=sample.id)
}

# autosomes
seqSetFilterChrom(gds, include=c(1:22, 25, 0)) # autosomes, PAR, unmapped
hw <- hwe(gds, parallel=argv$num_cores)
hw.perm <- hwe(gds, permute=TRUE, parallel=argv$num_cores)

# X chrom
if (is(gds, "SeqVarData") && "sex" %in% varLabels(sampleData(gds))) {
    seqSetFilterChrom(gds, include=c(23, "X"))
    female <- sampleData(gds)$sex %in% c("F", 2)
    seqSetFilter(gds, sample.sel=female, action="intersect")
    hw.x <- hwe(gds, parallel=argv$num_cores)
    hw.perm.x <- hwe(gds, permute=TRUE, parallel=argv$num_cores)

    hw <- rbind(hw, hw.x)
    hw.perm <- rbind(hw.perm, hw.perm.x)
}

saveRDS(hw, paste0(argv$out_prefix, "hwe.rds"))

p <- data.frame(obs=sort(hw$p),
                exp=sort(hw.perm$p)) %>%
    ggplot(aes(-log10(exp), -log10(obs))) +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    theme_bw()
ggsave(paste0(argv$out_prefix, "hwe.png"), plot=p, width=6, height=6)

seqClose(gds)
