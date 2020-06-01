#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Variant QC") %>%
    add_argument("gds_file", help="GDS file") %>%
    add_argument("--out_prefix", help="Prefix for output files", default="") %>%
    add_argument("--sample_id", help="File with vector of sample IDs") %>%
    add_argument("--pheno_file", help="Phenotype file in annotated dataframe format (used to account for sex in allele frequency calculations)") %>%
    add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1)
argv <- parse_args(argp)

# load libraries
library(SeqArray)
library(SeqVarTools)
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

# missing call rate
# the SeqVarTools missingGenotypeRate is merely a wrapper around seqMissing
missing.rate <- seqMissing(gds, per.variant=TRUE, parallel=argv$num_cores)

# allele frequency
ref.freq <- alleleFrequency(gds, parallel=argv$num_cores)
maf <- pmin(ref.freq, 1-ref.freq)

variant.id <- seqGetData(gds, "variant.id")
var.df <- data.frame(variant.id, missing.rate, ref.freq, maf)
saveRDS(var.df, paste0(argv$out_prefix, "variant_metrics.rds"))

# plot
p <- ggplot(var.df, aes(missing.rate)) +
    geom_histogram(binwidth=0.01, boundary=0)
ggsave(paste0(argv$out_prefix, "missing_by_variant.pdf"), plot=p)

p <- ggplot(var.df, aes(maf)) +
    geom_histogram(binwidth=0.01, boundary=0)
ggsave(paste0(argv$out_prefix, "maf.pdf"), plot=p)

seqClose(gds)
