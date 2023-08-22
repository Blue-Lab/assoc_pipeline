#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Variant QC") %>%
    add_argument("gds_file", help="GDS file") %>%
    add_argument("--out_prefix", help="Prefix for output files", default="") %>%
    add_argument("--sample_id", help="File with vector of sample IDs") %>%
    add_argument("--pheno_file", help="Phenotype file in annotated dataframe format (used to account for sex in allele frequency calculations)") %>%
    add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1) %>%
    add_argument("--img_ext", help="File extension for plots", default = "png") %>%
    add_argument("--ms_bindwidth", help = "Bin width for missingness histogram", type = "numeric", default=0.005) %>%
    add_argument("--ms_xmin", help = "Variant-missingness plot x-axis lower limit", type = "numeric") %>%
    add_argument("--ms_xmax", help = "Variant-missingness plot x-axis upper limit", type = "numeric") %>%
    add_argument("--ms_ymin", help = "Variant-missingness plot y-axis lower limit", type = "numeric") %>%
    add_argument("--ms_ymax", help = "Variant-missingness plot y-axis upper limit", type = "numeric") %>%
    add_argument("--maf_bindwidth", help = "Bin width for MAF histogram", type = "numeric", default=0.01) %>%
    add_argument("--maf_xmin", help = "MAF plot x-axis lower limit", type = "numeric") %>%
    add_argument("--maf_xmax", help = "MAF plot x-axis upper limit", type = "numeric") %>%
    add_argument("--maf_ymin", help = "MAF plot y-axis lower limit", type = "numeric") %>%
    add_argument("--maf_ymax", help = "MAF plot y-axis upper limit", type = "numeric")
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
p_ms <- ggplot(var.df, aes(missing.rate)) +
    geom_histogram(binwidth=as.numeric(argv$ms_bindwidth), boundary=0) +
    # Need to coerce NAs to numeric because of ggplot2 bug:
    lims(x = c(as.numeric(argv$ms_xmin), as.numeric(argv$ms_xmax)),
         y = c(as.numeric(argv$ms_ymin), as.numeric(argv$ms_ymax)))
ggsave(paste0(argv$out_prefix, "missing_by_variant.", argv$img_ext), plot = p_ms)

p_maf <- ggplot(var.df, aes(maf)) +
    geom_histogram(binwidth=as.numeric(argv$maf_bindwidth), boundary=0) +
    lims(x = c(as.numeric(argv$maf_xmin), as.numeric(argv$maf_xmax)),
         y = c(as.numeric(argv$maf_ymin), as.numeric(argv$maf_ymax)))
ggsave(paste0(argv$out_prefix, "maf.", argv$img_ext), plot = p_maf)

seqClose(gds)
