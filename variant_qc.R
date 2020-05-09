#! /usr/bin/env Rscript
library(SeqArray)
library(SeqVarTools)
library(ggplot2)
library(argparser)
library(magrittr)
sessionInfo()

# read arguments
argp <- arg_parser("Variant QC") %>%
  add_argument("gds_file", help="GDS file") %>%
  add_argument("--out_prefix", help="Prefix for output files",
               default="") %>%
  add_argument("--sample_id", help="File with vector of sample IDs") %>%
  add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1)

argv <- parse_args(argp)
print(argv)

# open GDS file
gds <- seqOpen(argv$gds_file)

# select samples
if (!is.na(argv$sample_id)) {
    sample.id <- readRDS(argv$sample_id)
    seqSetFilter(gds, sample.id=sample.id)
}

# missing call rate
# the SeqVarTools missingGenotypeRate is merely a wrapper around seqMissing
missing.rate <- seqMissing(gds, per.variant=TRUE, parallel=argv$num_cores)
variant.id <- seqGetData(gds, "variant.id")

miss.df <- data.frame(variant.id, missing.rate, stringsAsFactors=FALSE)
saveRDS(miss.df, paste0(argv$out_prefix, "missing_by_variant.rds"))

# plot
ggplot(miss.df, aes(missing.rate)) +
    geom_histogram(binwidth=0.01, boundary=0)
ggsave(paste0(argv$out_prefix, "missing_by_variant.pdf"))

seqClose(gds)
