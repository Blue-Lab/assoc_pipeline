#! /usr/bin/env Rscript
library(SeqArray)
library(SeqVarTools)
library(ggplot2)
library(argparser)
library(magrittr)
sessionInfo()

# read arguments
argp <- arg_parser("Sample QC") %>%
  add_argument("gds_file", help="GDS file") %>%
  add_argument("--out_prefix", help="Prefix for output files",
               default="") %>%
  add_argument("--variant_id", help="File with vector of variant IDs") %>%
  add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1)

argv <- parse_args(argp)
print(argv)

# open GDS file
gds <- seqOpen(argv$gds_file)

# select variants
if (!is.na(argv$variant_id)) {
    variant.id <- readRDS(argv$variant_id)
    seqSetFilter(gds, variant.id=variant.id)
}

# missing call rate
# the SeqVarTools missingGenotypeRate is merely a wrapper around seqMissing
missing.rate <- seqMissing(gds, per.variant=FALSE, parallel=argv$num_cores)
sample.id <- seqGetData(gds, "sample.id")

miss.df <- data.frame(sample.id, missing.rate, stringsAsFactors=FALSE)
saveRDS(miss.df, paste0(argv$out_prefix, "missing_by_sample.rds"))

# plot
ggplot(miss.df, aes(missing.rate)) +
    geom_histogram(binwidth=0.01, boundary=0)
ggsave(paste0(argv$out_prefix, "missing_by_sample.pdf"))

seqClose(gds)
