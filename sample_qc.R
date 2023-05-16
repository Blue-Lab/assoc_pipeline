#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Sample QC") %>%
  add_argument("gds_file", help="GDS file") %>%
  add_argument("--out_prefix", help="Prefix for output files",
               default="") %>%
  add_argument("--variant_id", help="File with vector of variant IDs") %>%
  add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1) %>%
  add_argument("--img_ext", help="File extension for plot", default = "png") %>%
  add_argument("--bindwidth", help = "Bin width for missingness histogram", type = "numeric") %>%
  add_argument("--xmin", help = "X-axis lower limit", type = "numeric") %>%
  add_argument("--xmax", help = "X-axis upper limit", type = "numeric") %>%
  add_argument("--ymin", help = "Y-axis lower limit", type = "numeric") %>%
  add_argument("--ymax", help = "Y-axis upper limit", type = "numeric")
argv <- parse_args(argp)

# load libraries
library(SeqArray)
library(ggplot2)

# log versions and arguments for reproducibility
sessionInfo()
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
p <- ggplot(miss.df, aes(missing.rate)) +
    geom_histogram(binwidth=as.numeric(bindwidth), boundary=0) +
    # Need to coerce NAs to numeric because of ggplot2 bug:
    lims(x = c(as.numeric(argv$xmin), as.numeric(argv$xmax)),
         y = c(as.numeric(argv$ymin), as.numeric(argv$ymax)))

ggsave(paste0(argv$out_prefix, "missing_by_sample.", argv$img_ext), plot=p)

seqClose(gds)
