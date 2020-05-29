#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Variant QC") %>%
    add_argument("gds_file", help="GDS file") %>%
    add_argument("--out_prefix", help="Prefix for output files", default="") %>%
    add_argument("--sample_id", help="File with vector of sample IDs") %>%
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

# select samples
if (!is.na(argv$sample_id)) {
    sample.id <- readRDS(argv$sample_id)
    seqSetFilter(gds, sample.id=sample.id)
}

hw <- hwe(gds, parallel=argv$num_cores)
saveRDS(hw, paste0(argv$out_prefix, "hwe.rds"))

hw.perm <- hwe(gds, permute=TRUE, parallel=argv$num_cores)



p <- data.frame(obs=sort(hw$p),
                exp=sort(hw.perm$p)) %>%
    ggplot(aes(-log10(exp), -log10(obs))) +
    geom_point() +
    geom_abline(intercept=0, slope=1, color="red") +
    xlab(expression(paste(-log[10], "(expected P)"))) +
    ylab(expression(paste(-log[10], "(observed P)"))) +
    theme_bw()
ggsave(paste0(argv$out_prefix, "hwe.png"), plot=p, width=6, height=6)
