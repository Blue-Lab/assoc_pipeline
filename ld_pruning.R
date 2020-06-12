#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("LD pruning") %>%
    add_argument("gds_file", help="GDS file") %>%
    add_argument("--out_file", help="output file name", default="pruned_snps.rds") %>%
    add_argument("--sample_id", help="RDS file with vector of sample.id to include") %>%
    add_argument("--variant_id", help="RDS file with vector of variant.id to include") %>%
    add_argument("--maf", help="minimum MAF for variants to include", default=0.05) %>%
    add_argument("--missing", help="maximum missing call rate for variants to include", default=0.05) %>%
    add_argument("--r_threshold", help="r threshold for LD", default=sqrt(0.1)) %>%
    add_argument("--window_size", help="window size in Mb", default=10)
argv <- parse_args(argp)

library(SeqArray)
library(SNPRelate)
sessionInfo()
print(argv)

# parse file paths
gds.file <- argv$gds_file
out.file <- argv$out_file
sample.id <- if (!is.na(argv$sample_id)) readRDS(argv$sample_id) else NULL
variant.id <- if (!is.na(argv$variant_id)) readRDS(argv$variant_id) else NULL

# open GDS file
gds <- seqOpen(gds.file)

# run LD pruning
snpset <- snpgdsLDpruning(gds,
                          sample.id = sample.id,
                          snp.id = variant.id,
                          maf = argv$maf,
                          missing.rate = argv$missing,
                          method = "corr",
                          slide.max.bp = argv$window_size * 1e6, 
                          ld.threshold = argv$r_threshold)

# convert list with one element per chrom to vector
pruned <- unlist(snpset, use.names=FALSE)

saveRDS(pruned, file=out.file)
