#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("LD pruning") %>%
    add_argument("gds_file", help="GDS file") %>%
    add_argument("--chromosome", help="chromosome") %>%
    add_argument("--out_prefix", help="Prefix for output files", default="") %>%
    add_argument("--sample_id", help="RDS file with vector of sample.id to include") %>%
    add_argument("--variant_id", help="RDS file with vector of variant.id to include") %>%
    add_argument("--maf", help="minimum MAF for variants to include", default=0.05) %>%
    add_argument("--missing", help="maximum missing call rate for variants to include", default=0.05) %>%
    add_argument("--r_threshold", help="r threshold for LD", default=sqrt(0.1)) %>%
    add_argument("--window_size", help="window size in Mb", default=10) %>%
    add_argument("--exclude_pca_corr_from_build", help="exclude variants in regions with high correlation with PCs (HLA, LCT, inversions), using positions from this genome build (hg18, hg19, or hg38)")
argv <- parse_args(argp)

library(SeqArray)
library(SeqVarTools)
library(SNPRelate)
library(dplyr)
sessionInfo()
print(argv)

# parse file paths
sample.id <- if (!is.na(argv$sample_id)) readRDS(argv$sample_id) else NULL
variant.id <- if (!is.na(argv$variant_id)) readRDS(argv$variant_id) else NULL

# open GDS file
gds <- seqOpen(argv$gds_file)

# select chromosome
if (!is.na(argv$chromosome)) {
    seqSetFilterChrom(gds, include=argv$chromosome)
}

var.info <- variantInfo(gds, alleles=FALSE)
if (!is.na(argv$variant_id)) {
    var.incl <- readRDS(argv$variant_id)
    var.info <- filter(var.info, variant.id %in% var.incl)
}

# filter PCA correlated regions
if (!is.na(argv$exclude_pca_corr_from_build)) {
    build <- argv$exclude_pca_corr_from_build
    filt <- get(data(list=paste("pcaSnpFilters", build, sep="."), package="GWASTools"))
    pca.filt <- rep(TRUE, nrow(var.info))
    for (f in 1:nrow(filt)) {
        pca.filt[var.info$chr == filt$chrom[f] &
                 filt$start.base[f] < var.info$pos &
                 var.info$pos < filt$end.base[f]] <- FALSE
    }
    var.info <- var.info[pca.filt,]
}

# run LD pruning
snpset <- snpgdsLDpruning(gds,
                          sample.id = sample.id,
                          snp.id = var.info$variant.id,
                          maf = argv$maf,
                          missing.rate = argv$missing,
                          method = "corr",
                          slide.max.bp = argv$window_size * 1e6,
                          ld.threshold = argv$r_threshold)

# convert list with one element per chrom to vector
pruned <- unlist(snpset, use.names=FALSE)

if (!is.na(argv$chromosome)) {
    out_file <- paste0(argv$out_prefix, "pruned_snps_chr", argv$chromosome, ".rds")
} else {
    out_file <- paste0(argv$out_prefix, "pruned_snps.rds")
}
saveRDS(pruned, file=out_file)

seqClose(gds)
