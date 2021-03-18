#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("coxmeg") %>%
    add_argument("--gds_file", help = "GDS file") %>%
    add_argument("--pheno_file", help = "Phenotype file in annotated dataframe format") %>%
    add_argument("--grm_file", help = "GRM file") %>%
    add_argument("--grm_type", help = "matrix type: 'bd' (block diagonal), 'sparse', or 'dense'", default = 'bd') %>%
    add_argument("--grm_spd", help = "is matrix SPD?", flag = TRUE) %>%
    add_argument("--time", help = "time variable name") %>%
    add_argument("--status", help = "status variable name") %>%
    add_argument("--out_prefix", help = "Prefix for output files", default = "") %>%
    add_argument("--covars", help = "Covariate variable names (space-separated)", nargs = Inf) %>%
    add_argument("--sample_id", help = "File with vector of sample IDs") %>%
    add_argument("--variant_id", help = "File with vector of variant IDs") %>%
    add_argument("--maf", help = "minimum MAF of variants to test", default = 0.01) %>%
    add_argument("--score", help = "perform score test", flag = TRUE) %>%
    add_argument("--threshold", help = "reestimate HRs for variants with a p-value<threshold by first estimating a variant-specific variance component", default = 0) %>%
    add_argument("--chromosome", help = "chromosome number")

argv <- parse_args(argp)

library(SeqArray)
library(Biobase)
library(Matrix)
library(coxmeg)
library(dplyr)
sessionInfo()
print(argv)

intToChr <- function(chr) {
    if (is.na(chr)) return(NA)
    if (chr == 23) return("X")
    if (chr == 24) return("Y")
    if (chr == 25) return("XY")
    if (chr == 26) return("MT")
    as.character(chr)
}
chrom <- intToChr(argv$chromosome)

phen <- readRDS(argv$pheno_file)
phen_covars <- pData(phen) %>%
    mutate(family=sample.id)

if (!is.na(argv$sample_id)) {
    sample_id <- readRDS(argv$sample_id)
    phen_covars <- phen_covars %>%
    filter(sample.id %in% sample_id)
}

pheno <- phen_covars %>%
    select(family, sample.id, !!argv$time, !!argv$status)

if (is.na(argv$covars[1])) {
    cov <- NULL
} else {
    cov <- phen_covars %>%
        select(family, sample.id, one_of(argv$covars))
}

grm <- readRDS(argv$grm_file)
keep <- as.character(phen_covars$sample.id)
grm <- grm[keep, keep]

gds <- seqOpen(argv$gds_file)

if (!is.na(chrom)) seqSetFilterChrom(gds, chrom)

if (!is.na(argv$variant_id)) {
    variant_id <- readRDS(argv$variant_id)
    seqSetFilter(gds, variant.id = variant_id, action = "intersect")
}

res <- coxmeg_gds(gds=gds, pheno=pheno, corr=grm, type=argv$grm_type,
                  cov=cov, maf=argv$maf, score=argv$score,
                  threshold=argv$threshold, spd=argv$grm_spd)

if (!is.na(chrom)) {
    out_file <- paste0(argv$out_prefix, "coxmeg_chr", chrom, ".rds")
} else {
    out_file <- paste0(argv$out_prefix, "coxmeg.rds")
}
saveRDS(res, file=out_file)

seqClose(gds)
