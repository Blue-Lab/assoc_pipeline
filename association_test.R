#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("Run association test") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pheno_file",
               help = "Phenotype file in annotated dataframe format") %>%
  add_argument("--out_prefix", help = "Prefix for output files", default = "") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs") %>%
  add_argument("--chromosome", help = "chromosome number") %>%
  add_argument("--null_model", help = "null model object (.rds)") %>%
  add_argument("--genome_build", help = 'Passed to assocTestSingle. Can be "hg19" or "hg38"') %>%
  add_argument("--dosage", help = "read dosage from DS node", flag = TRUE)

argv <- parse_args(argp)

library(SeqArray)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
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

gds <- seqOpen(argv$gds_file)
pheno <- readRDS(argv$pheno_file)

if (!is.na(argv$variant_id)) {
  variant_id <- readRDS(argv$variant_id)
} else {
  variant_id <- NULL
}
if (!is.na(argv$sample_id)) {
  sample_id <- readRDS(argv$sample_id)
} else {
  sample_id <- NULL
}
if (is.na(argv$genome_build)) {
  genome_build <- NULL
} else {
  genome_build <- argv$genome_build
}

gds.id <- seqGetData(gds, "sample.id")
seqData <- SeqVarData(gds, sampleData = pheno)
if (!is.na(chrom)) seqSetFilterChrom(gds, chrom)
seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id, action = "intersect")
iterator <- SeqVarBlockIterator(seqData, verbose=TRUE)

nullmod <- readRDS(argv$null_model)

assoc <- assocTestSingle(iterator, nullmod, imputed=argv$dosage,
                         genome.build = genome_build)

if (!is.na(chrom)) {
    out_file <- paste0(argv$out_prefix, "assoc_chr", chrom, ".rds")
} else {
    out_file <- paste0(argv$out_prefix, "assoc.rds")
}
saveRDS(assoc, out_file)
seqClose(gds)
