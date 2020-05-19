#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("Run association test") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pheno_file",
               help = "Phenotype file in annotated dataframe format") %>%
  add_argument("grm_file", help = "GRM file") %>%
  add_argument("outcome", help = "Outcome variable name") %>%
  add_argument("family", help = "Distribution family",
               default = "gaussian") %>%
  add_argument("--out_file", help="output file name",
               default = "assoc.rds") %>%
  add_argument("--covars",
               help = "Covariate variable names (space-separated)") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs")

argv <- parse_args(argp)

library(SeqArray)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()
print(argv)

gds <- seqOpen(argv$gds_file)
pheno <- readRDS(argv$pheno_file)

if (!is.na(argv$variant_id)) {
  variant_id <- readRDS(argv$variant_id)
} else {
  variant_id <- NULL
}
if (!is.na(argv$variant_id)) {
  sample_id <- readRDS(argv$sample_id)
} else {
  sample_id <- NULL
}

if (!is.na(argv$covars)) {
  covars <- strsplit(argv$covars, " ") %>% unlist
} else {
  covars <- NULL
}

gds.id <- seqGetData(gds, "sample.id")
seqData <- SeqVarData(gds, sampleData = pheno)
seqSetFilter(gds, variant.id = variant_id, sample.id = sample_id)
iterator <- SeqVarBlockIterator(seqData, verbose=TRUE)

grm <- readRDS(argv$grm_file)
nullmod <- fitNullModel(pheno, outcome = argv$outcome, covars = covars,
                        cov.mat = grm, family = argv$family, verbose=FALSE,
                        sample.id = sample_id)

message("Null model fixed effects:")
message(nullmod$fixef)


assoc <- assocTestSingle(iterator, nullmod)
saveRDS(assoc, argv$out_file)
seqClose(gds)
