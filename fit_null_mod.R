#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("Fit null model") %>%
  add_argument("pheno_file",
               help = "Phenotype file in annotated dataframe format") %>%
  add_argument("--grm_file", help = "GRM file") %>%
  add_argument("--outcome", help = "Outcome variable name") %>%
  add_argument("--family", help = "Distribution family",
               default = "gaussian") %>%
  add_argument("--out_file", help="output file name",
               default = "nullmod.rds") %>%
  add_argument("--covars",
               help = "Covariate variable names (space-separated)", nargs=Inf) %>%
  add_argument("--sample_id", help = "File with vector of sample IDs")

argv <- parse_args(argp)

library(SeqArray)
library(SeqVarTools)
library(Biobase)
library(GENESIS)
sessionInfo()
print(argv)

pheno <- readRDS(argv$pheno_file)

if (!is.na(argv$sample_id)) {
  sample_id <- readRDS(argv$sample_id)
} else {
  sample_id <- NULL
}

if (is.na(argv$covars[1])) {
  argv$covars <- NULL
}

if (!is.na(argv$grm_file)) {
    grm <- readRDS(argv$grm_file)
} else {
    grm <- NULL
}

nullmod <- fitNullModel(pheno, outcome = argv$outcome, covars = argv$covars,
                        cov.mat = grm, family = argv$family, verbose=FALSE,
                        sample.id = sample_id)

message("Null model fixed effects:")
message(nullmod$fixef)

saveRDS(nullmod, argv$out_file)
