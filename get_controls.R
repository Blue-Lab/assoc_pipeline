#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Get IDs of controls") %>%
  add_argument("pheno_file", help="Phenotype file in annotated dataframe format (.rds)") %>%
  add_argument("status", help = "Case/control status variable in pheno_file (binary)") %>%
  add_argument("--out_file", help = "Output filename", default = "controls.rds") %>%
  add_argument("--sample_id", help = "File with vector of sample IDs to keep (.rds)")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(Biobase)
library(dplyr)


pheno <- readRDS(argv$pheno_file) %>% pData()
controls <- filter(pheno, .data[[argv$status]] == 0)$sample.id

if (!is.na(argv$sample_id)) {
  sample_id <- readRDS(argv$sample_id)
  keep <- controls[controls %in% sample_id]
} else {
  keep <- controls
}

saveRDS(keep, argv$out_file)
