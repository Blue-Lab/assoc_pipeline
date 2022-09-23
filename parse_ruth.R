#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Parse RUTH output") %>%
  add_argument("ruth_file", help = "RUTH output file (.vcf)") %>%
  add_argument("--out_file", help = "Output filename (.rds)",
               default = "ruth.rds")
argv <- parse_args(argp)

sessionInfo()
print(argv)

library(tidyr)
library(dplyr)

ruth_res <- read.table(argv$ruth_file, sep = "\t", comment.char = "#",
                       na.strings = ".",
                       col.names = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                     "FILTER", "INFO")) %>%
  separate(INFO, sep = ";",
           into = c("PR", "FIBC_P", "HWE_SLP_P", "FIBC_I", "HWE_SLP_I",
                    "MAX_IF", "MIN_IF", "LLK0", "BETA_IF")) %>%
  mutate(across(FIBC_P:BETA_IF, ~ sub("^.*=", "", .x)))
# the BETA_IF field will be standard for all observations. Count the number of
# fields.
n_sep <- ruth_res$BETA_IF[1] %>% gregexpr(",", .) %>% extract2(1) %>% length()
out <- separate(ruth_res, BETA_IF, paste0("BETA_IF_", seq(n_sep + 1)), ",") %>%
  mutate(across(c(FIBC_P:LLK0, matches("BETA_IF")), as.numeric),
         ABS_PVAL = abs(HWE_SLP_I))

saveRDS(out, argv$out_file)
