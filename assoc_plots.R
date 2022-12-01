#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Association plots") %>%
  add_argument("assoc_file", help = "Association test results file (.rds)") %>%
  add_argument("--maxP_m", help = "Maximum p-value for Mannhattan plot",
               type = "numeric", default = NULL) %>%
  add_argument("--maxP_q", help = "Maximum p-value for QQ plot",
               type = "numeric", default = 14) %>%
  add_argument("--out_prefix", help = "Prefix for output files", default = "")
argv <- parse_args(argp)

library(dplyr)
library(fastman)
sessionInfo()
print(argv)

assoc <- readRDS(argv$assoc_file)

## qq plot
paste0(argv$out_prefix, "qq.png") %>% png(width=7, height=7, units="in", res=300)
fastqq(assoc$Score.pval, maxP = argv$maxP_q)
dev.off()

## manhattan plot
assoc <- mutate(assoc, chr = factor(chr, levels = c(0:26,"X","XY","Y")))

paste0(argv$out_prefix, "manh.png") %>% png(width=10, height=6, units="in", res=300)
fastman(assoc, "chr", "pos", "Score.pval", maxP = argv$maxP_m)
dev.off()
