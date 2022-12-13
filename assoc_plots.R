#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Association plots") %>%
  add_argument("assoc_file", help = "Association test results file (.rds)") %>%
  add_argument("--maxP_m", help = "Maximum p-value for Mannhattan plot",
               type = "numeric", default = NULL) %>%
  add_argument("--maxP_q", help = "Maximum p-value for QQ plot",
               type = "numeric", default = 14) %>%
  add_argument("--img_ext", help="File extension for plots. Must be one of c('png', 'pdf', jpeg)", default = "png") %>%
  add_argument("--out_prefix", help = "Prefix for output files", default = "")
argv <- parse_args(argp)

library(dplyr)
library(fastman)
sessionInfo()
print(argv)

img_ext <- argv$img_ext

if (img_ext == "png") {
  plot_type <- png
} else if (img_ext == "pdf") {
  plot_type <- pdf
} else if (img_ext == "jpeg") {
  plot_type <- jpeg
} else {
  stop("--img_ext must be one of c('png', 'pdf', jpeg)")
}

assoc <- readRDS(argv$assoc_file)

## qq plot
paste0(argv$out_prefix, "qq.", img_ext) %>% plot_type(width=7, height=7, units="in", res=300)
fastqq(assoc$Score.pval, maxP = argv$maxP_q)
dev.off()

## manhattan plot
assoc <- mutate(assoc, chr = factor(chr, levels = c(0:26,"X","XY","Y")))

paste0(argv$out_prefix, "manh.", img_ext) %>% plot_type(width=10, height=6, units="in", res=300)
fastman(assoc, "chr", "pos", "Score.pval", maxP = argv$maxP_m)
dev.off()
