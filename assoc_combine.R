#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Association plots") %>%
    add_argument("assoc_prefix",
                     help = "Prefix for association test results file (_chr*.rds)") %>%
    add_argument("--delete_indiv", help = "delete individual chromosome files", flag = TRUE)
argv <- parse_args(argp)

library(dplyr)
sessionInfo()
print(argv)

## find all files and order by chromosome
file.pattern <- paste0(basename(argv$assoc_prefix), "assoc_chr[[:alnum:]]+.rds")
files <- list.files(path=dirname(argv$assoc_prefix), pattern=file.pattern, full.names=TRUE)
chrs <- sub("_chr", "", regmatches(basename(files), regexpr("_chr[[:alnum:]]+", basename(files))))
chrs <- factor(chrs, levels=c(0:26,"X","XY","Y"))
files <- files[order(chrs)]

assoc <- lapply(files, readRDS) %>%
    bind_rows()

outfile <- paste0(argv$assoc_prefix, "assoc.rds")
saveRDS(assoc, file=outfile)

## delete individual files
if (argv$delete_indiv) unlink(files)
