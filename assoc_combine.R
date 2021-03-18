#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("Association plots") %>%
    add_argument("assoc_prefix",
                     help = "Prefix for association test results file (_chr*.rds)") %>%
    add_argument("--delete_indiv", help = "delete individual chromosome files", flag = TRUE) %>%
    add_argument("--coxmeg", help = "association files were produced by coxmeg", flag = TRUE)
argv <- parse_args(argp)

library(dplyr)
sessionInfo()
print(argv)

suffix <- if (argv$coxmeg) "coxmeg" else "assoc"
    
## find all files and order by chromosome
file.pattern <- paste0(basename(argv$assoc_prefix), suffix, "_chr[[:alnum:]]+.rds")
files <- list.files(path=dirname(argv$assoc_prefix), pattern=file.pattern, full.names=TRUE)
chrs <- sub("_chr", "", regmatches(basename(files), regexpr("_chr[[:alnum:]]+", basename(files))))
chrs <- factor(chrs, levels=c(0:26,"X","XY","Y"))
files <- files[order(chrs)]
print(files)

if (argv$coxmeg) {
    assoc <- lapply(files, function(x) readRDS(x)$summary)
} else {
    assoc <- lapply(files, readRDS)
}
assoc <- bind_rows(assoc)

outfile <- paste0(argv$assoc_prefix, suffix, ".rds")
saveRDS(assoc, file=outfile)

## delete individual files
if (argv$delete_indiv) unlink(files)
