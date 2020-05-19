#! /usr/bin/env Rscript
library(argparser)

# read arguments
argp <- arg_parser("Convert BED to GDS")
argp <- add_argument(argp, "bed_prefix", help="prefix for bed/bim/fam")
argp <- add_argument(argp, "--out_file", help="output file name (default is <bed_prefix>.gds)")
argv <- parse_args(argp)
library(SeqArray)
sessionInfo()
print(argv)

bed.prefix <- argv$bed_prefix
out.file <- if (is.na(argv$out_file)) paste0(bed.prefix, ".gds") else argv$out_file

seqBED2GDS(bed.fn=paste0(bed.prefix, ".bed"),
           fam.fn=paste0(bed.prefix, ".fam"),
           bim.fn=paste0(bed.prefix, ".bim"),
           out.gdsfn=out.file)
