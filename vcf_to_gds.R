#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

# read arguments
argp <- arg_parser("Convert VCF to GDS") %>%
    add_argument("vcf_file", help="vcf file name") %>%
    add_argument("out_file", help="output file name") %>%
    add_argument("--dosage", help="import dosage (DS)", flag=TRUE) %>%
    add_argument("--num_cores", help="Number of cores to utilize for parallel processing", default=1)
argv <- parse_args(argp)

library(SeqArray)
sessionInfo()
print(argv)

if (argv$dosage) {
    seqVCF2GDS(argv$vcf_file, argv$out_file,
               fmt.import="DS", scenario="imputation",
               parallel=argv$num_cores)
} else {
    seqVCF2GDS(argv$vcf_file, argv$out_file,
               fmt.import=character(0),
               parallel=argv$num_cores)
}
