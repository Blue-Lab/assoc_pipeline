#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("Correlation of variants with PCs") %>%
  add_argument("gds_file", help = "GDS file") %>%
  add_argument("pcair_file", help = "File with pcair object") %>%
  add_argument("--n_pcs", default = 32,
               "Number of PCs to calculate correlation") %>%
  add_argument("--out_file", default = "pca_corr.gds",
               help = "output filename") %>%
  add_argument("--variant_id", help = "File with vector of variant IDs") %>%
  add_argument("--num_core", help = "number of cores", default = 6)
argv <- parse_args(argp)

library(SeqVarTools)
library(SNPRelate)
library(gdsfmt)
sessionInfo()
print(argv)

n_pcs <- argv$n_pcs
outfile <- argv$out_file

gds <- seqOpen(argv$gds_file)
pcaobj <- readRDS(argv$pcair_file)

# Workaround so that the structure of our PC file matches what snpgdsPCACorr()
# expects.

if (!(inherits(pcaobj, "snpgdsPCAClass") | inherits(pcaobj, "snpgdsEigMixClass"))) {
  warning('pcaobj does inherit expected class. Attempting to coerce to "snpgdsEigMixClass"')
  class(pcaobj) <- "snpgdsEigMixClass"
}

if (!("eigenvect" %in% names(pcaobj))) {
  warning("`pcaobj$eigenvect` is not present. Attempting to use `pcaobj$vectors`")
  pcaobj$eigenvect <- pcaobj$vectors
}

variant.id <- if (!is.na(argv$variant_id)) readRDS(argv$variant_id) else NULL
sample_include <- c(pcaobj$rels, pcaobj$unrels)
seqSetFilter(gds, variant.id, sample.id = sample_include)

variant.id <- seqGetData(gds, "variant.id")
snpgdsPCACorr(pcaobj, gdsobj = gds, snp.id = variant.id, eig.which = 1:n_pcs,
              num.thread = argv$num_core, outgds = outfile)

# TOPMed Pipeline implementation included a second `seqSetFilter()` call.
seqSetFilter(gds, variant.id = variant.id)
chromosome <- seqGetData(gds, "chromosome")
position <- seqGetData(gds, "position")

seqClose(gds)

pca.corr <- openfn.gds(outfile, readonly=FALSE)
add.gdsn(pca.corr, "chromosome", chromosome, compress="LZMA_RA")
add.gdsn(pca.corr, "position", position, compress="LZMA_RA")
closefn.gds(pca.corr)
cleanup.gds(outfile)

# Memory stats
ms <- gc()
cat(">>> Max memory: ", ms[1,6]+ms[2,6], " MB\n")
