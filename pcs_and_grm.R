#! /usr/bin/env Rscript
library(argparser)
library(SeqArray)
library(GENESIS)
library(SeqVarTools)

library(SNPRelate)
library(MASS)
library(ggplot2)

sessionInfo()

# read arguments
argp <- arg_parser("Generate PCs and GRM")
argp <- add_argument(argp, "gds_file", help = "GDS file")
argp <- add_argument(argp, "snpset", help = "File with vector of variant IDs")
argp <- add_argument(argp, "sampset", help = "File with vector of sample IDs")
argp <- add_argument(argp, "kin.thresh",
                     help = "Kinship threshold for pcair (2 ^ -kin.thresh)",
                     default = 5.5)
argv <- parse_args(argp)
if (!is.na(argv$snpset)) {
  snpset <- readRDS(argv$snpset)
} else {
  snpset <- NULL
}
if (!is.na(argv$snpset)) {
  sampset <- readRDS(argv$sampset)
} else {
  sampset <- NULL
}
kin.thresh <- 2 ^ (-argv$kin.thresh)
gds <- seqOpen(argv$gds_file)

king <- snpgdsIBDKING(gds, snp.id = snpset, sample.id = sampset)
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id
mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = kin.thresh,
                 divobj = kingMat, snp.include = snpset, sample.include = sampset)

seqSetFilter(gds, variant.id = snpset, sample.id = sampset)
seqData <- SeqVarData(gds)
print("1st iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, 1:3],
                    training.set = mypcair$unrels)
pcrelate_matrix <- pcrelateToMatrix(mypcrel, scaleKin=2, thresh = kin.thresh)

pca <- pcair(seqData, kinobj = pcrelate_matrix, kin.thresh = kin.thresh,
             divobj = kingMat, snp.include = snpset, scan.include = sampset)
