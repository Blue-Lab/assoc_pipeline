#! /usr/bin/env Rscript
library(argparser)
library(SeqArray)
library(GENESIS)
library(SeqVarTools)
library(magrittr)

sessionInfo()

# read arguments
argp <- arg_parser("Generate PCs and GRM")
argp %<>% add_argument("gds_file", help = "GDS file")
argp %<>% add_argument("--out_prefix", help = "Prefix for output files")
argp %<>% add_argument("--snpset", help = "File with vector of variant IDs")
argp %<>% add_argument("--sampset", help = "File with vector of sample IDs")
argp %<>% add_argument("--kin.thresh", default = 5.5,
                       help = "Kinship threshold for pcair (2 ^ -kin.thresh)")
argp %<>% add_argument("--keep_king", flag = TRUE,
                       help = "Save KING-robust GRM")
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
out_prefix <- ifelse(!is.na(argv$out_prefix), argv$out_prefix, "")
gds <- seqOpen(argv$gds_file)

king <- snpgdsIBDKING(gds, snp.id = snpset, sample.id = sampset)
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

if (argv$keep_king) {
  kingMat_temp <- kingMat * 2 # Scaled to match pc-relate GRM
  # coerces low values in matrix to 0
  kingMat_temp[kingMat_temp <= kin.thresh] <- 0
  saveRDS(kingMat_temp, paste0(out_prefix, "king_robust_grm.rds"))
  rm(kingMat_temp)
}


mypcair <- pcair(gds, kinobj = kingMat, kin.thresh = kin.thresh,
                 divobj = kingMat, snp.include = snpset,
                 sample.include = sampset)

seqSetFilter(gds, variant.id = snpset, sample.id = sampset)
seqData <- SeqVarData(gds)
print("1st iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
mypcrel <- pcrelate(iterator, pcs = mypcair$vectors[, 1:3],
                    training.set = mypcair$unrels)
pcrelate_matrix <- pcrelateToMatrix(mypcrel, scaleKin=2, thresh = kin.thresh)

pca <- pcair(seqData, kinobj = pcrelate_matrix, kin.thresh = kin.thresh,
             divobj = kingMat, snp.include = snpset, scan.include = sampset)

resetIterator(iterator, verbose = TRUE)

print("2nd iteration PC-relate")
iterator <- SeqVarBlockIterator(seqData, verbose = FALSE)
pcrel2 <- pcrelate(iterator, pcs = pca$vectors[, 1:3],
                   training.set = pca$unrels)


pcrelate_matrix <- pcrelateToMatrix(pcrel2, scaleKin = 2, thresh = kin.thresh)
saveRDS(pca$vectors, paste0(out_prefix, "pcs.rds"))
saveRDS(pcrelate_matrix, paste0(out_prefix, "pcr_grm.rds"))
