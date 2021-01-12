#! /usr/bin/env Rscript
library(argparser)
library(magrittr)
argp <- arg_parser("Generate PC plots") %>%
  add_argument("pcair_file", help = "PC-AiR file (.rds)") %>%
  add_argument("--out_prefix", help = "Prefix for output files",
               default = "") %>%
  add_argument("--phenotype_file", help = "Phenotype file (.rds)") %>%
  add_argument("--group", help = "grouping variable") %>%
  add_argument("--n_pairs", help = "number of pairwise plots", default = 6)

argv <- parse_args(argp)

library(Biobase)
library(dplyr)
library(ggplot2)
library(GGally)
sessionInfo()
print(argv)

out_prefix <- argv$out_prefix

pca <- readRDS(argv$pcair_file)
pcs <- as.data.frame(pca$vectors[pca$unrels,])
n <- ncol(pcs)
names(pcs) <- paste0("PC", 1:n)
pcs$sample.id <- as.integer(row.names(pcs))

## scree plot
dat <- data.frame(pc = seq(n), varprop=pca$varprop[seq(n)])
p <- ggplot(dat, aes(x=factor(pc), y=100*varprop)) +
  geom_point() + theme_bw() +
  xlab("PC") + ylab("Percent of variance accounted for")
ggsave(paste0(out_prefix, "pc_scree.png"), plot=p, width=6, height=6)

## color by group
if (!is.na(argv$phenotype_file) & !is.na(argv$group)) {
    group <- argv$group
    annot <- readRDS(argv$phenotype_file)
	 #pData(annot)$sample.id <- as.character(pData(annot)$sample.id)
    if ("AnnotatedDataFrame" %in% class(annot)) annot %<>% pData()
    stopifnot(group %in% names(annot))
    annot %<>% select(sample.id, !!enquo(group))
    pcs <- left_join(pcs, annot, by="sample.id")
} else {
    ## make up dummy group
    group <- "group"
    pcs$group <- "NA"
}

p12 <- ggplot(pcs, aes_string("PC1", "PC2", color=as.character(group))) + geom_point() +
    guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(paste0(out_prefix, "pc12.png"), plot=p12, width=7, height=6)


p34 <- ggplot(pcs, aes_string("PC3", "PC4", color=as.character(group))) + geom_point() +
    guides(colour=guide_legend(override.aes=list(alpha=1)))
ggsave(paste0(out_prefix, "pc34.png"), plot=p34, width=7, height=6)



npr <- min(argv$n_pairs, n)
p <- ggpairs(pcs, mapping=aes_string(color=as.character(group)), columns=1:npr,
             lower=list(continuous=wrap("points")),
             diag=list(continuous="densityDiag"),
             upper=list(continuous="blank"))
png(paste0(out_prefix, "pairs.png"), width=8, height=8, units="in", res=150)
print(p)
dev.off()


pc2 <- pcs
names(pc2)[1:ncol(pc2)] <- sub("PC", "", names(pc2)[1:ncol(pc2)])

p <- ggparcoord(pc2, columns=1:n, groupColumn=as.character(group), scale="uniminmax") +
    guides(colour=guide_legend(override.aes=list(size=2))) +
    xlab("PC") + ylab("")
ggsave(paste0(out_prefix, "parcoord.png"), plot=p, width=10, height=5)
