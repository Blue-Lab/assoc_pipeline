#! /usr/bin/env Rscript
library(argparser)
library(magrittr)

argp <- arg_parser("PCA correlation plots") %>%
  add_argument("corr_file", help = "PCA correlation file") %>%
  add_argument("--n_pcs", default = 32,
               "Number of PC correlation plots to draw") %>%
  add_argument("--n_perpage", help = "Number of plots per page",default = 4) %>%
  add_argument("--out_prefix", help = "prefix for output file names", default =
               "pca_corr") %>%
  add_argument("--img_ext", help="File extension for plots", default = "png") %>%
  add_argument("--dense",
               help = "Skip sparsifying plots? This will increase computation time.",
               flag = TRUE) 

argv <- parse_args(argp)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gdsfmt)
sessionInfo()
print(argv)

# Function definitions. These are some internal TOPMed Analysis Pipeline
# functions replicated here.
thinPoints <- function(dat, value, n=3000, nbins=200, groupBy=NULL) {
  if (!(value %in% colnames(dat)) || !is.numeric(dat[[value]])) {
    stop("value should be name of a numeric column.")
  }
  if (!is.numeric(n) || !is.numeric(nbins)) stop("n and nbins should numbers.")
  if (!is.character(value) || (length(value) != 1)) stop("value should be a string.")

  if (!is.null(groupBy))
  {
    if (!is.character(groupBy) || (length(groupBy) != 1)) {
      stop("groupBy should be a character.")
    } else if (!(groupBy %in% colnames(dat))) {
      stop(sprintf("Column %s does not exist in dat.", groupBy))
    }
    group <- integer(nrow(dat))
    group[order(dat[[groupBy]])] <- unlist(tapply(dat[[value]], dat[[groupBy]], function(x) cut(x, breaks = nbins, labels = FALSE)))
    group <- interaction(dat[[groupBy]], group)
  } else {
    group <- factor(cut(dat[[value]], breaks = nbins, labels = FALSE))
  }

  indices <- unlist(tapply(1:nrow(dat), group, FUN = function(x) sample_vec(x, n)), use.names = FALSE)
  dat[indices,]
}

sample_vec <- function(x, n) {
  if (length(x) == 1) {
    return(x)
  } else if (length(x) <= n) {
    return(x)
  } else {
    return(sample(x, size = n, replace = FALSE))
  }
}

n_pcs <- argv$n_pcs

cf <- openfn.gds(argv$corr_file)
corr <- t(read.gdsn(index.gdsn(cf, "correlation")))
corr <- corr[,1:n_pcs]
ms <- rowSums(is.na(corr)) == n_pcs # monomorphic variants
corr <- corr[!ms,]
colnames(corr) <- paste0("PC", 1:n_pcs)
corr <- data.frame(corr,
                   chr=readex.gdsn(index.gdsn(cf, "chromosome"), sel=!ms),
                   pos=readex.gdsn(index.gdsn(cf, "position"), sel=!ms),
                   stringsAsFactors=FALSE)
closefn.gds(cf)

corr %<>% gather(PC, value, -chr, -pos) %>%
    filter(!is.na(value)) %>%
    mutate(value=abs(value)) %>%
    mutate(PC=factor(PC, levels=paste0("PC", 1:n_pcs)))

n_pcs <- min(as.integer(n_pcs), ncol(corr))
if (!argv$dense) corr %<>% thinPoints("value", n=10000, nbins=10, groupBy="PC")

corr %<>% mutate(chr=factor(chr, levels=c(1:22, "X")))
chr <- levels(corr$chr)
cmap <- setNames(rep_len(brewer.pal(8, "Dark2"), length(chr)), chr)

n_pcs <- length(unique(corr$PC))
n_plots <- ceiling(n_pcs/argv$n_perpage)
bins <- as.integer(cut(1:n_pcs, n_plots))
for (i in 1:n_plots) {
    bin <- paste0("PC", which(bins == i))
    dat <- filter(corr, PC %in% bin)

    p <- ggplot(dat, aes(chr, value, group=interaction(chr, pos), color=chr)) +
        geom_point(position=position_dodge(0.8)) +
        facet_wrap(~PC, scales="free", ncol=1) +
        scale_color_manual(values=cmap, breaks=names(cmap)) +
        ylim(0,1) +
        theme_bw() +
        theme(legend.position="none") +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        xlab("chromosome") + ylab("abs(correlation)")
    ggsave(paste0(argv$out_prefix , i, ".", argv$img_ext), plot=p, width=10, height=15)
}
