# Benjamin J.M. Tremblay
#
# Figure 5 and Supplementary Figures 5, 6
#
# Package requirements: rtracklayer, GenomicFeatures, matrixStats,
# EnrichedHeatmap, ComplexHeatmap, BSgenome.Athaliana.TAIR.TAIR9,
# universalmotif, ggplot2, dplyr, readr, zoo, viridisLite, genomation, rstatix,
# cowplot

library(rtracklayer)
library(GenomicFeatures)
library(matrixStats)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(ggplot2)
library(dplyr)

source("Functions.R")

dir.create("Figure5", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

ARAPORT11_GFF <- "Araport11_novel.gff3"
TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_BED <- "tss.size1.bed"
ANTISENSES <- "Antisenses.txt"
GENES_BED <- "genes.bed"
CLUSTERS_TABLE <- "csRNAseq_clustering/Clusters_anno.tsv"

# Read density
ATAC_57H_BW <- "atacseq_bw/L57.bw"
RNAPII_BW <- "external_bigWigs/RNAPII_14D.bw"
SER5P_BW <- "external_bigWigs/RNAPIISer5P_14D.bw"
SER2P_BW <- "external_bigWigs/RNAPIISer2P_14D.bw"
H3K4ME3_BW <- "external_bigWigs/H3K4me3_10D.bw"
H3K9AC_BW <- "external_bigWigs/H3K9ac_12D.bw"
PHYLOP_BW <- "PhyloP.bw"
PHASTCONS_BW <- "PhastCons.bw"
MNASE_BW <- "Ath_leaf_NPS.bw"

# Shared data ------------------------------------------------------------------

seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

a11 <- makeTxDbFromGFF(ARAPORT11_GFF)

tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tss1 <- import(TSS_BED)

tssM <- as.matrix(tssA[, 10:21])
rownames(tssM) <- tssA$TSS
tssMm <- as.matrix(data.frame(row.names = rownames(tssM),
  DS = rowMeans(tssM[, 1:2]),
  S24 = rowMeans(tssM[, 3:4]),
  S72 = rowMeans(tssM[, 5:6]),
  L6 = rowMeans(tssM[, 7:8]),
  L26 = rowMeans(tssM[, 9:10]),
  L57 = rowMeans(tssM[, 11:12])
))

pc_bed <- import(GENES_BED)

tssAcls <- readr::read_tsv(CLUSTERS_TABLE)

theme_5pt <- theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    legend.key.size = unit(5, "pt"),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

# Antisense pairs --------------------------------------------------------------

pair_dist <- function(x, y) {
  distance(tss1[tss1$name==x], tss1[tss1$name==y], ignore.strand=TRUE)
}

a11_g <- genes(a11)
g2len <- structure(width(a11_g), names = a11_g$gene_id)
final_as <- readr::read_tsv(ANTISENSES)
final_as$SenseLen <- g2len[final_as$PcGene]
final_as$AntisenseLen <- g2len[final_as$AntisenseGene]
final_as$TSS_dist <- NA_integer_
for (i in 1:nrow(final_as)) {
  final_as$TSS_dist[i] <- pair_dist(final_as$PcTSS[i], final_as$AntisenseTSS[i])
}
final_as$TTS_dist <- final_as$SenseLen - final_as$TSS_dist

final_as$Dist_ratio <- with(final_as, pmin(1,(SenseLen - TTS_dist) / SenseLen))

# Proximal: within the first half of the gene and more than 1 kbp from the TTS
final_as$Proximal <- final_as$Dist_ratio<=0.5 & final_as$TTS_dist>1000

readr::write_lines(final_as$PcGene[final_as$Proximal], "Figure5/antisenses_proximal.txt")
readr::write_lines(final_as$PcGene[!final_as$Proximal], "Figure5/antisenses_distal.txt")

# Sense, antisense and no-antisense TSS sets
tss1_as <- tss1[tss1$name %in% final_as$AntisenseTSS]
tss1_as$Proximal <- tss1_as$name %in% final_as$AntisenseTSS[final_as$Proximal]
tss1_pc <- tss1[tss1$name %in% final_as$PcTSS]
tss1_pc$Proximal <- tss1_pc$name %in% final_as$PcTSS[final_as$Proximal]

non_as_pc <- tssA$TSS[!is.na(tssA$NearestGene) & tssA$NearestGeneType == "protein_coding"]
non_as_pc <- non_as_pc[!non_as_pc %in% tss1_pc$name]
tss1_pc_non_as <- tss1[tss1$name %in% non_as_pc]

# Expression-matched random sample of TSSs without an antisense
dd <- density(bw=1,log2(1+tssMm[tss1_pc$name[!tss1_pc$Proximal],"L57"]))
tss1_pc_non_as2 <- tss1_pc_non_as[sample(1:length(tss1_pc_non_as),length(tss1_pc),prob=approxfun(dd$x,dd$y,yleft=0,yright=10)(log2(1+tssMm[tss1_pc_non_as$name,"L57"])))]

tss1_a <- tss1[tss1$name %in% non_as_pc]
tss1_b <- tss1[tss1$name %in% tss1_pc$name[tss1_pc$Proximal]]
tss1_c <- tss1[tss1$name %in% tss1_pc$name[!tss1_pc$Proximal]]
tss1_d <- tss1[tss1$name %in% tss1_as$name[tss1_as$Proximal]]
tss1_e <- tss1[tss1$name %in% tss1_as$name[!tss1_as$Proximal]]

zm <- function(x) zoo::rollmean(x, 10, na.pad=T)

# Figure 5a --------------------------------------------------------------------

length(unique(tssA$NearestGene[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]))
# 19,260 total genes detected, 3985 with an antisense

ggplot(data.frame(type=c("no", "yes"),N=c(19260-3985,3985)), aes(x="",y=N,fill=type)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c("grey", "black")) +
  coord_polar("y", start = 0) +
  theme_void()
ggsave("Figure5/Fig5a_pie3_gene_antisense.pdf", width = 3, height = .7)

# Figure 5b --------------------------------------------------------------------

ggplot(final_as, aes(pmin(7500,SenseLen), pmin(7500,TTS_dist), colour = Proximal)) +
  geom_point(size = .2, pch = 19) +
  scale_colour_manual(values = c("black", "red"), name = element_blank(),
    labels = c("Distal", "Proximal")) +
  scale_y_continuous(name = "Distance between antisense\nTSS and gene TTS (bp)",
    limits = c(-800, 8000),
    expand = c(0, 0),
    breaks = seq(0, 7500, 2500), labels = c("0", "2500", "5000", ">7500")) +
  scale_x_continuous(name = "Gene length (bp)",
    limits = c(-800, 8000),
    expand = c(0, 0),
    breaks = seq(0, 7500, 2500), labels = c("0", "2500", "5000", ">7500")) +
  theme_5pt +
  theme(aspect.ratio=1)
ggsave("Figure5/Fig5b_antisense_dist.pdf", width = 2.52, height = 1.33)

# Figure 5d --------------------------------------------------------------------

tss1_pc$AntiCorr <- structure(final_as$Correlation, names = final_as$PcTSS)[tss1_pc$name]
tss1_pc$MaxCPM <- structure(final_as$MaxPcCPM, names = final_as$PcTSS)[tss1_pc$name]

tss1_pcDf <- as.data.frame(mcols(tss1_pc))

ggplot() +
  geom_smooth(data=filter(tss1_pcDf,Proximal),aes(AntiCorr,log2(1+MaxCPM)),
    colour="red") +
  geom_smooth(data=filter(tss1_pcDf,!Proximal),aes(AntiCorr,log2(1+MaxCPM)),
    colour="black") +
  scale_x_continuous(limits = c(-1, 1), name = "Correlation between sense/antisense") +
  ylab("Log2 Max CPM") +
  theme_5pt
ggsave("Figure5/Fig5d_corr_sense_antisense_maxcpm.pdf", width = 1.5, height = 1.35)

# Figure 5e --------------------------------------------------------------------

mnase3 <- import(MNASE_BW)
seqlevels(mnase3) <- as.character(1:5)

tss2len <- structure(tssA$NearestGene, names = tssA$TSS)
tss2len <- tss2len[!is.na(tss2len)]
tss2len <- tss2len[unname(tss2len) %in% names(g2len)]
tss2len <- structure(unname(g2len[unname(tss2len)]), names = names(tss2len))

tss2dist <- structure(final_as$TSS_dist, names = final_as$PcTSS)

tss1_a2 <- tss1_a[sample(1:length(tss1_a), 4000)]
no_as_win <- promoters(tss1_a2[order(tss2len[tss1_a2$name])], 500, 4000)
prox_as_win <- promoters(tss1_b[order(tss2dist[tss1_b$name])], 500, 4000)
dist_as_win <- promoters(tss1_c[order(tss2dist[tss1_c$name])], 500, 4000)

# MNase-seq heatmaps: no antisense (sorted by gene length), proximal and distal
# antisense (sorted by distance)
# Note: the original code for these heatmaps was lost, this is a recreation to
# reproduce the published heatmaps near exactly

no_as_map <- get_unstranded_signal(mnase3, no_as_win, keep = c(.00, .95), return.matrix = TRUE)
prox_as_map <- get_unstranded_signal(mnase3, prox_as_win, keep = c(.00, .95), return.matrix = TRUE)
dist_as_map <- get_unstranded_signal(mnase3, dist_as_win, keep = c(.00, .95), return.matrix = TRUE)

mnase_maps <- list(
  Fig5e_heatmap_mnase3_no_antisense_4k = list(m = no_as_map, h = 2.5),
  Fig5e_heatmap_mnase3_proxi_antisense = list(m = prox_as_map, h = 0.71),
  Fig5e_heatmap_mnase3_distal_antisense = list(m = dist_as_map, h = 1.72)
)
for (i in names(mnase_maps)) {
  pdf(paste0("Figure5/", i, ".pdf"), width = 2, height = mnase_maps[[i]]$h)
  draw(Heatmap(mnase_maps[[i]]$m,
    col=viridisLite::plasma(100,begin=0),
    show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
    cluster_rows = FALSE, cluster_columns= FALSE
  ))
  dev.off()
}

# Gene lengths
pdf("Figure5/Fig5e_heatmap_mnase3_proxi_antisense_lengths.pdf", width = 4, height = 3)
plot(1:length(tss1_b), tss2len[tss1_b$name][order(tss2dist[tss1_b$name])],
  pch = 19, cex = .4, type = "h", ylim = c(0, 20000)
)
dev.off()

pdf("Figure5/Fig5e_heatmap_mnase3_distal_antisense_lengths.pdf", width = 4, height = 3)
plot(1:length(tss1_c), tss2len[tss1_c$name][order(tss2dist[tss1_c$name])],
  pch = 19, cex = .4, type = "h", ylim = c(0, 20000)
)
dev.off()

# Figure 5f, g -----------------------------------------------------------------

tss1_pc$cls <- structure(tssAcls$Cluster, names = tssAcls$TSS)[tss1_pc$name]
tss1_pc$antiTSS <- structure(final_as$AntisenseTSS, names = final_as$PcTSS)[tss1_pc$name]
tss1_pc$CorrType <- "."
tss1_pc$CorrType[tss1_pc$AntiCorr>=.5] <- "+ve"
tss1_pc$CorrType[tss1_pc$AntiCorr<=-.25] <- "-ve"

tssMmZ <- t(scale(t(log2(1+tssMm))))

f5 <- tss1_pc$CorrType=="+ve" & tss1_pc$cls!=0
f6 <- tss1_pc$CorrType=="-ve" & tss1_pc$cls!=0

z_sets <- list(
  Fig5f_z_poscorr = f5,
  Fig5g_z_negcorr = f6
)
for (i in names(z_sets)) {
  fi <- z_sets[[i]]
  pdf(paste0("Figure5/", i, ".pdf"), width = 2., height = 1.7)
  draw(Heatmap(tssMmZ[tss1_pc$name[fi], ],
    row_order = order(tss1_pc$cls[fi]), name = "Z-score",
    row_title_rot = 0, row_title_gp = gpar(fontsize = unit(6, "pt")),
    column_names_rot = 0, column_names_gp = gpar(fontsize = unit(6, "pt")),
    row_gap = unit(0.2, "mm"),
    row_split = paste0("C", tss1_pc$cls[fi]),
    use_raster=T,raster_quality=3,
    cluster_columns = FALSE, show_row_names = FALSE
  ) +
  Heatmap(tssMmZ[tss1_pc$antiTSS[fi], ],
    show_heatmap_legend = FALSE,
    column_names_rot = 0, column_names_gp = gpar(fontsize = unit(6, "pt")),
    use_raster=T,raster_quality=3,
    cluster_columns = FALSE, show_row_names = FALSE
  ))
  dev.off()
}

getg <- function(x) {
  y <- tssA$NearestGene[tssA$TSS %in% x]
  y <- y[!is.na(y)]
  y <- y[grepl("^AT", y)]
  y
}

f1 <- !tss1_pc$Proximal & tss1_pc$CorrType=="+ve" & tss1_pc$cls != 0
f2 <- tss1_pc$Proximal & tss1_pc$CorrType=="+ve" & tss1_pc$cls != 0
f3 <- !tss1_pc$Proximal & tss1_pc$CorrType=="-ve" & tss1_pc$cls != 0
f4 <- tss1_pc$Proximal & tss1_pc$CorrType=="-ve" & tss1_pc$cls != 0

readr::write_lines(getg(tss1_pc$name[f1]), "Figure5/genes_distal_poscorr.txt")
readr::write_lines(getg(tss1_pc$name[f2]), "Figure5/genes_proxi_poscorr.txt")
readr::write_lines(getg(tss1_pc$name[f3]), "Figure5/genes_distal_negcorr.txt")
readr::write_lines(getg(tss1_pc$name[f4]), "Figure5/genes_proxi_negcorr.txt")

# csRNA-seq quantification of individual TSSs ----------------------------------

xa <- c(1.25, 1.75, 3.25, 3.75, 5.25, 5.75, 7.25, 7.75, 9.25, 9.75, 11.25, 11.75)
xb <- c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5)

plot_tss_expr <- function(tss, ylim, file) {
  pdf(file, width = 3.9, height = 3.5)
  plot(xa,
    tssM[tss, ], ylim = c(0, ylim),
    pch = 19, bty = "n", xaxt = "n", xlab = NA, ylab = "csRNA-seq (CPM)"
  )
  lines(xb, tssMm[tss, ])
  axis(1, xb, c("DS", "S24", "S72", "L6", "L26", "L57"))
  dev.off()
}

# Supplementary Figure 5a ------------------------------------------------------

# FLC and COOLAIR
plot_tss_expr("TSS5_3179424", 120, "Figure5/Supp5a_flc_expr_sense.pdf")
plot_tss_expr("TSS5_3173457", 60, "Figure5/Supp5a_flc_expr_antisense.pdf")

# Supplementary Figure 5c ------------------------------------------------------

tss1_pc$asTSS <- structure(final_as$AntisenseTSS,names=final_as$PcTSS)[tss1_pc$name]
tss1_pc$dist <- abs(start(tss1_pc) - as.integer(gsub("TSS\\d_","",tss1_pc$asTSS)))

dens <- density(pmin(2000,tss1_pc$dist[tss1_pc$Proximal]))
dens$x[which.max(dens$y)]

pdf("Figure5/Supp5c_proximal_density.pdf", width = 3.8,height=3)
plot(density(pmin(2000,tss1_pc$dist[tss1_pc$Proximal]),bw=50),
  xlab = "Distance from PcTSS (bp)", main = NA
)
abline(v=492,lty=2,col="grey",lwd=2)
dev.off()

# Supplementary Figure 5d ------------------------------------------------------

# Quantiles used for the dashed lines
2^quantile(log2((1+final_as$MaxPcCPM)/(1+final_as$MaxAntisenseCPM)),p=c(.025,.1,.25,.5,.75,.9,.975))

ggplot(final_as, aes(log2((1+final_as$MaxPcCPM)/(1+final_as$MaxAntisenseCPM)))) +
  geom_vline(xintercept = c(1.43,3.04,4.667), colour = "grey50", lty = 2) +
  geom_density(bw=.5) +
  coord_cartesian(xlim = c(-6, 12)) +
  ylab("Density") +
  scale_x_continuous(name = "Sense/antisense ratio",
    breaks = c(-3, 0, 3, 6, 9),
    labels = c("1 / 8",  "1", "8", "64", "512")
  ) +
  theme_5pt
ggsave("Figure5/Supp5d_antisense_ratio.pdf", height = 1.33, width = 1.5)

# Supplementary Figure 5e ------------------------------------------------------

phylop <- import(PHYLOP_BW)

no_as_phylop <- normalizeToMatrix(phylop, tss1_a, w = 10,
  extend = c(1000, 1000), value_column = "score")
prox_as_phylop <- normalizeToMatrix(phylop, tss1_b, w = 10,
  extend = c(1000, 1000), value_column = "score")
dist_as_phylop <- normalizeToMatrix(phylop, tss1_c, w = 10,
  extend = c(1000, 1000), value_column = "score")

prox_as_phylop2 <- normalizeToMatrix(phylop, tss1_d, w = 10,
  extend = c(1000, 1000), value_column = "score")
dist_as_phylop2 <- normalizeToMatrix(phylop, tss1_e, w = 10,
  extend = c(1000, 1000), value_column = "score")

pdf("Figure5/Supp5e_antisense_fwd_phylop.pdf", width = 2.5, height = 2.6)
plot(1:ncol(no_as_phylop), colMeans(no_as_phylop), type = "l", lwd=2,
  col = "grey50", ylim = c(.1, .8), bty = "n", ylab="PhyloP")
lines(1:ncol(no_as_phylop), colMeans(prox_as_phylop), lwd=2, col="red")
lines(1:ncol(no_as_phylop), colMeans(dist_as_phylop), lwd=2,col="black")
dev.off()

pdf("Figure5/Supp5e_antisense_rev_phylop.pdf", width = 2.5, height = 2.6)
plot(1:ncol(no_as_phylop), rev(colMeans(prox_as_phylop2)), type = "l", lwd=2,
  col = "red", ylim = c(0.1, .8),bty="n", ylab="PhyloP",xaxt="n")
lines(1:ncol(no_as_phylop), rev(colMeans(dist_as_phylop2)), lwd=2,col="black")
dev.off()

rm(phylop)

# Supplementary Figure 5f ------------------------------------------------------

phast <- import(PHASTCONS_BW)

no_as_phast <- normalizeToMatrix(phast, tss1_a, w = 10,
  extend = c(1000, 1000), value_column = "score")
prox_as_phast <- normalizeToMatrix(phast, tss1_b, w = 10,
  extend = c(1000, 1000), value_column = "score")
dist_as_phast <- normalizeToMatrix(phast, tss1_c, w = 10,
  extend = c(1000, 1000), value_column = "score")

prox_as_phast2 <- normalizeToMatrix(phast, tss1_d, w = 10,
  extend = c(1000, 1000), value_column = "score")
dist_as_phast2 <- normalizeToMatrix(phast, tss1_e, w = 10,
  extend = c(1000, 1000), value_column = "score")

pdf("Figure5/Supp5f_antisense_fwd_phastcons.pdf", width = 2.5, height = 2.6)
plot(1:ncol(no_as_phast), colMeans(no_as_phast), type = "l", lwd=2,
  col = "grey50", ylim = c(.1, .6), bty = "n", ylab="PhastCons")
lines(1:ncol(no_as_phast), colMeans(prox_as_phast), lwd=2, col="red")
lines(1:ncol(no_as_phast), colMeans(dist_as_phast), lwd=2,col="black")
dev.off()

pdf("Figure5/Supp5f_antisense_rev_phastcons.pdf", width = 2.5, height = 2.6)
plot(1:ncol(no_as_phast), rev(colMeans(prox_as_phast2)), type = "l", lwd=2,
  col = "red", ylim = c(0.1, .6),bty="n", ylab="PhastCons",xaxt="n")
lines(1:ncol(no_as_phast), rev(colMeans(dist_as_phast2)), lwd=2,col="black")
dev.off()

rm(phast)

# Supplementary Figure 5g-j ----------------------------------------------------

get_freqs <- function(x) {
  list(
    A = apply(x, 2, function(x) sum(x=="A")) / nrow(x),
    C = apply(x, 2, function(x) sum(x=="C")) / nrow(x),
    G = apply(x, 2, function(x) sum(x=="G")) / nrow(x),
    T = apply(x, 2, function(x) sum(x=="T")) / nrow(x)
  )
}

pc_dist <- final_as$PcTSS[!final_as$Proximal]
nc_dist <- final_as$AntisenseTSS[!final_as$Proximal]
pc_prox <- final_as$PcTSS[final_as$Proximal]
nc_prox <- final_as$AntisenseTSS[final_as$Proximal]

tss1_seqs <- tss1[tss1$name %in% c(pc_dist, nc_dist, pc_prox, nc_prox)]
tss1fa2kb <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, promoters(tss1_seqs, 500, 500))
names(tss1fa2kb) <- tss1_seqs$name
tss1fa2kbM <- as.matrix(t(as.data.frame(strsplit(as.character(tss1fa2kb), "", fixed=T))))
rownames(tss1fa2kbM) <- names(tss1fa2kb)

base_freqs <- list(
  Supp5g_tss_antisense_pc_distal_base_freq = get_freqs(tss1fa2kbM[pc_dist,,drop=F]),
  Supp5h_tss_antisense_nc_distal_base_freq = get_freqs(tss1fa2kbM[nc_dist,,drop=F]),
  Supp5i_tss_antisense_pc_proximal_base_freq = get_freqs(tss1fa2kbM[pc_prox,,drop=F]),
  Supp5j_tss_antisense_nc_proximal_base_freq = get_freqs(tss1fa2kbM[nc_prox,,drop=F])
)
for (i in names(base_freqs)) {
  f <- base_freqs[[i]]
  pdf(paste0("Figure5/", i, ".pdf"), width = 3.8, height = 3.5)
  plot(c((-499):0, 1:500), f$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
    ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
    col = universalmotif:::DNA_COLOURS["A"], bty = "n")
  axis(1, c(-500, 0, 500))
  axis(2, c(0, .5, 1))
  lines(c((-499):0, 1:500), f$C, lwd = 2,
    col = universalmotif:::DNA_COLOURS["C"])
  lines(c((-499):0, 1:500), f$G, lwd = 2,
    col = universalmotif:::DNA_COLOURS["G"])
  lines(c((-499):0, 1:500), f$T, lwd = 2,
    col = universalmotif:::DNA_COLOURS["T"])
  legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
    legend = c("A", "C", "G", "T"), bty="n")
  dev.off()
}

# Supplementary Figure 6a ------------------------------------------------------

atac_57h <- import(ATAC_57H_BW)

atac_57h_pc_d2 <- get_unstranded_signal(atac_57h, promoters(tss1_pc[!tss1_pc$Proximal],2000,1000),
  row.norm=F)
atac_57h_pc_p2 <- get_unstranded_signal(atac_57h, promoters(tss1_pc[tss1_pc$Proximal],2000,1000),
  row.norm=F)
atac_57h_pc_non_as2 <- get_unstranded_signal(atac_57h, promoters(tss1_pc_non_as2,2000,1000),
  row.norm=F)

atac_57h_as_d2 <- get_unstranded_signal(atac_57h, promoters(tss1_as[!tss1_as$Proximal],2000,1000),
  row.norm=F)
atac_57h_as_p2 <- get_unstranded_signal(atac_57h, promoters(tss1_as[tss1_as$Proximal],2000,1000),
  row.norm=F)

pdf("Figure5/Supp6a_atac_pc.pdf", width = 3.0, height = 3)
plot(atac_57h_pc_d2$Coord,atac_57h_pc_d2$Avg,type="l",lwd=2,ylim=c(0.9,2.1),
  bty="n",ylab=NA,xlab=NA,xaxt="n",yaxt="n"
)
lines(atac_57h_pc_p2$Coord,atac_57h_pc_p2$Avg,col="red",lwd=2)
lines(atac_57h_pc_non_as2$Coord,atac_57h_pc_non_as2$Avg,col="grey50",lwd=2)
axis(2, c(1, 1.5, 2))
axis(1, c(1, 1000, 2000, 3000), c(-2000, -1000, 0, 1000))
dev.off()

pdf("Figure5/Supp6a_atac_as.pdf", width = 3.0, height = 3)
plot(atac_57h_as_d2$Coord,rev(atac_57h_as_d2$Avg),type="l",lwd=2,ylim=c(0.9,2.1),
  bty="n",ylab=NA,xlab=NA,xaxt="n",yaxt="n"
)
lines(atac_57h_as_p2$Coord,rev(atac_57h_as_p2$Avg),col="red",lwd=2)
axis(1, c(1, 1000, 2000, 3000), rev(c(-2000, -1000, 0, 1000)))
dev.off()

rm(atac_57h)

# Supplementary Figure 6b ------------------------------------------------------

h3k4me3 <- import(H3K4ME3_BW)

h3k4me3_pc_non_as <- normalizeToMatrix(h3k4me3,
  tss1_pc_non_as2,
  c(2000,2000),w=10,value_column="score")
h3k4me3_pc_p <- normalizeToMatrix(h3k4me3,
  tss1_pc[tss1_pc$Proximal],
  c(2000,2000),w=10,value_column="score")
h3k4me3_pc_d <- normalizeToMatrix(h3k4me3,
  tss1_pc[!tss1_pc$Proximal],
  c(2000,2000),w=10,value_column="score")

h3k4me3_as_p <- normalizeToMatrix(h3k4me3,
  tss1_as[tss1_as$Proximal],
  c(2000,2000),w=10,value_column="score")
h3k4me3_as_d <- normalizeToMatrix(h3k4me3,
  tss1_as[!tss1_as$Proximal],
  c(2000,2000),w=10,value_column="score")

pdf("Figure5/Supp6b_h3k4me3_pc.pdf", width = 2.5, height = 2.7)
plot(1:ncol(h3k4me3_pc_p),zm(colMeans(h3k4me3_pc_d)), type = "l", lwd = 2,
  ylim = c(3, 25),
  bty = "n", ylab = NA, xlab = NA, xaxt = "n", yaxt = "n"
)
lines(1:ncol(h3k4me3_pc_p),zm(colMeans(h3k4me3_pc_non_as)), col = "grey50", lwd = 2)
axis(1, c(1, ncol(h3k4me3_pc_p)/2,ncol(h3k4me3_pc_p)),
  c(-2, 0, 2))
axis(2, c(5, 15, 25))
lines(1:ncol(h3k4me3_pc_p),zm(colMeans(h3k4me3_pc_p)), col = "red", lwd = 2)
dev.off()

pdf("Figure5/Supp6b_h3k4me3_as.pdf", width = 2.5, height = 2.7)
plot(1:ncol(h3k4me3_pc_p),rev(zm(colMeans(h3k4me3_as_d))), type = "l", lwd = 2,
  ylim = c(3, 25), bty = "n", ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
lines(1:ncol(h3k4me3_pc_p),rev(zm(colMeans(h3k4me3_as_p))), col = "red", lwd = 2)
axis(1, c(1, ncol(h3k4me3_pc_p)/2,ncol(h3k4me3_pc_p)),
  c(-2, 0, 2))
axis(2, c(5, 15, 25))
dev.off()

rm(h3k4me3)

# Supplementary Figure 6c ------------------------------------------------------

h3k9ac <- import(H3K9AC_BW)

h3k9ac_pc_non_as <- normalizeToMatrix(h3k9ac,
  tss1_pc_non_as2,
  c(2000,2000),w=10,value_column="score")
h3k9ac_pc_p <- normalizeToMatrix(h3k9ac,
  tss1_pc[tss1_pc$Proximal],
  c(2000,2000),w=10,value_column="score")
h3k9ac_pc_d <- normalizeToMatrix(h3k9ac,
  tss1_pc[!tss1_pc$Proximal],
  c(2000,2000),w=10,value_column="score")

h3k9ac_as_p <- normalizeToMatrix(h3k9ac,
  tss1_as[tss1_as$Proximal],
  c(2000,2000),w=10,value_column="score")
h3k9ac_as_d <- normalizeToMatrix(h3k9ac,
  tss1_as[!tss1_as$Proximal],
  c(2000,2000),w=10,value_column="score")

pdf("Figure5/Supp6c_h3k9ac_pc.pdf", width = 2.5, height = 2.7)
plot(1:ncol(h3k9ac_pc_p),zm(colMeans(h3k9ac_pc_d)), type = "l", lwd = 2,
  ylim = c(5, 55),
  bty = "n", ylab = NA, xlab = NA, xaxt = "n", yaxt = "n"
)
lines(1:ncol(h3k9ac_pc_p),zm(colMeans(h3k9ac_pc_non_as)), col = "grey50", lwd = 2)
axis(1, c(1, ncol(h3k9ac_pc_p)/2,ncol(h3k9ac_pc_p)),
  c(-2, 0, 2))
axis(2, c(5, 25, 55))
lines(1:ncol(h3k9ac_pc_p),zm(colMeans(h3k9ac_pc_p)), col = "red", lwd = 2)
dev.off()

pdf("Figure5/Supp6c_h3k9ac_as.pdf", width = 2.5, height = 2.7)
plot(1:ncol(h3k9ac_pc_p),rev(zm(colMeans(h3k9ac_as_d))), type = "l", lwd = 2,
  ylim = c(5, 55), bty = "n", ylab = NA, xlab = NA, xaxt = "n", yaxt = "n")
lines(1:ncol(h3k9ac_pc_p),rev(zm(colMeans(h3k9ac_as_p))), col = "red", lwd = 2)
axis(1, c(1, ncol(h3k9ac_pc_p)/2,ncol(h3k9ac_pc_p)),
  c(-2, 0, 2))
axis(2, c(5, 25, 55))
dev.off()

rm(h3k9ac)

# Supplementary Figure 6d ------------------------------------------------------

tss1_pc$cv_csrna <- rowSds(tssM[tss1_pc$name,])/rowMeans(tssM[tss1_pc$name,])
noas_pc_cv <- tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding" & !tssA$TSS %in% tss1_pc$name]
noas_pc_cv <- rowSds(tssM[noas_pc_cv,])/rowMeans(tssM[noas_pc_cv,])

cv_df <- data.frame(
  Which = c(rep("NoAS", length(noas_pc_cv)), tss1_pc$Proximal),
  CV = c(noas_pc_cv, tss1_pc$cv_csrna))

rstatix::wilcox_test(cv_df, CV~Which)

ggplot(cv_df, aes(CV, colour = Which)) +
  geom_density(bw = .2, size = .75) +
  xlab("CV (csRNA-seq)") +
  ylab("Density") +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8)], name = element_blank()) +
  theme_5pt
ggsave("Figure5/Supp6d_antisense_cv.pdf", width = 2.75, height = 1.16)

# Supplementary Figure 6e-g ----------------------------------------------------

g_ends <- a11_g
start(g_ends) <- end(g_ends)
end(g_ends)[as.character(strand(a11_g))=="-"] <- start(g_ends)[as.character(strand(a11_g))=="-"] <-
  start(a11_g)[as.character(strand(a11_g))=="-"]

g_ends_as_p <- g_ends[g_ends$gene_id %in% final_as$PcGene[final_as$Proximal]]
g_ends_as_d <- g_ends[g_ends$gene_id %in% final_as$PcGene[!final_as$Proximal]]
g_ends_non <- g_ends[g_ends$gene_id %in% pc_bed$name & !g_ends$gene_id %in% final_as$PcGene]
g_ends_non <- g_ends_non[g_ends_non$gene_id %in% tssA$NearestGene]

g_p <- a11_g[a11_g$gene_id %in% g_ends_as_p$gene_id]
g_d <- a11_g[a11_g$gene_id %in% g_ends_as_d$gene_id]
g_n <- a11_g[a11_g$gene_id %in% g_ends_non$gene_id]

rnapii_signals <- list(
  Supp6e_antisense_rnapii = RNAPII_BW,
  Supp6f_antisense_rnapii_ser5p = SER5P_BW,
  Supp6g_antisense_rnapii_ser2p = SER2P_BW
)
for (i in names(rnapii_signals)) {
  rnapii <- import(rnapii_signals[[i]])

  a_p <- normalizeToMatrix(rnapii, g_p, w = 20,
    target_ratio = .6, value_column = "score")
  a_d <- normalizeToMatrix(rnapii, g_d, w = 20,
    target_ratio = .6, value_column = "score")
  a_n <- normalizeToMatrix(rnapii, g_n, w = 20,
    target_ratio = .6, value_column = "score")

  pdf(paste0("Figure5/", i, ".pdf"), width = 3, height = 2.6)
  plot(1:ncol(a_p), colMeans(a_n), type = "l", lwd=3, col = "grey50", bty = "n",
    ylim = c(0, 1.1*max(c(colMeans(a_n), colMeans(a_p), colMeans(a_d)))), xaxt = "n")
  lines(1:ncol(a_d), colMeans(a_p), lwd=3, col="red")
  lines(1:ncol(a_n), colMeans(a_d), lwd=3,col="black")
  axis(1, c(1, which(colnames(a_p)=="t1"),which(colnames(a_p)=="d1")-1,ncol(a_p)),
    c("-5kbp", "TSS", "TTS", "+5kbp"))
  dev.off()
}

rm(rnapii)

# Supplementary Figure 6h ------------------------------------------------------

no_as_map1 <- get_unstranded_signal(mnase3, no_as_win, keep = c(.00, .95),
  return.matrix = TRUE, row.norm = F)
no_as_map2 <- get_unstranded_signal(mnase3, prox_as_win, keep = c(.00, .95),
  return.matrix = TRUE, row.norm = F)
no_as_map3 <- get_unstranded_signal(mnase3, dist_as_win, keep = c(.00, .95),
  return.matrix = TRUE, row.norm = F)

a_none <- colMeans(no_as_map1)
b_proxi <- colMeans(no_as_map2)
c_distal <- colMeans(no_as_map3)

pdf("Figure5/Supp6h_mnase3.pdf", width = 2.6, height = 2.6)
plot(1:1500, c_distal[1:1500], type = "l", ylim = c(0, 60),
  lwd=2,bty="n")
lines(1:1500, b_proxi[1:1500], col = "red",lwd=2)
lines(1:1500, a_none[1:1500], col = "grey50",lwd=2)
dev.off()

rm(mnase3)

# Supplementary Figure 6i ------------------------------------------------------

table(tss1_pc$cls)

cc <- as.character(0:6)

exp_obs_plot <- function(sel) {
  ggplot(data.frame(
    Cluster = rep(paste0("C", 0:6), 2),
    Count = c(
      round(table(tss1_pc$cls) * (sum(sel) / length(tss1_pc)))[cc],
      round(table(tss1_pc$cls[sel]))[cc]),
    What = c(rep("Expected", 7), rep("Observed", 7))
  ), aes(Cluster, Count, fill = What)) +
    geom_bar(stat = "identity", width = 0.75, position = "dodge") +
    scale_fill_manual(values = c(colorBlindBlack8[1:2]), name = element_blank()) +
    scale_y_continuous(expand = c(0, 0), name = "# of pcTSSs") +
    xlab(element_blank()) +
    theme_5pt
}

# Chi-squared tests of observed versus expected counts per cluster
exp_obs_test <- function(sel, n = sum(tss1_pc$cls == 6)) {
  ex <- table(tss1_pc$cls) * (sum(sel) / length(tss1_pc))
  ob <- table(factor(tss1_pc$cls[sel], levels = names(ex)))
  sapply(as.character(1:6), function(i) {
    chisq.test(matrix(c(ob[i], n - ob[i], ex[i], n - ex[i]), ncol = 2))$p.value
  })
}

corr_sets <- list(
  "PCC > 0.25" = tss1_pc$AntiCorr >= .25,
  "PCC > 0.5" = tss1_pc$AntiCorr >= .5,
  "PCC > 0.75" = tss1_pc$AntiCorr >= .75,
  "PCC < -0.25" = tss1_pc$AntiCorr <= -.25,
  "PCC < -0.5" = tss1_pc$AntiCorr <= -.5,
  "PCC < -0.75" = tss1_pc$AntiCorr <= -.75
)

lapply(corr_sets, exp_obs_test)

cowplot::save_plot(
  "Figure5/Supp6i_expected_vs_observed_per_cluster.pdf",
  cowplot::plot_grid(plotlist = lapply(corr_sets, exp_obs_plot), nrow = 1),
  base_width = 13, base_height = 1)

# Supplementary Figure 6j, k ---------------------------------------------------

# PIN4
plot_tss_expr("TSS2_183328", 300, "Figure5/Supp6j_pin4_expr_sense.pdf")
plot_tss_expr("TSS2_181574", 15, "Figure5/Supp6j_pin4_expr_antisense.pdf")

# RD26
plot_tss_expr("TSS4_13709122", 340, "Figure5/Supp6k_rd26_expr_sense.pdf")
plot_tss_expr("TSS4_13707964", 10, "Figure5/Supp6k_rd26_expr_antisense.pdf")
