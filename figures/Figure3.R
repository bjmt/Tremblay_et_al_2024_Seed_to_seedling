# Benjamin J.M. Tremblay
#
# Figure 3 and Supplementary Figure 3
#
# Package requirements: universalmotif, matrixStats, Biostrings, rtracklayer,
# EnrichedHeatmap, ComplexHeatmap, BSgenome.Athaliana.TAIR.TAIR9,
# GenomicFeatures, ggplot2, dplyr, rstatix, reshape2, cowplot, circlize,
# viridisLite, readr, trackvizir, genomation

library(universalmotif)
library(matrixStats)
library(Biostrings)
library(rtracklayer)
library(EnrichedHeatmap)
library(ComplexHeatmap)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(GenomicFeatures)
library(ggplot2)
library(dplyr)
library(trackvizir)

source("Functions.R")

seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

dir.create("Figure3", showWarnings = FALSE)
dir.create("Figure3/motif_tss_dens", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

MOTIFS_MEME <- "all_final.meme"
TOMTOM_TSV <- "tomtom_final_motifs.tsv"
MOTIF_HITS_BED <- "final_streme_motifs.bed.gz"
ENRICHMENT_RDS <- "enrichment_matrix_final.RDS"
CLUSTER_FASTA_DIR <- "fasta2"
TSS_BED <- "tss.size1.bed"
TSS_500BP_BED <- "tss.size1_u400d100.bed"
ATAC_500BP_BED <- "peaks_final_500bp.bed"
TSS_CLUSTERS <- "csRNAseq_clustering/Clusters_anno.tsv"
ATAC_CLUSTERS <- "ATACseq_clustering3/Clusters_anno.tsv"
ARAPORT11_GFF <- "Araport11_novel.gff3"
ABI5_PEAKS_BED <- "ABI5.bed"
TAIR10_FASTA <- "genome.fa"

# Read density
ATAC_DS_BW <- "atacseq_bw/DS.bw"
ATAC_6H_BW <- "atacseq_bw/L6.bw"
ATAC_26H_BW <- "atacseq_bw/L26.bw"
ATAC_57H_BW <- "atacseq_bw/L57.bw"
ABI5_DAP_BW <- "ABI5-SRX670509.bw"
MNASE_BW <- "Ath_leaf_NPS.bw"
CSRNA_ALL_PLUS_BG <- "csrnaseq_bedGraphs/all.plus.bedGraph.gz"
CSRNA_ALL_MINUS_BG <- "csrnaseq_bedGraphs/all.minus.bedGraph.gz"

# Load motifs, motif instances and clusters ------------------------------------

mots <- read_meme(MOTIFS_MEME)
names(mots) <- sapply(mots, function(x) x["name"])

# Same motifs but named by their M number
mots0 <- mots
names(mots0) <- structure(to_df(mots)$altname, names = to_df(mots)$name)[names(mots0)]
switch_name <- function(x) { a=x@name;b=x@altname;x@name=b;x@altname=a;x }
for (i in seq_along(mots0)) mots0[[i]] <- switch_name(mots0[[i]])

# Best matching known motifs, with the M numbers added
tom_all <- readr::read_tsv(TOMTOM_TSV)
tom_all <- cbind(alt = structure(names(mots0), names = names(mots))[tom_all[[1]]], tom_all)
readr::write_tsv(tom_all, "Figure3/tomtom_final_motifs2.tsv")

all_mots <- import(MOTIF_HITS_BED,
  extraCols = c(raw.score = "numeric", pct.score = "numeric", pvalue = "numeric",
    qvalue = "numeric"))
all_mots$alt <- structure(names(mots0), names = names(mots))[all_mots$name]

tss1 <- import(TSS_BED)

tssA <- readr::read_tsv(TSS_CLUSTERS)
atacA <- readr::read_tsv(ATAC_CLUSTERS)

tssM <- as.matrix(tssA[, 11:22])
rownames(tssM) <- tssA$TSS
tssMm <- as.matrix(data.frame(row.names = rownames(tssM),
  DS = rowMeans(tssM[, 1:2]),
  S24 = rowMeans(tssM[, 3:4]),
  S72 = rowMeans(tssM[, 5:6]),
  L6 = rowMeans(tssM[, 7:8]),
  L26 = rowMeans(tssM[, 9:10]),
  L57 = rowMeans(tssM[, 11:12])
))

atacM <- as.matrix(atacA[, 13:20])
rownames(atacM) <- atacA$PeakID
atacMm <- as.matrix(data.frame(row.names = rownames(atacM),
  DS = rowMeans(atacM[, 1:2]),
  L6 = rowMeans(atacM[, 3:4]),
  L26 = rowMeans(atacM[, 5:6]),
  L57 = rowMeans(atacM[, 7:8])
))

get_gr <- function(x, y) {
  z <- y[y$name %in% x[[1]]]
  z$group <- structure(x[[2]], names = x[[1]])[z$name]
  z
}

# TSSs used for motif enrichment
fa2k <- lapply(list.files(CLUSTER_FASTA_DIR, full = TRUE), readDNAStringSet)
names(fa2k) <- gsub("[.]fa", "", list.files(CLUSTER_FASTA_DIR))
tss2k <- fa2k[c("tss_C1", "tss_C2", "tss_C3", "tss_C4", "tss_C5", "tss_C6")]
names(tss2k) <- gsub("tss_", "", names(tss2k))

tss2kgr <- get_gr(
  do.call(rbind, lapply(names(tss2k), function(x) data.frame(names(tss2k[[x]]), x))),
  tss1)

tss500 <- import(TSS_500BP_BED)
atac500 <- import(ATAC_500BP_BED)
seqlengths(atac500) <- seqlengths(BSgenome.Athaliana.TAIR.TAIR9)[1:5]
atac500 <- trim(atac500)
atac500 <- atac500[width(atac500) == 500]

enr <- readRDS(ENRICHMENT_RDS)

# Figure 3a --------------------------------------------------------------------

enr2 <- cbind(enr[[1]], enr[[2]])
enr2 <- enr2[rowMaxs(enr[[1]]) > 10 & rowMaxs(enr[[2]]) > 10, ]
enr2 <- enr2[c(1, 3, 4, 6, 5, 7, 2), ]

enr0 <- cbind(enr[[1]], enr[[2]])
rownames(enr0) <- structure(to_df(mots)$altname, names = to_df(mots)$name)[rownames(enr0)]

enr22 <- enr2
rownames(enr22) <- structure(names(mots0), names = names(mots))[rownames(enr2)]

p <- view_motifs(mots0[rownames(enr22)], names.pos = "r", use.type = "ICM",
  min.overlap = 99, tryRC = FALSE)
ggsave("Figure3/Fig3a_motif_logos.pdf", width = 3, height = 3, plot = p)

pdf("Figure3/Fig3a_motif_heatmap.pdf", width = 2.9, height = 1.6)
Heatmap(enr22, cluster_rows = FALSE, cluster_columns = FALSE,
  name = "-log10(Adj. P-value)",
  col = circlize::colorRamp2(0:30, viridisLite::magma(31)),
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_rot = 0,
  column_gap = unit(0.2, "mm"),
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6)),
  column_split = c(rep("A", 5), rep("C", 6))
)
dev.off()

# Figure 3a, Supplementary Figure 3a, b ----------------------------------------

# Positional density of each motif around the TSSs

make_mot_dens <- function(m, gr, extend = 1000, colMeans = TRUE) {
  s <- all_mots[all_mots$alt %in% m]
  if (colMeans) {
    lapply(m, function(x) colMeans(normalizeToMatrix(s[s$alt == x], gr, extend = extend)))
  } else {
    lapply(m, function(x) normalizeToMatrix(s[s$alt == x], gr, extend = extend)[])
  }
}
for (i in 1:nrow(enr0)) {
  message(rownames(enr0)[i])
  d <- make_mot_dens(rownames(enr0)[i], tss2kgr, colMeans = FALSE)[[1]]
  pdf(width = 6, height = 3.5, paste0("Figure3/motif_tss_dens/", rownames(enr0)[i], ".pdf"))
  plot(xaxt = "n", lwd = 2, density(bw = 2,
    unlist(sapply(1:ncol(d), function(x) rep(x, sum(d[, x]))), use.names = FALSE)
  ), xlab = NA, main = paste0(rownames(enr0)[i], " (n=", sum(d), ")"))
  axis(1, c(1, ncol(d) / 4, ncol(d) / 2, ncol(d) * (3/4), ncol(d)),
    c("-1000", "-500", 0, "+500", "+1000"))
  abline(v = ncol(d) / 2, lty = 2, lwd = 2, col = "grey50")
  dev.off()
}

# Supplementary Figure 3a ------------------------------------------------------

enr0a <- enr0[!rownames(enr0) %in% rownames(enr22), ]
enr0a <- enr0a[rowMaxs(enr0a[,1:5])>3 & rowMaxs(enr0a[,6:11])<=3, ]
enr0a <- enr0a[c("M11", "M3", "M32", "M7", "M31", "M8", "M10", "M12"),]

pdf("Figure3/SuppFig3a_motif_heatmap_atac.pdf", width = 2.2, height = 1.6)
Heatmap(enr0a[,1:5], cluster_rows = FALSE, cluster_columns = FALSE,
  name = "-log10(Adj. P-value)",
  col = circlize::colorRamp2(0:16, viridisLite::magma(17)),
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_rot = 0,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6))
)
dev.off()

p <- view_motifs(mots0[rownames(enr0a)], names.pos = "r", use.type = "ICM",
  min.overlap = 99, tryRC = FALSE)
ggsave("Figure3/SuppFig3a_motif_logos_atac.pdf", width = 3, height = 3.2, plot = p)

# Supplementary Figure 3b ------------------------------------------------------

enr0c <- enr0[!rownames(enr0) %in% rownames(enr22), ]
enr0c <- enr0c[rowMaxs(enr0c[,6:11])>10, ]

pdf("Figure3/SuppFig3b_motif_heatmap_csrna.pdf", width = 2.4, height = 1.45)
Heatmap(enr0c[,6:11], cluster_rows = FALSE, cluster_columns = FALSE,
  name = "-log10(Adj. P-value)",
  col = circlize::colorRamp2(0:16, viridisLite::magma(17)),
  column_names_gp = gpar(fontsize = 6),
  row_names_gp = gpar(fontsize = 6),
  column_names_rot = 0,
  heatmap_legend_param = list(title_gp = gpar(fontsize = 6))
)
dev.off()

p <- view_motifs(mots0[rownames(enr0c)], names.pos = "r", use.type = "ICM",
  min.overlap = 99, tryRC = FALSE)
ggsave("Figure3/SuppFig3b_motif_logos_csrna.pdf", width = 3, height = 3, plot = p)

# Figure 3b, c -----------------------------------------------------------------

ABI5_peaks <- import(ABI5_PEAKS_BED)

ABI5_targets <- tssA$TSS[tssA$Cluster == 1]
ABI5_targets_v2 <- rownames(tssMm[ABI5_targets,])[tssMm[ABI5_targets, 1] > 50]
ABI5_targets_v2_motif <- ABI5_targets_v2[ABI5_targets_v2 %in%
  tss500$name[overlapsAny(tss500, ABI5_peaks)]]
ABI5_v2_atac_motif <- atac500$name[overlapsAny(atac500, tss1[tss1$name %in% ABI5_targets_v2_motif])]

M2_targets <- tssA$TSS[tssA$Cluster %in% 4:5]
M2_targets <- rownames(tssMm[M2_targets,])[rowMaxs(tssMm[M2_targets, 4:5]) > 50]
M2_targets_motif <- M2_targets[M2_targets %in%
  tss500$name[overlapsAny(tss500, all_mots[all_mots$alt == "M2"])]]
M2_motif_atac <- atac500$name[overlapsAny(atac500, tss1[tss1$name %in% M2_targets_motif])]

melt2 <- function(x, colnames = NULL) {
  x <- reshape2::melt(x)
  x[[1]] <- as.character(x[[1]])
  x[[2]] <- as.character(x[[2]])
  if (!is.null(colnames) && length(colnames) == 3) {
    colnames(x) <- colnames
  }
  x
}

abi5_csrna_targets <- melt2(t(scale(t(tssMm[ABI5_targets_v2_motif, ]))), c("TSS", "Sample", "z-score"))
abi5_atac_targets <- melt2(t(scale(t(atacMm[ABI5_v2_atac_motif, ]))), c("Peak", "Sample", "z-score"))
m2_csrna_targets <- melt2(t(scale(t(tssMm[M2_targets_motif,]))),c("TSS","Sample","z-score"))
m2_atac_targets <- melt2(t(scale(t(atacMm[M2_motif_atac,]))),c("Peak", "Sample","z-score"))

abi5_csrna_targets$Sample <- factor(abi5_csrna_targets$Sample, levels = colnames(tssMm))
abi5_atac_targets$Sample <- factor(abi5_atac_targets$Sample, levels = colnames(atacMm))
m2_csrna_targets$Sample <- factor(m2_csrna_targets$Sample, levels = colnames(tssMm))
m2_atac_targets$Sample <- factor(m2_atac_targets$Sample, levels = colnames(atacMm))

# Pairwise two-sided Mann-Whitney tests
abi5_csrna_targets[,-1] %>%
  rename(zscore = `z-score`) %>%
  rstatix::wilcox_test(zscore~Sample)

abi5_atac_targets[,-1] %>%
  rename(zscore = `z-score`) %>%
  rstatix::wilcox_test(zscore~Sample)

m2_atac_targets[,-1] %>%
  rename(zscore = `z-score`) %>%
  rstatix::wilcox_test(zscore~Sample)

m2_csrna_targets[,-1] %>%
  rename(zscore = `z-score`) %>%
  rstatix::wilcox_test(zscore~Sample)

p1 <- ggplot(abi5_csrna_targets, aes(Sample, `z-score`)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, fill = "grey75", size = .3, width = .2, colour = "black") +
  coord_cartesian(ylim = c(-2, 2)) +
  ylab("Z-score") +
  xlab(element_blank()) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, colour = "grey50") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p2 <- ggplot(abi5_atac_targets, aes(Sample, `z-score`)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, fill = "grey75", width = .2, size = .3) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, colour = "grey50") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p3 <- ggplot(m2_csrna_targets, aes(Sample, `z-score`)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, fill = "grey75", width = .2, size = .3) +
  coord_cartesian(ylim = c(-2, 2)) +
  ylab("Z-score") +
  xlab(element_blank()) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, colour = "grey50") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

p4 <- ggplot(m2_atac_targets, aes(Sample, `z-score`)) +
  geom_violin() +
  geom_boxplot(outlier.shape = NA, fill = "grey75", width = .2, size = .3) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  xlab(element_blank()) +
  ylab(element_blank()) +
  geom_hline(yintercept = 0, lwd = 0.5, lty = 2, colour = "grey50") +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )

cowplot::save_plot("Figure3/Fig3bc_m1_m2_z.pdf", base_width = 2.75, base_height = 1.65,
cowplot::plot_grid(p1, p2, p3, p4,
  align = "hv",
  rel_widths = c(1.25, 1, 1.25, 1),
  ncol = 2, byrow = TRUE
)
)

# Figure 3d --------------------------------------------------------------------

# Gene, TSS and promoter ACR IDs
ABI5 <- c("AT2G36270", "TSS2_15207491", "ATAC2_15207685")
EM1 <- c("AT3G51810", "TSS3_19214746", "ATAC3_19214649")
EM6 <- c("AT2G40170", "TSS2_16780256", "ATAC2_16780364")

txdb <- makeTxDbFromGFF(ARAPORT11_GFF)
abi5_gr <- GRanges("2", IRanges(15203928, 15208926))
abi5_hl <- GRanges("2", IRanges(15207433, 15207899))
em1_gr <- GRanges("3", IRanges(19214088-1250, 19216586+1250))
em1_hl <- GRanges("3", IRanges(19214509, 19214780))
em6_gr <- GRanges("2", IRanges(16778975-1500, 16780973+1500))
em6_hl <- GRanges("2", IRanges(16780151, 16780727))

tp <- trackvizir(abi5_gr, highlight.ranges = abi5_hl, compact = TRUE, y.axis.hide=T) +
  track_area(read_bw(ATAC_DS_BW, ranges = abi5_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_6H_BW, ranges = abi5_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_26H_BW, ranges = abi5_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_57H_BW, ranges = abi5_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ABI5_DAP_BW, ranges = abi5_gr), ylims = c(0, 75), bin.size=10,colour.pos="forestgreen") +
  track_anno(txdb) +
  track_coords()
save_tracks(build_tracks(tp), "Figure3/Fig3d_ABI5_browser.pdf", base_width = 2.0, base_height = 1.9)

tp <- trackvizir(em1_gr, highlight.ranges = em1_hl, compact = TRUE, y.axis.hide=T) +
  track_area(read_bw(ATAC_DS_BW, ranges = em1_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_6H_BW, ranges = em1_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_26H_BW, ranges = em1_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_57H_BW, ranges = em1_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ABI5_DAP_BW, ranges = em1_gr), ylims = c(0, 75), bin.size=10,colour.pos="forestgreen") +
  track_anno(txdb) +
  track_coords()
save_tracks(build_tracks(tp), "Figure3/Fig3d_EM1_browser.pdf", base_width = 2.2, base_height = 1.9)

tp <- trackvizir(em6_gr, highlight.ranges = em6_hl, compact = TRUE, y.axis.hide=T) +
  track_area(read_bw(ATAC_DS_BW, ranges = em6_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_6H_BW, ranges = em6_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_26H_BW, ranges = em6_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ATAC_57H_BW, ranges = em6_gr), ylims = c(0, 30), bin.size=10,colour.pos="grey25") +
  track_area(read_bw(ABI5_DAP_BW, ranges = em6_gr), ylims = c(0, 75), bin.size=10,colour.pos="forestgreen") +
  track_anno(txdb) +
  track_coords()
save_tracks(build_tracks(tp), "Figure3/Fig3d_EM6_browser.pdf", base_width = 2.2, base_height = 1.9)

# csRNA-seq (black) and ATAC-seq (blue, right axis) quantification
abi5_expr <- tssM[ABI5[2],]
abi5_acc <- atacM[ABI5[3],]

em1_expr <- tssM[EM1[2],]
em1_acc <- atacM[EM1[3],]

em6_expr <- tssM[EM6[2],]
em6_acc <- atacM[EM6[3],]

pdf("Figure3/Fig3d_abi5_point_line.pdf", width = 3.8, height = 3.4)
plot(rep(1:6, each = 2)-.1, abi5_expr,
  pch = 19, bty = "n", xlim = c(0.8, 6.2), xlab = NA, ylab = NA,
  xaxt = "n", yaxt = "n"
)
lines(1:6 - .1, tapply(abi5_expr, rep(1:6, each = 2), mean), lwd = 1.5)
points(rep(c(1, 4, 5, 6), each = 2)+.1, abi5_acc*3, col = "blue", pch = 19)
lines(c(1, 4, 5, 6)+.1, tapply(abi5_acc, rep(1:4, each = 2), mean)*3, col = "blue", lwd=1.5)
axis(1, 1:6, colnames(tssMm))
axis(2, c(0, 300, 600))
axis(4, c(0, 300, 600), c(0, 100, 200))
dev.off()

pdf("Figure3/Fig3d_em1_point_line.pdf", width = 3.8, height = 3.4)
plot(rep(1:6, each = 2)-.1, em1_expr,
  pch = 19, bty = "n", xlim = c(0.8, 6.2), xlab = NA, ylab = NA,
  xaxt = "n", yaxt = "n", ylim = c(0, 8000)
)
lines(1:6 - .1, tapply(em1_expr, rep(1:6, each = 2), mean), lwd=1.5)
points(rep(c(1, 4, 5, 6), each = 2)+.1, em1_acc*260, col = "blue", pch = 19)
lines(c(1, 4, 5, 6)+.1, tapply(em1_acc, rep(1:4, each = 2), mean)*260, col = "blue", lwd=1.5)
axis(1, 1:6, colnames(tssMm))
axis(2, c(0, 4000, 8000))
axis(4, c(0, 4000, 8000), c(0, 15, 30))
dev.off()

pdf("Figure3/Fig3d_em6_point_line.pdf", width = 3.8, height = 3.4)
plot(rep(1:6, each = 2)-.1, em6_expr,
  pch = 19, bty = "n", xlim = c(0.8, 6.2), xlab = NA, ylab = NA,
  xaxt = "n", yaxt = "n", ylim = c(0, 25000)
)
lines(1:6 - .1, tapply(em6_expr, rep(1:6, each = 2), mean), lwd=1.5)
points(rep(c(1, 4, 5, 6), each = 2)+.1, em6_acc*180, col = "blue", pch = 19)
lines(c(1, 4, 5, 6)+.1, tapply(em6_acc, rep(1:4, each = 2), mean)*180, col = "blue", lwd=1.5)
axis(1, 1:6, colnames(tssMm))
axis(2, c(0, 12500, 25000))
axis(4, c(0, 12500, 25000), c(0, 70, 140))
dev.off()

# Supplementary Figure 3d, e, g ------------------------------------------------

# Average leaf MNase-seq, csRNA-seq and ATAC-seq signal around the motif instances

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

m <- read_meme(MOTIFS_MEME)

atac500_mnase <- import(ATAC_500BP_BED)

seqs <- lapply(list.files(CLUSTER_FASTA_DIR, full=TRUE), readDNAStringSet)
names(seqs) <- basename(list.files(CLUSTER_FASTA_DIR))
names(seqs)[names(seqs) == "atac_bkg.fa"] <- "Abkg"
names(seqs)[names(seqs) == "tss_bkg.fa"] <- "Cbkg"
names(seqs) <- gsub("[.]fa", "", gsub("atac_", "", gsub("tss_", "", names(seqs))))

seqs_A <- c(seqs[["A1"]], seqs[["A2"]], seqs[["A3"]], seqs[["A4"]], seqs[["A5"]])
seqs_C <- c(seqs[["C1"]], seqs[["C2"]], seqs[["C3"]], seqs[["C4"]], seqs[["C5"]], seqs[["C6"]])

tair <- Biostrings::readDNAStringSet(TAIR10_FASTA)
names(tair) <- paste0("Chr", c(1:5, "M", "C"))

gr_i <- c(tss500[tss500$name %in% names(seqs_C)], atac500_mnase[atac500_mnase$name %in% names(seqs_A)])
seqlevels(gr_i) <- paste0("Chr", 1:5)

csrna_p <- import(CSRNA_ALL_PLUS_BG)
csrna_m <- import(CSRNA_ALL_MINUS_BG)
strand(csrna_p) <- "+"
strand(csrna_m) <- "-"
seqlevels(csrna_p) <- paste0("Chr", c(1:5, "M", "C"))
seqlevels(csrna_m) <- paste0("Chr", c(1:5, "M", "C"))
csrna <- sort(c(csrna_p, csrna_m))
rm(csrna_p, csrna_m)

atac_ds <- import(ATAC_DS_BW)
seqlevels(atac_ds) <- paste0("Chr", 1:5)
atac_l6 <- import(ATAC_6H_BW)
seqlevels(atac_l6) <- paste0("Chr", 1:5)
atac_l26 <- import(ATAC_26H_BW)
seqlevels(atac_l26) <- paste0("Chr", 1:5)
atac_l57 <- import(ATAC_57H_BW)
seqlevels(atac_l57) <- paste0("Chr", 1:5)

mnase_panels <- c(M2 = "SuppFig3d", M23 = "SuppFig3e", M1 = "SuppFig3g")

for (i in which(sapply(m, function(x) x["altname"]) %in% names(mnase_panels))) {

  message(i)

  i_motif <- scan_sequences(m[[i]], tair,
    threshold = 0.0001, RC = TRUE, return.granges = TRUE)

  i_motif <- resize(i_motif[overlapsAny(i_motif, gr_i)], 2000, "center")

  i_mark <- get_unstranded_signal(MNASE_BW, i_motif, row.norm = FALSE)$Avg
  i_csrna <- get_stranded_signal(csrna, i_motif, row.norm=F, keep=c(0,.99))
  i_atac_ds <- get_unstranded_signal(atac_ds, i_motif, row.norm = FALSE)$Avg
  i_atac_l6 <- get_unstranded_signal(atac_l6, i_motif, row.norm = FALSE)$Avg
  i_atac_l26 <- get_unstranded_signal(atac_l26, i_motif, row.norm = FALSE)$Avg
  i_atac_l57 <- get_unstranded_signal(atac_l57, i_motif, row.norm = FALSE)$Avg

  pdf(paste0("Figure3/", mnase_panels[m[[i]]["altname"]], "_MNase_csRNA_ATAC_",
      m[[i]]["altname"], "_", m[[i]]["name"], ".pdf"),
    width = 6.5, height = 7.5)
  par(mfrow = c(3,1), mai = c(0.4, 0.8, 0.2, 0.1))
  plot(seq(-1000, 999), i_mark, ylim = c(0, max(i_mark)*1.1),
    type = "l", lwd = 3, xlab = NA, ylab = "MNase-seq", bty = "n",
    main = paste0("Motif ", m[[i]]["altname"]))
  abline(v=0, lty = 2, lwd = 2, col = "grey45")
  plot(seq(-1000, 999), i_csrna$PosAvg, col = "blue", type = "l",
    xlab = NA, ylab = "csRNA-seq", bty = "n", main = NA,
    ylim=c(0, max(c(i_csrna$PosAvg,i_csrna$NegAvg))*1.1))
  lines(seq(-1000, 999), i_csrna$NegAvg, col = "forestgreen")
  legend("topright", bty = "n", col = c("blue", "forestgreen"), lty = 1, lwd = 2,
    legend = c("Sense", "Antisense"))
  plot(seq(-1000, 999), i_atac_ds, type = "l",
    xlab = NA, ylab = "ATAC-seq", bty = "n", main = NA, col = colorBlindBlack8[1],
    ylim=c(0, max(c(i_atac_ds,i_atac_l6,i_atac_l26,i_atac_l57))*1.1))
  lines(seq(-1000, 999), i_atac_l6, col = colorBlindBlack8[4])
  lines(seq(-1000, 999), i_atac_l26, col = colorBlindBlack8[5])
  lines(seq(-1000, 999), i_atac_l57, col = colorBlindBlack8[6])
  legend("topright", bty = "n", col = colorBlindBlack8[c(1,4:6)], lty = 1, lwd = 2,
    legend = c("DS", "L6", "L26", "L57"))
  dev.off()

}

# Supplementary Figure 3f ------------------------------------------------------

# TSSs making up the top 50% of the cumulative csRNA-seq signal per time-point
tss1$CS_DS <- cumsum(sort(tssMm[,"DS"] / sum(tssMm[,"DS"])))[tss1$name]
tss1$CS_L6 <- cumsum(sort(tssMm[,"L6"] / sum(tssMm[,"L6"])))[tss1$name]
tss1$CS_L26 <- cumsum(sort(tssMm[,"L26"] / sum(tssMm[,"L26"])))[tss1$name]
tss1$CS_L57 <- cumsum(sort(tssMm[,"L57"] / sum(tssMm[,"L57"])))[tss1$name]

top50pct <- list(
  DS = tss1$name[tss1$CS_DS>=.5],
  L6 = tss1$name[tss1$CS_L6>=.5],
  L26 = tss1$name[tss1$CS_L26>=.5],
  L57 = tss1$name[tss1$CS_L57>=.5]
)
print(lengths(top50pct))

tss1_ucsc <- tss1
seqlevels(tss1_ucsc) <- paste0("Chr", 1:5)

mn_ds <- get_unstranded_signal(MNASE_BW,
  resize(tss1_ucsc[tss1_ucsc$name %in% top50pct$DS], 1000, "center"),row.norm=F)
mn_l6 <- get_unstranded_signal(MNASE_BW,
  resize(tss1_ucsc[tss1_ucsc$name %in% top50pct$L6], 1000, "center"),row.norm=F)
mn_l26 <- get_unstranded_signal(MNASE_BW,
  resize(tss1_ucsc[tss1_ucsc$name %in% top50pct$L26], 1000, "center"),row.norm=F)
mn_l57 <- get_unstranded_signal(MNASE_BW,
  resize(tss1_ucsc[tss1_ucsc$name %in% top50pct$L57], 1000, "center"),row.norm=F)

pdf("Figure3/SuppFig3f_mnase_top_timepoint_genes.pdf", width = 3.75, height = 3.25)
plot(mn_ds$Coord, mn_ds$Avg, type = "l", lwd = 3, bty = "n",
  xlab = "Distance from TSS (bp)", ylab = "Leaf MNase-seq (PlantDHS)",
  xaxt = "n", col = colorBlindBlack8[1], ylim = c(0, 100)
)
axis(1, c(1, 501, 1000), c("-500", 0, "500"))
lines(mn_l6$Coord, mn_l6$Avg, lwd = 3, col = colorBlindBlack8[4])
lines(mn_l26$Coord, mn_l26$Avg, lwd = 3, col = colorBlindBlack8[5])
lines(mn_l57$Coord, mn_l57$Avg, lwd = 3, col = colorBlindBlack8[6])
dev.off()
