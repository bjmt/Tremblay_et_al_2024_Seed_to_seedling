# Benjamin J.M. Tremblay
#
# Figure 1 and Supplementary Figures 1, 11, 12
#
# Package requirements: matrixStats, cowplot, ggplot2, rtracklayer,
# EnrichedHeatmap, GenomicFeatures, ComplexHeatmap, circlize, pals, readr,
# BSgenome.Athaliana.TAIR.TAIR9, universalmotif, viridisLite, dplyr, genomation

library(matrixStats)
library(cowplot)
library(ggplot2)
library(rtracklayer)
library(EnrichedHeatmap)
library(GenomicFeatures)
library(ComplexHeatmap)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

source("Functions.R")

dir.create("Figure1", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_CLUSTERS <- "csRNAseq_clustering/Clusters_anno.tsv"
ACR_CLUSTERS <- "ATACseq_clustering3/Clusters_anno.tsv"
RNA_TPM <- "RNAseq.quantifications.tpm.tsv"
ATAC_CPM <- "ATAC_peaks_cpm_fixed.txt"
TSS_BED <- "tss.size1.bed"
ARAPORT11_GFF <- "Araport11_novel.gff3"
READ_COUNTS_DIR <- "read-counts"
CSRNA_TAGDIRS <- "csrnaseq_tagDirs"
SRNA_TAGDIRS <- "srnaseq_tagDirs"

# Read density
CSRNA_57H_PLUS_BW <- "L57.merged.csRNA.plus.bw"
CSRNA_57H_MINUS_BW <- "L57.merged.csRNA.minus.bw"
RNA_57H_PLUS_BW <- "L57.merged.RNA.plus.bw"
RNA_57H_MINUS_BW <- "L57.merged.RNA.minus.bw"
SRNA_57H_PLUS_BG <- "srnaseq_bedGraphs/L57.merged.plus.bedGraph.gz"
SRNA_57H_MINUS_BG <- "srnaseq_bedGraphs/L57.merged.minus.bedGraph.gz"
SRNA_57H_TSS_PLUS_BG <- "srnaseq_bedGraphs/L57.merged.tss.plus.bedGraph.gz"
SRNA_57H_TSS_MINUS_BG <- "srnaseq_bedGraphs/L57.merged.tss.minus.bedGraph.gz"
GROCAP_PLUS_BG <- "external_bedGraphs/GROcap.5p.RPM.Signal.Unique.str1.out.bedGraph.gz"
GROCAP_MINUS_BG <- "external_bedGraphs/GROcap.5p.RPM.Signal.Unique.str2.out.bedGraph.gz"
CAGE_BG <- "external_bedGraphs/CAGE_14D_merged.bedGraph.gz"
ATAC_DS_BW <- "atacseq_bw/DS.bw"
ATAC_L6_BW <- "atacseq_bw/L6.bw"
ATAC_L26_BW <- "atacseq_bw/L26.bw"
ATAC_57H_BW <- "atacseq_bw/L57.bw"

# Figure 1b --------------------------------------------------------------------

tssA <- readr::read_tsv(TSS_TABLE)
atac <- as.matrix(read.table(ATAC_CPM, check.names = FALSE))

tssM <- as.matrix(tssA[, 10:26])
totM <- as.matrix(tssA[, 28:45])

# 2000 most variable TSSs
tssM2 <- tssM[, -17]
tssM2 <- tssM2[order(rowVars(tssM2), decreasing = TRUE)[1:2000], ]
tssM_pc <- prcomp(t(log2(1+tssM2)), scale = F)

# 2000 most variable transcripts
totM2 <- totM[, -(17:18)]
totM2 <- totM2[order(rowVars(totM2), decreasing = TRUE)[1:2000], ]
totM_pc <- prcomp(t(log2(1+totM2)), scale = F)

atac2 <- atac[order(rowVars(atac), decreasing = TRUE)[1:2000], ]
atac_pc <- prcomp(t(log2(1+atac)), scale = F)

tssM_df <- data.frame(row.names = NULL,
  Sample = rownames(tssM_pc$x),
  Group = gsub("[1257]$", "", rownames(tssM_pc$x)),
  PC1 = tssM_pc$x[, 1],
  PC2 = tssM_pc$x[, 2]
)
totM_df <- data.frame(row.names = NULL,
  Sample = gsub("RNAseq_", "", rownames(totM_pc$x)),
  Group = gsub("[1257]$", "", gsub("RNAseq_", "", rownames(totM_pc$x))),
  PC1 = totM_pc$x[, 1],
  PC2 = totM_pc$x[, 2]
)
atac_df <- data.frame(row.names = NULL,
  Sample = rownames(atac_pc$x),
  Group = gsub("[12]$", "", rownames(atac_pc$x)),
  PC1 = atac_pc$x[, 1],
  PC2 = atac_pc$x[, 2]
)

# Variance explained, for the axis labels
summary(tssM_pc)
summary(totM_pc)
summary(atac_pc)

# Nature Methods 8:441
colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p1 <- ggplot(tssM_df, aes(PC1, PC2, fill = Group)) +
  geom_point(pch = 21, size = 1) +
  xlab("PC1 (74.4%)") +
  ylab("PC2 (12.4%)") +
  scale_fill_manual(values = colorBlindBlack8,
    name = element_blank(),
    labels = c("DS", "S24", "S72", "L6", "L26", "L57", expression(italic("hen2-4")), expression(italic("rrp4-2"))),
    limits = c("DS", "I24", "I72", "6H", "26H", "57H", "H","R")) +
  theme_bw() +
  theme(
    panel.grid = element_line(linetype = 2, colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = "bottom",
    legend.text.align = 0,
    legend.spacing.x = unit(0, "pt"),
    legend.text = element_text(colour = "black", size = unit(5, "pt"), margin = margin(l = 0)),
    panel.border = element_rect(colour = "black", size = 1),
    axis.text.x = element_text(size = unit(5, "pt"), colour = "black"),
    axis.text.y = element_text(size = unit(5, "pt"), colour = "black"),
    axis.title.x = element_text(size = unit(6, "pt"), colour = "black"),
    axis.title.y = element_text(size = unit(6, "pt"), colour = "black")
  ) + guides(fill = guide_legend(nrow = 1))

p2 <- ggplot(totM_df, aes(PC1, PC2, fill = Group)) +
  geom_point(pch = 21, size = 1) +
  xlab("PC1 (75.3%)") +
  ylab("PC2 (16.4%)") +
  scale_fill_manual(values = colorBlindBlack8,
    limits = c("DS", "I24", "I72", "6H", "26H", "57H", "H","R")) +
  theme_bw() +
  theme(
    panel.grid = element_line(linetype = 2, colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    axis.text.x = element_text(size = unit(5, "pt"), colour = "black"),
    axis.text.y = element_text(size = unit(5, "pt"), colour = "black"),
    axis.title.x = element_text(size = unit(6, "pt"), colour = "black"),
    axis.title.y = element_text(size = unit(6, "pt"), colour = "black")
  )

p3 <- ggplot(atac_df, aes(PC1, PC2, fill = Group)) +
  geom_point(pch = 21, size = 1) +
  xlab("PC1 (38.0%)") +
  ylab("PC2 (30.5%)") +
  scale_fill_manual(values = colorBlindBlack8[c(1, 4:6)],
    limits = c("DS", "6H", "26H", "57H")) +
  theme_bw() +
  theme(
    panel.grid = element_line(linetype = 2, colour = "grey50", size = 0.2),
    panel.grid.minor = element_blank(),
    aspect.ratio = 1,
    legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    axis.text.x = element_text(size = unit(5, "pt"), colour = "black"),
    axis.text.y = element_text(size = unit(5, "pt"), colour = "black"),
    axis.title.x = element_text(size = unit(6, "pt"), colour = "black"),
    axis.title.y = element_text(size = unit(6, "pt"), colour = "black")
  )

lgnd <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

save_plot("Figure1/Fig1b_pca.pdf",
  plot_grid(rel_heights = c(1, 0.1),
    plot_grid(p1, p2, p3, nrow = 1, align = "hv", axis = "tblr"), lgnd, ncol = 1),
  base_height = grid::convertUnit(unit(40, "mm"), "in"),
  base_width = grid::convertUnit(unit(90, "mm"), "in"))

# Figure 1d --------------------------------------------------------------------

csrna <- read_signal(CSRNA_57H_PLUS_BW, CSRNA_57H_MINUS_BW)
rna <- read_signal(RNA_57H_PLUS_BW, RNA_57H_MINUS_BW)
srna <- read_signal(SRNA_57H_PLUS_BG, SRNA_57H_MINUS_BG)
grocap <- read_signal(GROCAP_PLUS_BG, GROCAP_MINUS_BG)
cage <- import(CAGE_BG)
strand(cage[[1]]) <- "+"
strand(cage[[2]]) <- "-"
cage <- sort(c(cage[[1]], cage[[2]]))

atac <- import(ATAC_57H_BW)

tss1 <- import(TSS_BED)

# Keep the most expressed TSS of each mRNA
tssMe <- tssM
rownames(tssMe) <- tssA$TSS
tssMe <- tssMe[rownames(tssMe) %in%
  tssA$TSS[!is.na(tssA$NearestFeatureType) & tssA$NearestFeatureType == "mRNA"], ]

tssMe <- tssMe[order(rowMedians(tssMe, na.rm = TRUE), decreasing = TRUE), ]
tssMeG <- structure(tssA$NearestFeature, names = tssA$TSS)[rownames(tssMe)]
tssMe <- tssMe[!duplicated(unname(tssMeG)), ]
tssMeG <- tssMeG[rownames(tssMe)]
tssMeR <- structure(names(tssMeG), names = unname(tssMeG))

txdb <- makeTxDbFromGFF(ARAPORT11_GFF)
tx <- transcripts(txdb)
tx <- tx[tx$tx_name %in% unname(tssMeG)]

# From the TSS to the annotated transcript end
tssMeP1 <- structure(start(tss1), names = tss1$name)[rownames(tssMe)]
tssMeP1 <- tssMeP1[tssMeR[tx$tx_name]]
txEnd <- end(tx)
txEnd[as.character(strand(tx)) == "-"] <- start(tx)[as.character(strand(tx)) == "-"]

tx2 <- GRanges(seqnames(tx), IRanges(pmin(tssMeP1, txEnd), pmax(tssMeP1, txEnd)),
    strand = strand(tx), Transcript = tx$tx_name, TSS = tssMeR[tx$tx_name])

csrna_tx <- normalizeToMatrix(csrna, tx2, extend = 2000, value_column = "score",
  w = 20, keep = c(0.05, .95))
grocap_tx <- normalizeToMatrix(grocap, tx2, extend = 2000, value_column = "score",
  w = 20, keep = c(0.05, .95))
cage_tx <- normalizeToMatrix(cage, tx2, extend = 2000, value_column = "score",
  w = 20, keep = c(0.05, .95))
rna_tx <- normalizeToMatrix(rna, tx2, extend = 2000, value_column = "score",
  w = 20, keep = c(0.05, .95))
srna_tx <- normalizeToMatrix(srna, tx2, extend = 2000, value_column = "score",
  w = 20, keep = c(0.05, .95))
atac_tx <- normalizeToMatrix(atac, tx2, extend = 2000, value_column = "score",
  w = 20, keep = c(0.01, .99))

watCols <- pals::watlington()

pdf("Figure1/Fig1d_heatmaps.pdf",
  height = grid::convertUnit(unit(80, "mm"), "in"),
  width = grid::convertUnit(unit(90, "mm"), "in"))
draw(heatmap_legend_side = "bottom",
EnrichedHeatmap(csrna_tx, column_title = "csRNA-seq",
  column_title_gp = gpar(fontsize = unit(6, "pt")),
  axis_name = c("", "TSS", "TTS", ""),
  axis_name_gp = gpar(fontsize = unit(5, "pt")),
  top_annotation = HeatmapAnnotation(
    border = FALSE,
    enriched = anno_enriched(gp = gpar(col = "black"), yaxis = FALSE,
      height = unit(10, "mm"))
  ),
  heatmap_legend_param = list(title = "", labels = c("Min", "Max"),
    legend_width = unit(10, "mm"), legend_height = unit(0.5, "mm"),
    at = range(csrna_tx),
    labels_gp = gpar(fontsize = unit(5, "pt")), direction = "horizontal"),
  col = c("white", watCols[3]),
  raster_quality = 3,
  pos_line = FALSE
) +
EnrichedHeatmap(srna_tx, column_title = "sRNA-seq",
  column_title_gp = gpar(fontsize = unit(6, "pt")),
  axis_name = c("", "TSS", "TTS", ""),
  axis_name_gp = gpar(fontsize = unit(5, "pt")),
  top_annotation = HeatmapAnnotation(
    border = FALSE,
    enriched = anno_enriched(gp = gpar(col = "black"), yaxis = FALSE,
      height = unit(10, "mm"))
  ),
  heatmap_legend_param = list(title = "", labels = c("Min", "Max"),
    legend_width = unit(10, "mm"), legend_height = unit(0.5, "mm"),
    at = range(srna_tx),
    labels_gp = gpar(fontsize = unit(5, "pt")), direction = "horizontal"),
  col = c("white", watCols[4]),
  raster_quality = 3,
  pos_line = FALSE
) +
EnrichedHeatmap(rna_tx, column_title = "RNA-seq",
  column_title_gp = gpar(fontsize = unit(6, "pt")),
  axis_name = c("", "TSS", "TTS", ""),
  axis_name_gp = gpar(fontsize = unit(5, "pt")),
  top_annotation = HeatmapAnnotation(
    border = FALSE,
    enriched = anno_enriched(gp = gpar(col = "black"), yaxis = FALSE,
      height = unit(10, "mm"))
  ),
  heatmap_legend_param = list(title = "", labels = c("Min", "Max"),
    legend_width = unit(10, "mm"), legend_height = unit(0.5, "mm"),
    at = range(rna_tx),
    labels_gp = gpar(fontsize = unit(5, "pt")), direction = "horizontal"),
  col = c("white", watCols[5]),
  raster_quality = 3,
  pos_line = FALSE
) +
EnrichedHeatmap(atac_tx, column_title = "ATAC-seq",
  column_title_gp = gpar(fontsize = unit(6, "pt")),
  axis_name = c("", "TSS", "TTS", ""),
  axis_name_gp = gpar(fontsize = unit(5, "pt")),
  top_annotation = HeatmapAnnotation(
    border = FALSE,
    enriched = anno_enriched(gp = gpar(col = "black"), yaxis = FALSE,
      height = unit(10, "mm"))
  ),
  heatmap_legend_param = list(title = "", labels = c("Min", "Max"),
    legend_width = unit(10, "mm"), legend_height = unit(0.5, "mm"),
    at = range(atac_tx),
    labels_gp = gpar(fontsize = unit(5, "pt")), direction = "horizontal"),
  col = circlize::colorRamp2(c(0, 3), c("white", watCols[6])),
  raster_quality = 3,
  pos_line = FALSE
) +
EnrichedHeatmap(grocap_tx, column_title = "GRO-cap (Hetzel et al.)",
  column_title_gp = gpar(fontsize = unit(6, "pt")),
  axis_name = c("", "TSS", "TTS", ""),
  axis_name_gp = gpar(fontsize = unit(5, "pt")),
  top_annotation = HeatmapAnnotation(
    border = FALSE,
    enriched = anno_enriched(gp = gpar(col = "black"), yaxis = FALSE,
      height = unit(10, "mm"))
  ),
  heatmap_legend_param = list(title = "", labels = c("Min", "Max"),
    legend_width = unit(10, "mm"), legend_height = unit(0.5, "mm"),
    at = range(grocap_tx),
    labels_gp = gpar(fontsize = unit(5, "pt")), direction = "horizontal"),
  col = c("white", watCols[7]),
  raster_quality = 3,
  pos_line = FALSE
) +
EnrichedHeatmap(cage_tx, column_title = "CAGE (Thieffry et al.)",
  column_title_gp = gpar(fontsize = unit(6, "pt")),
  axis_name = c("", "TSS", "TTS", ""),
  axis_name_gp = gpar(fontsize = unit(5, "pt")),
  top_annotation = HeatmapAnnotation(
    border = FALSE,
    enriched = anno_enriched(gp = gpar(col = "black"), yaxis = FALSE,
      height = unit(10, "mm"))
  ),
  heatmap_legend_param = list(title = "", labels = c("Min", "Max"),
    legend_width = unit(10, "mm"), legend_height = unit(0.5, "mm"),
    at = range(cage_tx),
    labels_gp = gpar(fontsize = unit(5, "pt")), direction = "horizontal"),
  col = c("white", watCols[12]),
  raster_quality = 3,
  pos_line = FALSE
)
)
dev.off()

# Supplementary Figure 1a ------------------------------------------------------

tssA <- readr::read_tsv(TSS_CLUSTERS)
rna <- readr::read_tsv(RNA_TPM)
atac <- readr::read_tsv(ACR_CLUSTERS)

tssM <- as.matrix(tssA[, 11:26])
rownames(tssM) <- tssA$TSS

tssR <- as.matrix(tssA[, 29:44])
rownames(tssR) <- tssA$TSS

rnaM <- as.matrix(rna[, -1])
rownames(rnaM) <- rna$Transcript
rnaM <- rnaM[apply(rnaM, 1, function(x) sum(x>=2)>=2), ]

atacM <- as.matrix(atac[, 13:20])
rownames(atacM) <- atac$PeakID

tssMm <- as.matrix(data.frame(row.names = rownames(tssM),
  DS = rowMeans(tssM[, 1:2]),
  S24 = rowMeans(tssM[, 3:4]),
  S72 = rowMeans(tssM[, 5:6]),
  L6 = rowMeans(tssM[, 7:8]),
  L26 = rowMeans(tssM[, 9:10]),
  L57 = rowMeans(tssM[, 11:12]),
  hen2 = rowMeans(tssM[, 13:14]),
  rrp4 = rowMeans(tssM[, 15:16])
))

snames <- c("DS.1", "DS.2", "S24.1", "S24.2", "S72.1", "S72.2",
  "L6.1", "L6.2", "L26.1", "L26.2", "L57.1", "L57.2", "hen2-4.1",
  "hen2-4.2", "rrp4-2.1", "rrp4-2.2")
snames2 <- c("DS.1", "DS.2", "L6.1", "L6.2", "L26.1", "L26.2", "L57.1", "L57.2")

pdf("Figure1/SuppFig1a_corr_matrix_csrna.pdf", width = 4, height = 3.5)
draw(Heatmap(cor(log2(1+tssM)),
  name = "PCC",
  row_labels = snames, column_labels = snames,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = viridisLite::magma(100)
))
dev.off()

# Supplementary Figure 1b ------------------------------------------------------

pdf("Figure1/SuppFig1b_corr_matrix_rna.pdf", width = 4, height = 3.5)
draw(Heatmap(cor(log2(1+rnaM)),
  name = "PCC",
  row_labels = snames, column_labels = snames,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = viridisLite::magma(100)
))
dev.off()

# Supplementary Figure 1c ------------------------------------------------------

pdf("Figure1/SuppFig1c_corr_matrix_atac.pdf", width = 2.5, height = 2)
draw(Heatmap(cor(log2(1+atacM)),
  name = "PCC",
  row_labels = snames2, column_labels = snames2,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = viridisLite::magma(100)
))
dev.off()

# Supplementary Figure 1f ------------------------------------------------------

tss1 <- import(TSS_BED)

v4 <- viridisLite::viridis(4)

atac_ds_0_10 <- get_unstranded_signal(ATAC_DS_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"DS"]<=10]],4000,"center"))
atac_ds_10_100 <- get_unstranded_signal(ATAC_DS_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"DS"]>10&tssMm[,"DS"]<=100]],4000,"center"))
atac_ds_100_1000 <- get_unstranded_signal(ATAC_DS_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"DS"]>100&tssMm[,"DS"]<=1000]],4000,"center"))
atac_ds_1000 <- get_unstranded_signal(ATAC_DS_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"DS"]>1000]],4000,"center"))

pdf("Figure1/SuppFig1f_atac_tss_pileup_DS.pdf", width = 3.75, height = 3.25)
plot(atac_ds_0_10$Coord, atac_ds_0_10$Avg, type = "l", lwd = 3, bty = "n",
  xlab = "Distance to TSS", col = v4[1], ylab = "ATAC-seq RPM",
  xaxt = "n", ylim = c(0.7, 7))
axis(1, c(1, 2000, 4000), c("-2 kb", "0", "+2 kb"))
lines(atac_ds_10_100$Coord, atac_ds_10_100$Avg, col = v4[2], lwd = 3)
lines(atac_ds_100_1000$Coord, atac_ds_100_1000$Avg, col = v4[3], lwd = 3)
lines(atac_ds_1000$Coord, atac_ds_1000$Avg, col = v4[4], lwd = 3)
legend("topleft", col = v4, lwd = 3, legend = c("", "", "", ""),bty="n")
dev.off()

# Supplementary Figure 1g ------------------------------------------------------

atac_l6_0_10 <- get_unstranded_signal(ATAC_L6_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L6"]<=10]],4000,"center"))
atac_l6_10_100 <- get_unstranded_signal(ATAC_L6_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L6"]>10&tssMm[,"DS"]<=100]],4000,"center"))
atac_l6_100_1000 <- get_unstranded_signal(ATAC_L6_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L6"]>100&tssMm[,"DS"]<=1000]],4000,"center"))
atac_l6_1000 <- get_unstranded_signal(ATAC_L6_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L6"]>1000]],4000,"center"))

pdf("Figure1/SuppFig1g_atac_tss_pileup_L6.pdf", width = 3.75, height = 3.25)
plot(atac_l6_0_10$Coord, atac_l6_0_10$Avg, type = "l", lwd = 3, bty = "n",
  xlab = "Distance to TSS", col = v4[1], ylab = "ATAC-seq RPM",
  xaxt = "n", ylim = c(0.7, 7))
axis(1, c(1, 2000, 4000), c("-2 kb", "0", "+2 kb"))
lines(atac_l6_10_100$Coord, atac_l6_10_100$Avg, col = v4[2], lwd = 3)
lines(atac_l6_100_1000$Coord, atac_l6_100_1000$Avg, col = v4[3], lwd = 3)
lines(atac_l6_1000$Coord, atac_l6_1000$Avg, col = v4[4], lwd = 3)
legend("topleft", col = v4, lwd = 3, legend = c("", "", "", ""),bty="n")
dev.off()

# Supplementary Figure 1h ------------------------------------------------------

atac_l26_0_10 <- get_unstranded_signal(ATAC_L26_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L26"]<=10]],4000,"center"))
atac_l26_10_100 <- get_unstranded_signal(ATAC_L26_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L26"]>10&tssMm[,"DS"]<=100]],4000,"center"))
atac_l26_100_1000 <- get_unstranded_signal(ATAC_L26_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L26"]>100&tssMm[,"DS"]<=1000]],4000,"center"))
atac_l26_1000 <- get_unstranded_signal(ATAC_L26_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L26"]>1000]],4000,"center"))

pdf("Figure1/SuppFig1h_atac_tss_pileup_L26.pdf", width = 3.75, height = 3.25)
plot(atac_l26_0_10$Coord, atac_l26_0_10$Avg, type = "l", lwd = 3, bty = "n",
  xlab = "Distance to TSS", col = v4[1], ylab = "ATAC-seq RPM",
  xaxt = "n", ylim = c(0.7, 7))
axis(1, c(1, 2000, 4000), c("-2 kb", "0", "+2 kb"))
lines(atac_l26_10_100$Coord, atac_l26_10_100$Avg, col = v4[2], lwd = 3)
lines(atac_l26_100_1000$Coord, atac_l26_100_1000$Avg, col = v4[3], lwd = 3)
lines(atac_l26_1000$Coord, atac_l26_1000$Avg, col = v4[4], lwd = 3)
legend("topleft", col = v4, lwd = 3, legend = c("", "", "", ""),bty="n")
dev.off()

# Supplementary Figure 1i ------------------------------------------------------

atac_l57_0_10 <- get_unstranded_signal(ATAC_57H_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L57"]<=10]],4000,"center"))
atac_l57_10_100 <- get_unstranded_signal(ATAC_57H_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L57"]>10&tssMm[,"DS"]<=100]],4000,"center"))
atac_l57_100_1000 <- get_unstranded_signal(ATAC_57H_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L57"]>100&tssMm[,"DS"]<=1000]],4000,"center"))
atac_l57_1000 <- get_unstranded_signal(ATAC_57H_BW, row.norm=F,
  resize(tss1[tss1$name %in% rownames(tssMm)[tssMm[,"L57"]>1000]],4000,"center"))

pdf("Figure1/SuppFig1i_atac_tss_pileup_L57.pdf", width = 3.75, height = 3.25)
plot(atac_l57_0_10$Coord, atac_l57_0_10$Avg, type = "l", lwd = 3, bty = "n",
  xlab = "Distance to TSS", col = v4[1], ylab = "ATAC-seq RPM",
  xaxt = "n", ylim = c(0.7, 7))
axis(1, c(1, 2000, 4000), c("-2 kb", "0", "+2 kb"))
lines(atac_l57_10_100$Coord, atac_l57_10_100$Avg, col = v4[2], lwd = 3)
lines(atac_l57_100_1000$Coord, atac_l57_100_1000$Avg, col = v4[3], lwd = 3)
lines(atac_l57_1000$Coord, atac_l57_1000$Avg, col = v4[4], lwd = 3)
legend("topleft", col = v4, lwd = 3, legend = c("", "", "", ""),bty="n")
dev.off()

# Supplementary Figure 1j, k ---------------------------------------------------

tssMr <- tssM[!is.na(tssR[,1]), ]
tssRr <- tssR[!is.na(tssR[,1]), ]

gnames <- gsub("[.]\\d", "", snames)

a1 <- apply(tssMr, 2, function(x) sum(x>=1))
a3 <- apply(tssMr, 2, function(x) sum(x>=5))
a5 <- apply(tssMr, 2, function(x) sum(x>=10))
a7 <- apply(tssMr, 2, function(x) sum(x>=20))
a2 <- apply(tssRr, 2, function(x) sum(x>=1))
a4 <- apply(tssRr, 2, function(x) sum(x>=5))
a6 <- apply(tssRr, 2, function(x) sum(x>=10))
a8 <- apply(tssRr, 2, function(x) sum(x>=20))
am <- data.frame(
  samples = factor(c(gnames, gnames, gnames, gnames), levels = unique(gnames)),
  counts = c(a1, a3, a5, a7),
  threshold = c(rep(1, 16), rep(5, 16), rep(10, 16), rep(20, 16))
)
ar <- data.frame(
  samples = factor(c(gnames, gnames, gnames, gnames), levels = unique(gnames)),
  counts = c(a2, a4, a6, a8),
  threshold = c(rep(1, 16), rep(5, 16), rep(10, 16), rep(20, 16))
)

p1 <- ggplot(am, aes(samples, counts, colour = factor(threshold), group = threshold)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 25000)) +
  geom_line(
    data = dplyr::summarize(dplyr::group_by(am, samples, threshold), counts = mean(counts)),
    aes(samples, counts, group = threshold)
  ) +
  scale_colour_manual(values = colorBlindBlack8[1:4], name = "CPM threshold") +
  xlab(element_blank()) +
  ylab("Number of detected TSSs") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    legend.title = element_text(colour = "black", size = unit(5, "pt")),
    panel.grid = element_blank()
  )

p2 <- ggplot(ar, aes(samples, counts, colour = factor(threshold), group = threshold)) +
  geom_point() +
  scale_y_continuous(limits = c(0, 25000)) +
  geom_line(
    data = dplyr::summarize(dplyr::group_by(ar, samples, threshold), counts = mean(counts)),
    aes(samples, counts, group = threshold)
  ) +
  scale_colour_manual(values = colorBlindBlack8[1:4], name = "TPM threshold") +
  xlab(element_blank()) +
  ylab("Number of detected Transcripts") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    legend.title = element_text(colour = "black", size = unit(5, "pt")),
    panel.grid = element_blank()
  )

save_plot("Figure1/SuppFig1jk_num_detected_features.pdf",
  plot_grid(p1, p2, nrow = 1),
  base_width = 8.5 * 2 / 3, base_height = 1.5)

# Supplementary Figure 1l ------------------------------------------------------

tss1$CS_DS <- cumsum(sort(tssMm[,"DS"] / sum(tssMm[,"DS"])))[tss1$name]
tss1$CS_S24 <- cumsum(sort(tssMm[,"S24"] / sum(tssMm[,"S24"])))[tss1$name]
tss1$CS_S72 <- cumsum(sort(tssMm[,"S72"] / sum(tssMm[,"S72"])))[tss1$name]
tss1$CS_L6 <- cumsum(sort(tssMm[,"L6"] / sum(tssMm[,"L6"])))[tss1$name]
tss1$CS_L26 <- cumsum(sort(tssMm[,"L26"] / sum(tssMm[,"L26"])))[tss1$name]
tss1$CS_L57 <- cumsum(sort(tssMm[,"L57"] / sum(tssMm[,"L57"])))[tss1$name]
tss1$CS_hen2 <- cumsum(sort(tssMm[,"hen2"] / sum(tssMm[,"hen2"])))[tss1$name]
tss1$CS_rrp4 <- cumsum(sort(tssMm[,"rrp4"] / sum(tssMm[,"rrp4"])))[tss1$name]

pdf("Figure1/SuppFig1l_tss_cumfrac_persample.pdf", width = 4.5, height = 4)
plot(1:length(tss1),sort(tss1$CS_DS),type="l",bty="n",col=colorBlindBlack8[1],lwd=2,
  xlab = "TSS rank", ylab = "Cumulative csRNA-seq fraction"
)
lines(1:length(tss1),sort(tss1$CS_S24),col=colorBlindBlack8[2],lwd=2)
lines(1:length(tss1),sort(tss1$CS_S72),col=colorBlindBlack8[3],lwd=2)
lines(1:length(tss1),sort(tss1$CS_L6),col=colorBlindBlack8[4],lwd=2)
lines(1:length(tss1),sort(tss1$CS_L26),col=colorBlindBlack8[5],lwd=2)
lines(1:length(tss1),sort(tss1$CS_L57),col=colorBlindBlack8[6],lwd=2)
lines(1:length(tss1),sort(tss1$CS_hen2),col=colorBlindBlack8[7],lwd=2)
lines(1:length(tss1),sort(tss1$CS_rrp4),col=colorBlindBlack8[8],lwd=2)
legend("topleft", lty = 1, col = colorBlindBlack8, bty = "n", lwd =2,
  legend = colnames(tssMm))
dev.off()

# Supplementary Figure 1m ------------------------------------------------------

ggplot(mapping = aes(Timepoint, Count),
data.frame(
  Timepoint = factor(rep(colnames(tssMm), each = 3), levels = colnames(tssMm)),
  Threshold = rep(c("25%", "50%", "75%"), times = 8),
  Count = c(
    sum(tss1$CS_DS>=.75),
    sum(tss1$CS_DS>=.5),
    sum(tss1$CS_DS>=.25),
    sum(tss1$CS_S24>=.75),
    sum(tss1$CS_S24>=.5),
    sum(tss1$CS_S24>=.25),
    sum(tss1$CS_S72>=.75),
    sum(tss1$CS_S72>=.5),
    sum(tss1$CS_S72>=.25),
    sum(tss1$CS_L6>=.75),
    sum(tss1$CS_L6>=.5),
    sum(tss1$CS_L6>=.25),
    sum(tss1$CS_L26>=.75),
    sum(tss1$CS_L26>=.5),
    sum(tss1$CS_L26>=.25),
    sum(tss1$CS_L57>=.75),
    sum(tss1$CS_L57>=.5),
    sum(tss1$CS_L57>=.25),
    sum(tss1$CS_hen2>=.75),
    sum(tss1$CS_hen2>=.5),
    sum(tss1$CS_hen2>=.25),
    sum(tss1$CS_rrp4>=.75),
    sum(tss1$CS_rrp4>=.5),
    sum(tss1$CS_rrp4>=.25)
  )
)) +
  geom_bar(stat = "identity", width = .8, fill = "black") +
  theme_bw() +
  scale_y_continuous(expand = c(0, 0), name = "Number of top expressed TSSs per cumulative fraction") +
  xlab(element_blank()) +
  theme(
    legend.text = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  facet_wrap(~Threshold, scales = "free_x") +
  coord_flip()
ggsave("Figure1/SuppFig1m_top_expressed_tss_counts.pdf", width = 3.65, height = 1.62)

# Supplementary Figure 1n ------------------------------------------------------

tss1fa2kb <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, promoters(tss1, 500, 500))
names(tss1fa2kb) <- tss1$name
tss1fa2kbM <- as.matrix(t(as.data.frame(strsplit(as.character(tss1fa2kb), "", fixed=T))))

get_freqs <- function(x) {
  list(
    A = apply(x, 2, function(x) sum(x=="A")) / nrow(x),
    C = apply(x, 2, function(x) sum(x=="C")) / nrow(x),
    G = apply(x, 2, function(x) sum(x=="G")) / nrow(x),
    T = apply(x, 2, function(x) sum(x=="T")) / nrow(x)
  )
}

f_pc <- get_freqs(tss1fa2kbM[tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType=="protein_coding"], ])
f_nc <- get_freqs(tss1fa2kbM[tssA$TSS[is.na(tssA$NearestGeneType) | tssA$NearestGeneType!="protein_coding"], ])

f_pc100 <- get_freqs(tss1fa2kbM[tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType=="protein_coding"], 400:599])
f_nc100 <- get_freqs(tss1fa2kbM[tssA$TSS[is.na(tssA$NearestGeneType) | tssA$NearestGeneType!="protein_coding"], 400:599])

pdf("Figure1/SuppFig1n_tss_pc_base_freq.pdf", width = 3.8, height = 3.5)
plot(c((-499):0, 1:500), f_pc$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-500, 0, 500))
axis(2, c(0, .5, 1))
lines(c((-499):0, 1:500), f_pc$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-499):0, 1:500), f_pc$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-499):0, 1:500), f_pc$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 1o ------------------------------------------------------

pdf("Figure1/SuppFig1o_tss_nc_base_freq.pdf", width = 3.8, height = 3.5)
plot(c((-499):0, 1:500), f_nc$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-500, 0, 500))
axis(2, c(0, .5, 1))
lines(c((-499):0, 1:500), f_nc$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-499):0, 1:500), f_nc$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-499):0, 1:500), f_nc$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 1p ------------------------------------------------------

pdf("Figure1/SuppFig1p_tss_pc_base_freq_zoom.pdf", width = 3.8, height = 3.5)
plot(c((-99):0, 1:100), f_pc100$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-100, 0, 100))
axis(2, c(0, .5, 1))
lines(c((-99):0, 1:100), f_pc100$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-99):0, 1:100), f_pc100$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-99):0, 1:100), f_pc100$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 1q ------------------------------------------------------

pdf("Figure1/SuppFig1q_tss_nc_base_freq_zoom.pdf", width = 3.8, height = 3.5)
plot(c((-99):0, 1:100), f_nc100$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-100, 0, 100))
axis(2, c(0, .5, 1))
lines(c((-99):0, 1:100), f_nc100$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-99):0, 1:100), f_nc100$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-99):0, 1:100), f_nc100$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 1r ------------------------------------------------------

L57 <- read_signal(CSRNA_57H_PLUS_BW, CSRNA_57H_MINUS_BW)
L57s <- read_signal(SRNA_57H_TSS_PLUS_BG, SRNA_57H_TSS_MINUS_BG)
GROcap <- read_signal(GROCAP_PLUS_BG, GROCAP_MINUS_BG)

CAGE <- import(CAGE_BG)
strand(CAGE[[1]]) <- "+"
strand(CAGE[[2]]) <- "-"
CAGE <- sort(c(CAGE[[1]], CAGE[[2]]))

tss1pc <- tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType=="protein_coding"]]

tss1pc_gb <- resize(resize(tss1pc, 1200, "start"), 1000, "end")

L57_tss <- rowSums(get_stranded_signal(L57, resize(tss1pc, 100, "center"), return.matrix=T,row.norm=F)$Sense)
L57s_tss <- rowSums(get_stranded_signal(L57s, resize(tss1pc, 100, "center"), return.matrix=T,row.norm=F)$Sense)
GROcap_tss <- rowSums(get_stranded_signal(GROcap, resize(tss1pc, 100, "center"), return.matrix=T,row.norm=F)$Sense)
CAGE_tss <- rowSums(get_stranded_signal(CAGE, resize(tss1pc, 100, "center"), return.matrix=T,row.norm=F)$Sense)

L57_gb <- rowSums(get_stranded_signal(L57, tss1pc_gb, return.matrix=T,row.norm=F)$Sense)
L57s_gb <- rowSums(get_stranded_signal(L57s, tss1pc_gb, return.matrix=T,row.norm=F)$Sense)
GROcap_gb <- rowSums(get_stranded_signal(GROcap, tss1pc_gb, return.matrix=T,row.norm=F)$Sense)
CAGE_gb <- rowSums(get_stranded_signal(CAGE, tss1pc_gb, return.matrix=T,row.norm=F)$Sense)

pdf("Figure1/SuppFig1r_sn_ratio_cage_grocap.pdf", width = 3, height = 3.5)
boxplot(list(
  csRNA = log2((1+L57_tss) / (1+L57_gb)),
  sRNA = log2((1+L57s_tss) / (1+L57s_gb)),
  GROcap = log2((1+GROcap_tss) / (1+GROcap_gb)),
  CAGE = log2((1+CAGE_tss) / (1+CAGE_gb[names(CAGE_tss)]))
), outline=F)
dev.off()

wilcox.test(
  log2((1+L57_tss) / (1+L57_gb)),
  log2((1+L57s_tss) / (1+L57s_gb))
)$p.value

wilcox.test(
  log2((1+L57_tss) / (1+L57_gb)),
  log2((1+CAGE_tss) / (1+CAGE_gb[names(CAGE_tss)]))
)$p.value

# Supplementary Figure 11 ------------------------------------------------------

# Percent of reads per chromosome

s <- c("DS1", "DS2", "I241", "I242", "I725", "I727", "6H1", "6H2",
  "26H1", "26H2", "57H1", "57H2", "H1", "H2", "R1", "R2")

f_total <- structure(
  lapply(list.files(READ_COUNTS_DIR, pattern = "*total.tsv", full.names = TRUE),
    readr::read_tsv),
  names = gsub("_total[.]tsv$", "", list.files(READ_COUNTS_DIR, pattern = "*total.tsv")))

m_total <- matrix(0, nrow = length(s), ncol = 7,
  dimnames = list(s, c(1:5, "Mt", "Pt")))

for (i in seq_along(s)) {
  m_total[i, ] <- structure(f_total[[s[i]]]$count, names = f_total[[s[i]]]$chr)[colnames(m_total)]
}

p_total <- 100 * t(apply(m_total, 1, function(x) x / sum(x)))

write.table(p_total, "Figure1/total_reads.tab", quote = FALSE, sep = "\t")

pct_pt <- read.table("Figure1/total_reads.tab", check.names = FALSE)
pct_pt <- data.frame(row.names = NULL,
  Sample = rownames(pct_pt),
  Group = rep(c("DS", "S24", "S72", "L6", "L26", "L57", "H", "R"), each = 2),
  PctPt = pct_pt$Pt
)

p4 <- ggplot(pct_pt, aes(Group, PctPt, fill = Sample)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
  scale_x_discrete(limits = c("DS", "S24", "S72", "L6", "L26", "L57")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 50)) +
  scale_fill_manual(values = rep("black", 16)) +
  theme_bw() +
  ylab("% chloroplastic RNA") +
  xlab(element_blank()) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    axis.text = element_text(colour = "black", size = unit(6, "pt"))
  )

ggsave("Figure1/SuppFig11_pct_chloroplastic_rna.pdf", p4, width = 1.8, height = 1.15)

# Supplementary Figure 12 ------------------------------------------------------

samples <- c("DS1", "DS2", "S241", "S242", "S721", "S722", "L61", "L62",
  "L261", "L262", "L571", "L572", "hen21", "hen22", "rrp41", "rrp42")

LensCap <- lapply(samples,
  function(x) readr::read_tsv(paste0(CSRNA_TAGDIRS, "/", x, "-tagDir/tagLengthDistribution.txt")))
names(LensCap) <- samples

LensInp <- lapply(samples,
  function(x) readr::read_tsv(paste0(SRNA_TAGDIRS, "/", x, "-tagDir/tagLengthDistribution.txt")))
names(LensInp) <- samples

for (i in c("DS", "S24", "S72", "L6", "L26", "L57", "hen2", "rrp4")) {
  pdf(paste0("Figure1/SuppFig12_read_sizes_", i, ".pdf"), width = 3.5, height = 3)
  plot(LensInp[[paste0(i, 1)]][[1]], LensInp[[paste0(i, 1)]][[2]], type = "l", lwd = 3,
    col = "blue", ylim = c(0, 0.25), bty = "n", xlim = c(0, 80),
    xlab = "Read length", ylab = "Fraction of reads"
  )
  lines(LensInp[[paste0(i, 2)]][[1]], LensInp[[paste0(i, 2)]][[2]], lwd = 3, col = "blue")
  lines(LensCap[[paste0(i, 1)]][[1]], LensCap[[paste0(i, 1)]][[2]], lwd = 3, col = "black")
  lines(LensCap[[paste0(i, 2)]][[1]], LensCap[[paste0(i, 2)]][[2]], lwd = 3, col = "black")
  legend("topright", legend = c("Input", "Cap-enriched"), col = c("blue", "black"),
    lwd = 3, lty  =1, bty = "n")
  dev.off()
}
