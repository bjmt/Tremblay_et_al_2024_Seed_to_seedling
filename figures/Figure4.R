# Benjamin J.M. Tremblay
#
# Figure 4 and Supplementary Figures 4, 13
#
# Package requirements: rtracklayer, matrixStats, GenomicFeatures, universalmotif,
# EnrichedHeatmap, viridisLite, eulerr, BSgenome.Athaliana.TAIR.TAIR9, ggplot2,
# dplyr, readr, rstatix, zoo

library(rtracklayer)
library(matrixStats)
library(GenomicFeatures)
library(universalmotif)
library(EnrichedHeatmap)
library(eulerr)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(ggplot2)
library(dplyr)

dir.create("Figure4", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

A11_GFF <- "Araport11_novel.gff3"
TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_SIZE1_BED <- "tss.size1.bed"
ANTISENSE_PAIRS <- "all_antisense_pairs.txt"
GENES_BED <- "genes.bed"
PHYLOP_BW <- "PhyloP.bw"
PHASTCONS_BW <- "PhastCons.bw"

# ChIP-seq
RNAPII_BW <- "external_bigWigs/RNAPII_14D.bw"
H3K4ME3_BW <- "external_bigWigs/H3K4me3_10D.bw"
H3K9AC_BW <- "external_bigWigs/H3K9ac_12D.bw"

# CPC2 results
CPC2_PC <- "protein_coding_2k_cpc2.txt"
CPC2_LNC_ANNOTATED <- "lncrna_annotated_cpc2.txt"
CPC2_LNC_UNANNOTATED <- "lncrna_unannotated_cpc2.txt"

# Differential expression
TOTAL_FC <- "total-quant/Log2FoldChanges.txt"
TOTAL_FDR <- "total-quant/Log2FoldChangesFDR.txt"
TSS_FC <- "tss-quant-qnorm/Log2FoldChanges.txt"
TSS_FDR <- "tss-quant-qnorm/Log2FoldChangesFDR.txt"

# csRNA-seq tag directories (for tagGCcontent.txt)
CSRNA_TAGDIRS <- "csrnaseq_tagDirs"

# Setup ------------------------------------------------------------------------

seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

a11 <- makeTxDbFromGFF(A11_GFF)
a11_tx <- transcripts(a11)

tx_lens <- transcriptLengths(a11)
tx_seqs <- extractTranscriptSeqs(BSgenome.Athaliana.TAIR.TAIR9, a11, use.names = TRUE)

tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tss1 <- import(TSS_SIZE1_BED)

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
totalM <- as.matrix(tssA[, 28:39])
rownames(totalM) <- tssA$TSS
totalMm <- as.matrix(data.frame(row.names = rownames(totalM),
  DS = rowMeans(totalM[, 1:2]),
  S24 = rowMeans(totalM[, 3:4]),
  S72 = rowMeans(totalM[, 5:6]),
  L6 = rowMeans(totalM[, 7:8]),
  L26 = rowMeans(totalM[, 9:10]),
  L57 = rowMeans(totalM[, 11:12])
))

tssAdet <- tssA[!is.na(tssA$RNAseq_DS1), ]

# Figure 4a --------------------------------------------------------------------

tss_counts <- data.frame(
  Type1 = c("mRNA", "Other ncRNA", "lncRNA", "TSS only", "Putative lncRNA"),
  Type2 = c("Araport11", "Araport11", "Araport11", "This study", "This study"),
  Type3 = c("Coding", "Non-coding", "Non-coding", "Non-coding", "Non-coding"),
  N = c(20177, 185, 1108, 6225, 2588)
)

tss_counts$Type1 <- factor(tss_counts$Type1, levels = tss_counts$Type1)
ggplot(tss_counts, aes(x = "", y = N, fill = Type1)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = colorBlindBlack8[-(1:3)]) +
  coord_polar("y", start = 0) +
  theme_void()
ggsave("Figure4/Fig4a_pie1.pdf", width = 2, height = 1)

# Figure 4b --------------------------------------------------------------------

all_as <- readr::read_tsv(ANTISENSE_PAIRS)
tssA$IsAntisense <- tssA$TSS %in% all_as$asTSS

pc_bed <- import(GENES_BED)
tssA$OverlapsSensePC <- structure(overlapsAny(tss1, pc_bed), names = tss1$name)[tssA$TSS]

# Counts for the pie charts
# Num of antisenses
t(table(tssA$IsAntisense, tssA$NearestGeneType, useNA="ifany"))
# Num of sense overlapping
t(table(!tssA$IsAntisense&tssA$OverlapsSensePC, tssA$NearestGeneType, useNA="ifany"))
# Num of intergenic
t(table(!tssA$IsAntisense&!tssA$OverlapsSensePC, tssA$NearestGeneType, useNA="ifany"))

nc_counts_tss <- data.frame(
  type = c("Intergenic", "Sense overlapping", "Antisense overlapping"),
  N = c(2841, 682, 2702)
)

nc_counts_put <- data.frame(
  type = c("Intergenic", "Sense overlapping", "Antisense overlapping"),
  N = c(1689, 0, 899)
)

nc_counts_lnc <- data.frame(
  type = c("Intergenic", "Sense overlapping", "Antisense overlapping"),
  N = c(515, 444, 149)
)

ggplot(nc_counts_tss, aes(x = "", y = N, fill = type)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c("black", "grey50", "grey85")) +
  coord_polar("y", start = 0) +
  theme_void()
ggsave("Figure4/Fig4b_pie2_tss.pdf", height = 0.4, width = 3)

ggplot(nc_counts_put, aes(x = "", y = N, fill = type)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c("black", "grey50", "grey85")) +
  coord_polar("y", start = 0) +
  theme_void()
ggsave("Figure4/Fig4b_pie2_putative.pdf", height = 0.4, width = 3)

ggplot(nc_counts_lnc, aes(x = "", y = N, fill = type)) +
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = c("black", "grey50", "grey85")) +
  coord_polar("y", start = 0) +
  theme_void()
ggsave("Figure4/Fig4b_pie2_lncrna.pdf", height = 0.4, width = 3)

# Figure 4c --------------------------------------------------------------------

phylop <- import(PHYLOP_BW)

phylop_unknown <- colMeans(normalizeToMatrix(phylop, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[is.na(tssA$NearestGeneType)]]))
phylop_putative <- colMeans(normalizeToMatrix(phylop, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "putative_rna"]]))
phylop_other <- colMeans(normalizeToMatrix(phylop, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")]]))
phylop_lncrna <- colMeans(normalizeToMatrix(phylop, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")]]))
phylop_mrna <- colMeans(normalizeToMatrix(phylop, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]]))

pdf("Figure4/Fig4c_phylop.pdf", width = 3.2, height = 3.2)
plot(1:length(phylop_mrna), phylop_mrna, type = "l", lwd = 3,
  xaxt = "n", xlab = "Distance from TSS (bp)", ylab = "PhyloP score",
  ylim = c(0, 0.8), col = colorBlindBlack8[4], bty = "n")
axis(1, c(1, length(phylop_mrna)/2, length(phylop_mrna)), c("-500", "0", "+500"))
lines(1:length(phylop_unknown), phylop_unknown, lwd = 3, col = colorBlindBlack8[7])
lines(1:length(phylop_putative), phylop_putative, lwd = 3, col = colorBlindBlack8[8])
lines(1:length(phylop_lncrna), phylop_lncrna, lwd = 3, col = colorBlindBlack8[6])
dev.off()

rm(phylop)

# Figure 4d, e -----------------------------------------------------------------

tssType <- structure(tssA$NearestGeneType, names = tssA$TSS)
tssType[is.na(tssType)] <- "unknown"
tssType[tssType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")] <- "lncRNA"
tssType[tssType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")] <- "Other_ncRNA"

tss1$Type <- "ncTSS"
tss1$Type[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]] <- "pcTSS"
tss1$Type2 <- tssType[tss1$name]

tss1_pc <- tss1[tss1$Type == "pcTSS"]
tss1_nc <- tss1[tss1$Type == "ncTSS"]
tss1_nc2 <- tss1_nc[tss1_nc$Type2 != "Other_ncRNA"]

rnapii <- import(RNAPII_BW)
h3k4me3 <- import(H3K4ME3_BW)
h3k9ac <- import(H3K9AC_BW)

# Figure 4d
rnapii_pc <- normalizeToMatrix(rnapii, tss1_pc, 2000, value_column = "score", keep = c(0.0, .9))
h3k4me3_pc <- normalizeToMatrix(h3k4me3, tss1_pc, 2000, value_column = "score", keep = c(0, .9))
h3k9ac_pc <- normalizeToMatrix(h3k9ac, tss1_pc, 2000, value_column = "score", keep = c(0, .9))

pdf("Figure4/Fig4d_pc_marks.pdf")
EnrichedHeatmap(rnapii_pc, column_title = "RNAPII",
  col = viridisLite::magma(100),
  pos_line = FALSE,
  raster_quality = 2,
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(yaxis = FALSE, bty = "n",
      gp = gpar(lwd = 2, col = "black")))
) +
EnrichedHeatmap(h3k4me3_pc, column_title = "H3K4me3",
  col = viridisLite::magma(100),
  pos_line = FALSE,
  row_title_gp = gpar(col = colorBlindBlack8[1:3]),
  raster_quality = 2,
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(yaxis = FALSE, bty = "n",
      gp = gpar(lwd = 2, col = colorBlindBlack8[1:3])))
) +
EnrichedHeatmap(h3k9ac_pc, column_title = "H3K9ac",
  col = viridisLite::magma(100),
  pos_line = FALSE,
  raster_quality = 2,
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(yaxis = FALSE, bty = "n",
      gp = gpar(lwd = 2, col = "black")))
)
dev.off()

# Figure 4e
rnapii_nc <- normalizeToMatrix(rnapii, tss1_nc2, 2000, value_column = "score", keep = c(0.0, .9))
h3k4me3_nc <- normalizeToMatrix(h3k4me3, tss1_nc2, 2000, value_column = "score", keep = c(0, .9))
h3k9ac_nc <- normalizeToMatrix(h3k9ac, tss1_nc2, 2000, value_column = "score", keep = c(0, .9))

pdf("Figure4/Fig4e_nc_marks.pdf")
EnrichedHeatmap(rnapii_nc, column_title = "RNAPII", row_split = tss1_nc2$Type2,
  col = viridisLite::magma(100),
  pos_line = FALSE,
  row_title_gp = gpar(col = colorBlindBlack8[1:3]),
  raster_quality = 2,
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(yaxis = FALSE, bty = "n",
      gp = gpar(lwd = 2, col = colorBlindBlack8[1:3])))
) +
EnrichedHeatmap(h3k4me3_nc, column_title = "H3K4me3", row_split = tss1_nc2$Type2,
  col = viridisLite::magma(100),
  pos_line = FALSE,
  row_title_gp = gpar(col = colorBlindBlack8[1:3]),
  raster_quality = 2,
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(yaxis = FALSE, bty = "n",
      gp = gpar(lwd = 2, col = colorBlindBlack8[1:3])))
) +
EnrichedHeatmap(h3k9ac_nc, column_title = "H3K9ac", row_split = tss1_nc2$Type2,
  col = viridisLite::magma(100),
  pos_line = FALSE,
  row_title_gp = gpar(col = colorBlindBlack8[1:3]),
  raster_quality = 2,
  top_annotation = HeatmapAnnotation(
    enriched = anno_enriched(yaxis = FALSE, bty = "n",
      gp = gpar(lwd = 2, col = colorBlindBlack8[1:3])))
)
dev.off()

rm(rnapii, h3k4me3, h3k9ac)

# Supplementary Figure 4a ------------------------------------------------------

tssA_pc <- !is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"
tssC <- as.matrix(tssA[, c("H1", "H2", "R1", "R2", "57H1", "57H2", "GROcap", "CAGE")])
tssCdf <- data.frame(
  Sample = factor(colnames(tssC), levels = colnames(tssC)),
  PctNC = colSums(tssC[!tssA_pc,]) / colSums(tssC) * 100
)
ggplot(tssCdf, aes(Sample, PctNC)) +
  geom_bar(stat = "identity", fill = "black", width = 0.75) +
  scale_y_continuous(expand = c(0, 0), name = "Percent non-coding transcription") +
  xlab(element_blank()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
    legend.position = "none",
    axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.x = element_text(colour = "black", size = unit(5, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    legend.key.size = unit(5, "pt"),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  coord_flip()
ggsave("Figure4/SuppFig4a_percent_nc_transcription.pdf", width = 1.87, height = 1.18)

# Supplementary Figure 4b ------------------------------------------------------

countTATA <- function(x, prop = FALSE) {
  y <- length(unique(scan_sequences(create_motif("TATAWAW"),
    getSeq(BSgenome.Athaliana.TAIR.TAIR9, resize(promoters(x, 80, 0), 60)),
    threshold = 1, threshold.type = "logodds")$sequence))
  if (prop) y <- y / length(x)
  y
}

# unknown TSS: .2583
countTATA(tss1[tss1$name %in% tssA$TSS[is.na(tssA$NearestGeneType)]], T)
# putative lncRNA: .1932
countTATA(tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "putative_rna"]], T)
# Araport11 mRNA: .2377
countTATA(tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]], T)
# Araport11 lncRNA: .2545
countTATA(tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")]], T)
# Araport11 other RNA: .6
countTATA(tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")]], T)

tata_counts <- data.frame(
  Type1 = c("mRNA", "Other ncRNA", "lncRNA", "Unknown TSS", "Putative lncRNA"),
  Type2 = c("Araport11", "Araport11", "Araport11", "This study", "This study"),
  Type3 = c("Coding", "Non-coding", "Non-coding", "Non-coding", "Non-coding"),
  N = c(.2377, .6, .2545, .2583, .1932)
)
tata_counts$Type1 <- factor(tata_counts$Type1, levels = tata_counts$Type1)

ggplot(tata_counts, aes(Type1, N * 100, fill = Type1)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = colorBlindBlack8[-(1:3)]) +
  xlab(element_blank()) +
  scale_y_continuous(limits = c(0, 100), name = "% with TATA", expand = c(0, 0)) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black", size = unit(5, "pt")),
    axis.title.x = element_text(colour = "black", size = unit(5, "pt")),
    panel.grid = element_blank(),
    legend.position = "none"
  )
ggsave("Figure4/SuppFig4b_tata.pdf", width = 1.75, height = 1)

# Supplementary Figure 4c ------------------------------------------------------

tssBases <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, resize(tss1, 9, "center"))
names(tssBases) <- tss1$name

# unknown TSS
view_motifs(create_motif(
    tssBases[tssA$TSS[is.na(tssA$NearestGeneType)]],
  ), min.height = 0.001) +
  theme(axis.text = element_text(size = unit(5, "pt")), axis.title.y = element_text(size=unit(5,"pt")))
ggsave("Figure4/SuppFig4c_ic_unknown.pdf", width = 1, height = .75)
# putative lncRNA
view_motifs(create_motif(
    tssBases[tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "putative_rna"]],
  ), min.height = 0.001) +
  theme(axis.text = element_text(size = unit(5, "pt")), axis.title.y = element_text(size=unit(5,"pt")))
ggsave("Figure4/SuppFig4c_ic_putative.pdf", width = 1, height = .75)
# Araport11 mRNA
view_motifs(create_motif(
    tssBases[tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]],
  ), min.height = 0.001) +
  theme(axis.text = element_text(size = unit(5, "pt")), axis.title.y = element_text(size=unit(5,"pt")))
ggsave("Figure4/SuppFig4c_ic_mrna.pdf", width = 1, height = .75)
# Araport11 lncRNA
view_motifs(create_motif(
    tssBases[tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")]],
  ), min.height = 0.001) +
  theme(axis.text = element_text(size = unit(5, "pt")), axis.title.y = element_text(size=unit(5,"pt")))
ggsave("Figure4/SuppFig4c_ic_lncrna.pdf", width = 1, height = .75)
# Araport11 other RNA
view_motifs(create_motif(
    tssBases[tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")]],
  ), min.height = 0.001) +
  theme(axis.text = element_text(size = unit(5, "pt")), axis.title.y = element_text(size=unit(5,"pt")))
ggsave("Figure4/SuppFig4c_ic_other.pdf", width = 1, height = .75)

# Supplementary Figure 4d ------------------------------------------------------

phast <- import(PHASTCONS_BW)

phast_unknown <- colMeans(normalizeToMatrix(phast, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[is.na(tssA$NearestGeneType)]]))
phast_putative <- colMeans(normalizeToMatrix(phast, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "putative_rna"]]))
phast_other <- colMeans(normalizeToMatrix(phast, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")]]))
phast_lncrna <- colMeans(normalizeToMatrix(phast, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")]]))
phast_mrna <- colMeans(normalizeToMatrix(phast, extend = 500, w = 10, value_column = "score",
  tss1[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]]))

pdf("Figure4/SuppFig4d_phastcons.pdf", width = 3.2, height = 3.2)
plot(1:length(phast_mrna), phast_mrna, type = "l", lwd = 3,
  xaxt = "n", xlab = "Distance from TSS (bp)", ylab = "PhastCons score",
  ylim = c(0, 0.5), col = colorBlindBlack8[4], bty = "n")
axis(1, c(1, length(phast_mrna)/2, length(phast_mrna)), c("-500", "0", "+500"))
lines(1:length(phast_unknown), phast_unknown, lwd = 3, col = colorBlindBlack8[7])
lines(1:length(phast_putative), phast_putative, lwd = 3, col = colorBlindBlack8[8])
lines(1:length(phast_lncrna), phast_lncrna, lwd = 3, col = colorBlindBlack8[6])
dev.off()

rm(phast)

# TSS and transcript sets by type ----------------------------------------------

tss_unknown <- tss1$name[tss1$name %in% tssA$TSS[is.na(tssA$NearestGeneType)]]
tss_putative <- tss1$name[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "putative_rna"]]
tss_other <- tss1$name[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")]]
tss_lncrna <- tss1$name[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")]]
tss_mrna <- tss1$name[tss1$name %in% tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"]]

tx_putative <- unique(tssA$NearestFeature[tssA$TSS %in% tss_putative])
tx_lncrna <- unique(tssA$NearestFeature[tssA$TSS %in% tss_lncrna])
tx_mrna <- unique(tssA$NearestFeature[tssA$TSS %in% tss_mrna])

tx2len <- structure(tx_lens$tx_len, names = tx_lens$tx_name)

# Supplementary Figure 4e ------------------------------------------------------

tx_lengths <- data.frame(
  type = c(rep("Putative lncRNA", length(tx_putative)),
    rep("Araport11 lncRNA", length(tx_lncrna)), rep("Araport11 mRNA", length(tx_mrna))),
  id = c(tx_putative, tx_lncrna, tx_mrna), len = NA_integer_, row.names = NULL)
tx_lengths$len <- tx2len[tx_lengths$id]

# Two-sided Mann-Whitney tests
rstatix::wilcox_test(tx_lengths, len~type)

ggplot(tx_lengths, aes(len, colour = type)) +
  geom_density(alpha = 1, bw = 250, size = 0.75) +
  coord_cartesian(xlim = c(0, 5000)) +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8)], name = element_blank()) +
  xlab("Transcript size (bp)") +
  ylab("Density") +
  theme_bw() +
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
ggsave("Figure4/SuppFig4e_size_dens.pdf", width = 3, height = 1.16)

# Supplementary Figure 4f-i ----------------------------------------------------

tss_max <- data.frame(row.names = NULL,
  id = c(tss_unknown, tss_putative, tss_lncrna, tss_mrna),
  type = c(rep("Unknown TSS", length(tss_unknown)), rep("Putative lncRNA", length(tss_putative)),
    rep("Araport11 lncRNA", length(tss_lncrna)), rep("Araport11 mRNA", length(tss_mrna))),
  max_csrna = NA_real_,
  max_rna = NA_real_
)
tss_max$max_csrna <- rowMaxs(tssMm[tss_max$id, ])
tss_max$max_rna <- rowMaxs(totalMm[tss_max$id, ])
tss_max$cv2_csrna <- rowSds(tssM[tss_max$id,])/rowMeans(tssM[tss_max$id,])
tss_max$cv2_rna <- rowSds(totalM[tss_max$id,])/rowMeans(totalM[tss_max$id,])

total_max <- tss_max[order(tss_max$max_rna, decreasing = T), ]
total_max <- total_max[!is.na(total_max$max_rna), ]
total_max <- total_max[!duplicated(total_max$max_rna), ] # means 0 TPM genes are gone, remember this

# Supplementary Figure 4f
rstatix::wilcox_test(tss_max, max_csrna~type)

ggplot(tss_max, aes(log10(1+max_csrna), colour = type)) +
  geom_density(bw = .5, size = .75) +
  xlab("Max csRNA-seq CPM+1") +
  ylab("Density") +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8, 7)], name = element_blank()) +
  scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)),
    labels = c(1, 10, 100, 1000, 10000)) +
  theme_bw() +
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
ggsave("Figure4/SuppFig4f_tss_expr_dens.pdf", width = 2.75, height = 1.16)

# Supplementary Figure 4g
rstatix::wilcox_test(tss_max, cv2_csrna~type)

ggplot(tss_max, aes(cv2_csrna, colour = type)) +
  geom_density(bw = .2, size = .75) +
  xlab("CV (csRNA-seq)") +
  ylab("Density") +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8, 7)], name = element_blank()) +
  theme_bw() +
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
ggsave("Figure4/SuppFig4g_tss_expr_cv.pdf", width = 2.75, height = 1.16)

# Supplementary Figure 4h
ggplot(total_max, aes(log10(1+max_rna), colour = type)) +
  geom_density(bw = .5, size = .75) +
  xlab("Max RNA-seq TPM+1") +
  ylab("Density") +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8)], name = element_blank()) +
  scale_x_continuous(breaks = log10(c(1, 10, 100, 1000, 10000)),
    labels = c(1, 10, 100, 1000, 10000)) +
  theme_bw() +
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
ggsave("Figure4/SuppFig4h_total_expr_dens.pdf", width = 2.75, height = 1.16)

# Supplementary Figure 4i
rstatix::wilcox_test(total_max, cv2_rna~type)

ggplot(total_max, aes(cv2_rna, colour = type)) +
  geom_density(bw = .2, size = .75) +
  xlab("CV (RNA-seq)") +
  ylab("Density") +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8, 7)], name = element_blank()) +
  theme_bw() +
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
ggsave("Figure4/SuppFig4i_total_expr_cv.pdf", width = 2.75, height = 1.16)

# Supplementary Figure 4j ------------------------------------------------------

tssAdet$type2 <- "protein_coding"
tssAdet$type2[tssAdet$NearestGeneType != "protein_coding"] <- "long_non_coding"
tssAdet$type2[tssAdet$NearestGeneType %in% c(
  "miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")] <- "small_non_coding"
tssAdet$type2[tssAdet$NearestGeneType == "putative_rna"] <- "unannotated_long"

tssAdet2 <- filter(tssAdet, type2 != "small_non_coding")

writeXStringSet(
  tx_seqs[names(tx_seqs) %in% tssAdet2$NearestFeature[tssAdet2$type2=="unannotated_long"]],
  "Figure4/lncrna_unannotated.fa")
writeXStringSet(
  tx_seqs[names(tx_seqs) %in% tssAdet2$NearestFeature[tssAdet2$type2=="long_non_coding"]],
  "Figure4/lncrna_annotated.fa")
writeXStringSet(
  tx_seqs[sample(names(tx_seqs)[names(tx_seqs) %in% tssAdet2$NearestFeature[tssAdet2$type2=="protein_coding"]], 2000)],
  "Figure4/protein_coding_2k.fa")

# CPC2 was run on these files (https://cpc2.gao-lab.org)

cpc_pc <- readr::read_tsv(CPC2_PC)
cpc_ann <- readr::read_tsv(CPC2_LNC_ANNOTATED)
cpc_una <- readr::read_tsv(CPC2_LNC_UNANNOTATED)
cpc_pc$type <- "Araport11 mRNA"
cpc_ann$type <- "Araport11 lncRNA"
cpc_una$type <- "Putative lncRNA"
cpc_all <- rbind(cpc_pc, cpc_ann, cpc_una)

rstatix::wilcox_test(cpc_all, coding_probability~type)

ggplot(cpc_all, aes(coding_probability, colour = type)) +
  geom_density(alpha = 1, bw = .05, size = 0.75) +
  scale_colour_manual(values = colorBlindBlack8[c(6, 4, 8)], name = element_blank()) +
  xlab("Coding probability (CPC2)") +
  ylab("Density") +
  theme_bw() +
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
ggsave("Figure4/SuppFig4j_cpc2.pdf", width = 2.75, height = 1.16)

# Supplementary Figure 13a-d ---------------------------------------------------

total_fc <- as.matrix(read.table(TOTAL_FC))
total_fdr <- as.matrix(read.table(TOTAL_FDR))
tss_fc <- as.matrix(read.table(TSS_FC))
tss_fdr <- as.matrix(read.table(TSS_FDR))

# csRNA-seq, non-coding TSSs
cols1 <- rep("black", length(c(tss_lncrna,tss_unknown,tss_putative)))
cols1[tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"]>=1 &
    tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"]<.05] <- "red"
cols1[tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"]<=-1 &
    tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"]<.05] <- "blue"
table(cols1)

cols2 <- rep("black", length(c(tss_lncrna,tss_unknown,tss_putative)))
cols2[tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"]>=1 &
    tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"]<.05] <- "red"
cols2[tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"]<=-1 &
    tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"]<.05] <- "blue"
table(cols2)

# Supplementary Figure 13a
ggplot(data.frame(
  fc = tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"],
  fdr = pmin(10,-log10(tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"])),
  cols = cols2
), aes(fc, fdr, colour = cols)) +
  geom_point(size=.2) +
  scale_colour_identity() +
  xlab("Log2(Fold change WT:hen2)") +
  ylab("-log10(FDR)") +
  xlim(c(-10, 10)) +
  theme_bw() +
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
ggsave("Figure4/SuppFig13a_hen2_ncrna_de.pdf", width = 1.8, height = 1.25)

# Supplementary Figure 13b
ggplot(data.frame(
  fc = tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"],
  fdr = pmin(10,-log10(tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"])),
  cols = cols1
), aes(fc, fdr, colour = cols)) +
  geom_point(size=.2) +
  scale_colour_identity() +
  xlab("Log2(Fold change WT:rrp4)") +
  ylab("-log10(FDR)") +
  xlim(c(-10, 10)) +
  theme_bw() +
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
ggsave("Figure4/SuppFig13b_rrp4_ncrna_de.pdf", width = 1.8, height = 1.25)

# RNA-seq, Araport11 and putative lncRNAs
tt <- c(tx_lncrna, tx_putative)
tt <- tt[tt %in% rownames(total_fc)]
cols3 <- rep("black", length(tt))
cols3[total_fc[tt, "L57_R"]>=1 &
    total_fdr[tt, "L57_R"]<.05] <- "red"
cols3[total_fc[tt, "L57_R"]<=-1 &
    total_fdr[tt, "L57_R"]<.05] <- "blue"
table(cols3)

cols4 <- rep("black", length(tt))
cols4[total_fc[tt, "L57_H"]>=1 &
    total_fdr[tt, "L57_H"]<.05] <- "red"
cols4[total_fc[tt, "L57_H"]<=-1 &
    total_fdr[tt, "L57_H"]<.05] <- "blue"
table(cols4)

# Supplementary Figure 13c
ggplot(data.frame(
  fc = pmin(10, total_fc[tt, "L57_H"]),
  fdr = pmin(10,-log10(total_fdr[tt, "L57_H"])),
  cols = cols4
), aes(fc, fdr, colour = cols)) +
  geom_point(size=.2) +
  scale_colour_identity() +
  ylim(c(0, 10)) +
  xlim(c(-10, 10)) +
  xlab("Log2(Fold change WT:hen2)") +
  ylab("-log10(FDR)") +
  theme_bw() +
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
ggsave("Figure4/SuppFig13c_hen2_ncrna_de_rnaseq.pdf", width = 1.8, height = 1.25)

# Supplementary Figure 13d
ggplot(data.frame(
  fc = pmin(10, total_fc[tt, "L57_R"]),
  fdr = pmin(10,-log10(total_fdr[tt, "L57_R"])),
  cols = cols3
), aes(fc, fdr, colour = cols)) +
  geom_point(size=.2) +
  scale_colour_identity() +
  ylim(c(0, 10)) +
  xlim(c(-10, 10)) +
  xlab("Log2(Fold change WT:rrp4)") +
  ylab("-log10(FDR)") +
  theme_bw() +
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
ggsave("Figure4/SuppFig13d_rrp4_ncrna_de_rnaseq.pdf", width = 1.8, height = 1.25)

# Supplementary Figure 13e -----------------------------------------------------

phyper(969, 297+969, 522+969, 1788, lower.tail = F)
phyper(351, 107+351, 40+351, 498, lower.tail = F)
phyper(58, 911+58, 293+58, 1262, lower.tail = F)

table(cols3,cols4)

pdf("Figure4/SuppFig13e_hen_rrp4_ov_rnaseq.pdf", width = 1, height = 1)
plot(euler(c(h = 458-351, r = 391-351, "h&r" = 351)),
  edges = FALSE, fills = list(fill = c("grey25", "grey75")))
dev.off()

table(cols1,cols2)

pdf("Figure4/SuppFig13e_hen_rrp4_ov_csrnaseq.pdf", width = 1, height = 1)
plot(euler(c(h = 1266-969, r = 1491-969, "h&r" = 969)),
  edges = FALSE, fills = list(fill = c("grey25", "grey75")))
dev.off()

# Supplementary Figure 13f, g --------------------------------------------------

# Map each transcript to its highest-expressed (L57 RNA-seq) TSS
f2t <- structure(tssA$TSS[!is.na(tssA$NearestFeature)], names = tssA$NearestFeature[!is.na(tssA$NearestFeature)])
f2t <- f2t[unname(f2t)%in%rownames(totalMm)]
f2t <- f2t[order(totalMm[f2t,"L57"],decreasing=T)]
f2t <- f2t[!duplicated(names(f2t))]
f2t <- f2t[names(f2t)%in%rownames(total_fc)]

# Exosome-insensitive: log2FC < 1 and FDR >= 0.05 in both mutants
cols5 <- rep("black", length(c(tss_lncrna,tss_unknown,tss_putative)))
cols5[tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"]<1 &
    tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_R"]>=.05] <- "red"

cols6 <- rep("black", length(c(tss_lncrna,tss_unknown,tss_putative)))
cols6[tss_fc[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"]<1 &
    tss_fdr[c(tss_lncrna,tss_unknown,tss_putative), "L57_H"]>=.05] <- "red"

cols7 <- rep("black", length(tt))
cols7[total_fc[tt, "L57_R"]<1 &
    total_fdr[tt, "L57_R"]>=.05] <- "red"

cols8 <- rep("black", length(tt))
cols8[total_fc[tt, "L57_H"]<1 &
    total_fdr[tt, "L57_H"]>=.05] <- "red"

cc <- c(tss_lncrna,tss_unknown,tss_putative)

# Exosome-sensitive: upregulated in both mutants (red in both volcano plots)
tss_exo <- c(tss_lncrna,tss_unknown,tss_putative)[cols1=="red" & cols2=="red"]
total_exo <- tt[cols3=="red" & cols4=="red"]

tss_ins <- cc[cols5=="red" & cols6=="red"]
total_ins <- tt[cols7=="red" & cols8=="red"]

tss_ins3 <- tss_ins[tssMm[tss_ins, "L57"]>=1]

readr::write_lines(tss_exo, "Figure4/tss_exosome_sensitive.txt")
readr::write_lines(total_exo, "Figure4/transcripts_exosome_sensitive.txt")
readr::write_lines(tss_ins, "Figure4/tss_exosome_insensitive.txt")
readr::write_lines(total_ins, "Figure4/transcripts_exosome_insensitive.txt")

noNA <- function(x) x[!is.na(x)]

total_ins4 <- noNA(tssA$TSS[tssA$NearestFeature %in% total_ins])
total_ins4 <- total_ins4[totalMm[total_ins4, "L57"]>=.1]

getGC <- function(x) {
  s <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, resize(x, 2000, "center"))
  m <- t(do.call(cbind, strsplit(as.character(s), "", T)))
  m <- apply(m, 2, function(y) sum(toupper(y) %in% c("G", "C"))) / length(x)
  m
}

zm <- function(x) zoo::rollmean(x, 20, na.pad=T)

# Supplementary Figure 13f (csRNA-seq)
pdf("Figure4/SuppFig13f_tss_cg_exo.pdf", width = 3, height = 2.9)
plot(1:2000, zm(getGC(tss1[tss1$name%in%cc[cc%in%f2t]])), lwd=2,type = "l", ylim = c(.29, .47),xaxt="n",xlab="Distance from TSS",ylab="CG content",bty="n")
axis(1,c(1,1000,2000),c("-1kbp","0","+1kbp"))
lines(1:2000, zm(getGC(tss1[tss1$name%in%tss_exo[tss_exo%in%f2t]])), lwd=2,col = "red")
lines(1:2000, zm(getGC(tss1[tss1$name%in%tss_ins3[tss_ins3%in%f2t]])), lwd=2,col = "blue")
dev.off()

# Supplementary Figure 13g (RNA-seq)
pdf("Figure4/SuppFig13g_total_cg_exo.pdf", width = 3, height = 2.9)
plot(1:2000, zm(getGC(tss1[tss1$name%in%f2t[tt]])), lwd=2,type = "l", ylim = c(.29, .47),xaxt="n",xlab="Distance from TSS",ylab="CG content",bty="n")
axis(1,c(1,1000,2000),c("-1kbp","0","+1kbp"))
lines(1:2000, zm(getGC(tss1[tss1$name%in%f2t[total_exo]])), lwd=2,col = "red")
lines(1:2000, zm(getGC(tss1[tss1$name%in%total_ins4])), lwd=2,col = "blue")
dev.off()

# Supplementary Figure 13h -----------------------------------------------------

lib_samples <- c("DS1", "DS2", "S241", "S242", "S721", "S722", "L61", "L62",
  "L261", "L262", "L571", "L572", "hen21", "hen22", "rrp41", "rrp42")
GCs <- lapply(lib_samples,
  function(x) readr::read_tsv(file.path(CSRNA_TAGDIRS, paste0(x, "-tagDir"), "tagGCcontent.txt")))
names(GCs) <- lib_samples

pdf("Figure4/SuppFig13h_library_gc.pdf", width = 5.5, height = 4)
plot(GCs[["DS1"]][[1]], GCs[["DS1"]][[3]], type = "l", lwd = 3,
  col = colorBlindBlack8[1], ylim = c(0, 8), bty = "n",
  xlab = "GC %", ylab = "Normalized fraction (PDF)"
)
lines(GCs[["DS2"]][[1]], GCs[["DS2"]][[3]], lwd = 3, col = colorBlindBlack8[1])

lines(GCs[["S241"]][[1]], GCs[["S241"]][[3]], lwd = 3, col = colorBlindBlack8[2])
lines(GCs[["S242"]][[1]], GCs[["S242"]][[3]], lwd = 3, col = colorBlindBlack8[2])

lines(GCs[["S721"]][[1]], GCs[["S721"]][[3]], lwd = 3, col = colorBlindBlack8[3])
lines(GCs[["S722"]][[1]], GCs[["S722"]][[3]], lwd = 3, col = colorBlindBlack8[3])

lines(GCs[["L61"]][[1]], GCs[["L61"]][[3]], lwd = 3, col = colorBlindBlack8[4])
lines(GCs[["L62"]][[1]], GCs[["L62"]][[3]], lwd = 3, col = colorBlindBlack8[4])

lines(GCs[["L261"]][[1]], GCs[["L261"]][[3]], lwd = 3, col = colorBlindBlack8[5])
lines(GCs[["L262"]][[1]], GCs[["L262"]][[3]], lwd = 3, col = colorBlindBlack8[5])

lines(GCs[["L571"]][[1]], GCs[["L571"]][[3]], lwd = 3, col = colorBlindBlack8[6])
lines(GCs[["L572"]][[1]], GCs[["L572"]][[3]], lwd = 3, col = colorBlindBlack8[6])

lines(GCs[["hen21"]][[1]], GCs[["hen21"]][[3]], lwd = 3, col = colorBlindBlack8[7])
lines(GCs[["hen22"]][[1]], GCs[["hen22"]][[3]], lwd = 3, col = colorBlindBlack8[7])

lines(GCs[["rrp41"]][[1]], GCs[["rrp41"]][[3]], lwd = 3, col = colorBlindBlack8[8])
lines(GCs[["rrp42"]][[1]], GCs[["rrp42"]][[3]], lwd = 3, col = colorBlindBlack8[8])

legend("topright", legend = c("DS", "S24", "S72", "L6", "L26", "L57", "hen2-4", "rrp4-2"),
  col = colorBlindBlack8, lwd = 3, lty = 1, bty = "n")
dev.off()
