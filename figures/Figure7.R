# Benjamin J.M. Tremblay
#
# Figure 7 and Supplementary Figures 9, 10
#
# Package requirements: rtracklayer, matrixStats, ComplexHeatmap, circlize,
# viridisLite, magick, ggplot2, dplyr, tidyr, readxl, rstatix, readr, zoo,
# genomation, EnrichedHeatmap, reshape2, gprofiler2,
# BSgenome.Athaliana.TAIR.TAIR9, universalmotif

library(rtracklayer)
library(matrixStats)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(ggplot2)
library(dplyr)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

source("Functions.R")

dir.create("Figure7", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_SIZE1_BED <- "tss.size1.bed"
BIDIR_PAIRS <- "Bidir_all_pairs.tsv"
DETECTED_TX_BED <- "detected.transcripts.bed"
ACR_500BP_BED <- "peaks_final_500bp.bed"
ATAC_CLUSTERS <- "ATACseq_clustering3/Clusters_anno.tsv"
ENHANCER_ACTIVITY <- "Final_enhancer_activity.tsv"
ENHANCER_GENE_LINKS <- "Enhancer_gene_links.tsv"

# RT-qPCR and luciferase assays
QPCR_SLP2 <- "SLP2_qpcr.xlsx"
LUC_ASSAY <- "LUC_assay.xlsx"

# Read density
CSRNA_L57_PLUS_BW <- "L57.merged.csRNA.plus.bw"
CSRNA_L57_MINUS_BW <- "L57.merged.csRNA.minus.bw"
CSRNA_ALL_PLUS_BG <- "csrnaseq_bedGraphs/all.plus.bedGraph.gz"
CSRNA_ALL_MINUS_BG <- "csrnaseq_bedGraphs/all.minus.bedGraph.gz"
HEN2_RNA_PLUS_BW <- "hen2.merged.RNA.plus.bw"
HEN2_RNA_MINUS_BW <- "hen2.merged.RNA.minus.bw"
ATAC_DS_BW <- "atacseq_bw/DS.bw"
ATAC_L6_BW <- "atacseq_bw/L6.bw"
ATAC_L26_BW <- "atacseq_bw/L26.bw"
ATAC_L57_BW <- "atacseq_bw/L57.bw"
RNAPII_BW <- "external_bigWigs/RNAPII_14D.bw"
H3K4ME1_BW <- "external_bigWigs/H3K4me1_14D.bw"
H3K4ME3_BW <- "external_bigWigs/H3K4me3_10D.bw"
H3K9AC_BW <- "external_bigWigs/H3K9ac_12D.bw"
H3K27ME3_BW <- "external_bigWigs/H3K27me3_7D.bw"
H3K36ME3_BW <- "external_bigWigs/H3K36me3_10D.bw"
H2AZ_BW <- "external_bigWigs/H2AZ_10D.bw"
H2AK121UB_BW <- "external_bigWigs/H2AK121ub_7D.bw"
MNASE_BW <- "Ath_leaf_NPS.bw"
PHYLOP_BW <- "PhyloP.bw"
PHASTCONS_BW <- "PhastCons.bw"

# Shared data ------------------------------------------------------------------

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tss1 <- import(TSS_SIZE1_BED)
det_tx <- import(DETECTED_TX_BED)

tss_pc <- tssA$TSS[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType=="protein_coding"]

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

bd_all <- readr::read_tsv(BIDIR_PAIRS)

# Bidirectional non-coding (ncTSS-ncTSS) pairs within 1 kbp
ncnc <- bd_all[bd_all$TypeV2=="nc-nc" & bd_all$dist <= 1000, ]

# Midpoint between the two ncTSSs, sorted by inter-TSS distance
ncncBedMid <- with(ncnc, GRanges(seq1,
  IRanges(floor((pos2+pos1)/2), floor((pos2+pos1)/2)),
  len = pos2 - pos1,
  TSS1 = TSS1, TSS2 = TSS2,
))
ncncBedMid$len <- abs(ncncBedMid$len)
ncncBedMid <- ncncBedMid[order(ncncBedMid$len)]
ncncBedMid$Intragenic <- overlapsAny(ncncBedMid, det_tx)
ncncBedMid500 <- ncncBedMid[ncncBedMid$len<=500]

# Putative enhancers and their activities (csRNA-seq signal)
enhA <- readr::read_tsv(ENHANCER_ACTIVITY)
enhAll500 <- GRanges(enhA$seqnames, IRanges(enhA$start, enhA$end), name = enhA$name)
enhActivity <- as.matrix(enhA[, c("DS", "S24", "S72", "L6", "L26", "L57")])
rownames(enhActivity) <- enhA$name

# Figure 7a --------------------------------------------------------------------

# ATAC-seq (L57)

atac_57 <- import(ATAC_L57_BW)

atac_bd <- get_unstranded_signal(atac_57, promoters(ncncBedMid, 1000, 1000), keep=c(0.,.9),
  return.matrix=TRUE)

pdf("Figure7/Fig7a_bd_atac.pdf", width = 2, height = 0.8)
Heatmap(atac_bd,
  show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
  col=viridisLite::plasma(100,begin=0),
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

ncncBedMidS <- ncncBedMid
strand(ncncBedMidS) <- "+"

# RNA-seq (hen2-4)

hen2 <- read_signal(
  pos = HEN2_RNA_PLUS_BW,
  neg = HEN2_RNA_MINUS_BW)

hen2_bd <- get_stranded_signal(hen2, keep = c(0., .95),
  windows = promoters(ncncBedMidS, 1000, 1000), return.matrix = TRUE, merge = TRUE)

pdf("Figure7/Fig7a_bd_rna_hen2.pdf", width = 2, height = 0.8)
Heatmap(hen2_bd,
  show_row_names = FALSE, show_column_names = FALSE,
  use_raster=T,raster_quality=4,
  col = circlize::colorRamp2(c(-.1, 0, .1),
      c("red", "white", "blue")),
  raster_by_magick = TRUE,
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# MNase-seq

mnase <- import(MNASE_BW)
seqlevels(mnase) <- as.character(1:5)
mnase_bd <- get_unstranded_signal(mnase, keep = c(0., .90),
  windows = promoters(ncncBedMidS, 1000, 1000), return.matrix = TRUE)

pdf("Figure7/Fig7a_bd_mnase.pdf", width = 2, height = .8)
Heatmap(mnase_bd,
  show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
  col=viridisLite::plasma(100,begin=0),
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# RNAPII ChIP-seq

rnap <- import(RNAPII_BW)

rnap_bd <- get_unstranded_signal(rnap, promoters(ncncBedMidS, 1000, 1000), keep=c(0.,.9),
  return.matrix=TRUE)

pdf("Figure7/Fig7a_bd_rnapii.pdf", width = 2, height = .8)
Heatmap(rnap_bd,
  show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
  col=viridisLite::plasma(100,begin=0),
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# csRNA-seq (L57)

csrna_57 <- read_signal(
  pos = CSRNA_L57_PLUS_BW,
  neg = CSRNA_L57_MINUS_BW
)

csrna_bd <- get_stranded_signal(
  csrna_57,
  windows = promoters(ncncBedMidS, 1000, 1000),
  row.norm=T,
  return.matrix = T,
  merge = T
)

pdf("Figure7/Fig7a_bd_csrna.pdf", width = 2, height = .8)
Heatmap(t(apply(csrna_bd, 1, nn0)),
  show_row_names = FALSE, show_column_names = FALSE,
  use_raster=T,raster_quality=4,
  col = circlize::colorRamp2(c(-.1, 0, .1),
      c("red", "white", "blue")),
  raster_by_magick = TRUE,
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# Figure 7f --------------------------------------------------------------------

# Top and bottom 500 enhancers by L57 activity
enh57act <- resize(enhAll500[enhAll500$name%in%rownames(enhActivity[order(enhActivity[,"L57"],decreasing=T),])[1:500]],2000,"center")
enh57ina <- resize(enhAll500[enhAll500$name%in%rownames(enhActivity[order(enhActivity[,"L57"],decreasing=F),])[1:500]],2000,"center")
enh57act$Status <- "Active"
enh57ina$Status <- "Inactive"
enh57 <- c(enh57act,enh57ina)

rnapii <- import(RNAPII_BW)
h3k4me1 <- import(H3K4ME1_BW)
h3k4me3 <- import(H3K4ME3_BW)
h3k9ac <- import(H3K9AC_BW)
h3k27me3 <- import(H3K27ME3_BW)
h3k36me3 <- import(H3K36ME3_BW)
h2az <- import(H2AZ_BW)
h2ak121ub <- import(H2AK121UB_BW)
atac_l57 <- import(ATAC_L57_BW)
mnase <- import(MNASE_BW)
seqlevels(mnase) <- as.character(1:5)

rnapii_act <- normalizeToMatrix(rnapii, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h3k4me1_act <- normalizeToMatrix(h3k4me1, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h3k4me3_act <- normalizeToMatrix(h3k4me3, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h3k9ac_act <- normalizeToMatrix(h3k9ac, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h3k27me3_act <- normalizeToMatrix(h3k27me3, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h3k36me3_act <- normalizeToMatrix(h3k36me3, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h2az_act <- normalizeToMatrix(h2az, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
h2ak121ub_act <- normalizeToMatrix(h2ak121ub, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
atac_l57_act <- normalizeToMatrix(atac_l57, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))
mnase_act <- normalizeToMatrix(mnase, resize(enh57,1,"center"),
  extend = 2000, value_column = "score", keep = c(0.,.9))

rm(rnapii, h3k4me1, h3k4me3, h3k9ac, h3k27me3, h3k36me3, h2az, h2ak121ub, atac_l57)

pdf("Figure7/Fig7f_enhancer_marks.pdf", width = 10.85, height = 3.5)
EnrichedHeatmap(atac_l57_act, row_split = enh57$Status, column_title="Accessibility",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(rnapii_act, row_split = enh57$Status, column_title="RNAPII",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(mnase_act, row_split = enh57$Status, column_title="Nucleosome\noccupancy",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h3k9ac_act, row_split = enh57$Status, column_title="H3K9ac",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h3k4me1_act, row_split = enh57$Status, column_title="H3K4me1",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h3k4me3_act, row_split = enh57$Status, column_title="H3K4me3",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h3k36me3_act, row_split = enh57$Status, column_title="H3K36me3",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h2az_act, row_split = enh57$Status, column_title="H2AZ",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h2ak121ub_act, row_split = enh57$Status, column_title="H2AK121ub",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
) +
EnrichedHeatmap(h3k27me3_act, row_split = enh57$Status, column_title="H3K27me3",
  row_title_gp = gpar(col = 1:2), top_annotation=HeatmapAnnotation(enriched=anno_enriched(yaxis=F,height=unit(1,"cm"),gp=gpar(col=1:2,lwd=2))),
  use_raster=T,raster_quality=2,
  col = viridisLite::magma(100), pos_line = FALSE
)
dev.off()

# Supplementary Figure 9a ------------------------------------------------------

pdf("Figure7/Supp9a_bidir_size_density.pdf", width = 4, height = 3.5)
plot(density(ncncBedMid500$len,bw=25), lwd = 3, xlim = c(0, 500),
  ylim = c(0, 0.0045), main = NA, xlab = "Distance between ncTSSs (bp)",
  ylab = "Density of bidirectional promoter sizes"
)
abline(v=211,lty=2,lwd=2,col="grey50")
dev.off()

density(ncncBedMid500$len,bw=25)$x[which.max(density(ncncBedMid500$len,bw=25)$y)]

# Supplementary Figure 9b ------------------------------------------------------

# Promoter midpoints on both strands
bdMid <- ncncBedMid500
strand(bdMid) <- "+"
tmp <- bdMid
strand(tmp) <- "-"
bdMid <- sort(c(bdMid,tmp))
rm(tmp)

# Individual ncTSSs
bd_tss <- tss1[tss1$name %in% c(bdMid$TSS1, bdMid$TSS2)]
bd_tss$Intragenic <- bd_tss$name %in% c(bdMid$TSS1[bdMid$Intragenic],bdMid$TSS2[bdMid$Intragenic])

bd_atac_ds <- get_unstranded_signal(ATAC_DS_BW,promoters(bd_tss,1000,0))
bd_atac_l6 <- get_unstranded_signal(ATAC_L6_BW,promoters(bd_tss,1000,0))
bd_atac_l26 <- get_unstranded_signal(ATAC_L26_BW,promoters(bd_tss,1000,0))
bd_atac_l57 <- get_unstranded_signal(ATAC_L57_BW,promoters(bd_tss,1000,0))

which.max(bd_atac_ds$Avg)
which.max(bd_atac_l6$Avg)
which.max(bd_atac_l26$Avg)
which.max(bd_atac_l57$Avg)

pdf("Figure7/Supp9b_bidir_atac.pdf", width = 4, height = 3.5)
plot(-bd_atac_ds$Coord[501:1000],bd_atac_ds$Avg[501:1000],type="l",lwd=3,
  xaxt = "n", xlab = "Distance from ncTSS (bp)", ylab = "ATAC-seq signal probability",
  col=colorBlindBlack8[1], ylim=0:1)
axis(1,-c(500,600,700,800,900,1000),c(-500,-400,-300,-200,-100,0))
lines(-bd_atac_l6$Coord[501:1000], bd_atac_l6$Avg[501:1000],lwd=3,col=colorBlindBlack8[4])
lines(-bd_atac_l26$Coord[501:1000], bd_atac_l26$Avg[501:1000],lwd=3,col=colorBlindBlack8[5])
lines(-bd_atac_l57$Coord[501:1000], bd_atac_l57$Avg[501:1000],lwd=3,col=colorBlindBlack8[6])
abline(v=c(-911),lty=2,lwd=2,col="grey50")
dev.off()

# Supplementary Figure 9c ------------------------------------------------------

bd_inter_phylop <- get_unstranded_signal(PHYLOP_BW, promoters(bd_tss[!bd_tss$Intragenic],1000,500), row.norm=F)
bd_intra_phylop <- get_unstranded_signal(PHYLOP_BW, promoters(bd_tss[bd_tss$Intragenic],1000,500), row.norm=F)

pdf("Figure7/Supp9c_phylop.pdf", width = 4, height = 3.5)
plot(bd_inter_phylop$Coord,(bd_inter_phylop$Avg),type="l",lwd=3,
  xaxt = "n", xlab = "Distance from ncTSS (bp)", ylab = "PhyloP score",
  col="black", ylim = c(0.3,1)
)
axis(1,c(1,500,1000,1500),c(-1000,-500,0,500))
lines(bd_intra_phylop$Coord,(bd_intra_phylop$Avg),lwd=3,col="grey75")
abline(v=1000,lty=2,lwd=2)
dev.off()

# Supplementary Figure 9d ------------------------------------------------------

bd_inter_phast <- get_unstranded_signal(PHASTCONS_BW, promoters(bd_tss[!bd_tss$Intragenic],1000,500), row.norm=F)
bd_intra_phast <- get_unstranded_signal(PHASTCONS_BW, promoters(bd_tss[bd_tss$Intragenic],1000,500), row.norm=F)

pdf("Figure7/Supp9d_phastcons.pdf", width = 4, height = 3.5)
plot(bd_inter_phast$Coord,(bd_inter_phast$Avg),type="l",lwd=3,
  xaxt = "n", xlab = "Distance from ncTSS (bp)", ylab = "PhastCons score",
  col="black", ylim = c(0.15,0.45)
)
axis(1,c(1,500,1000,1500),c(-1000,-500,0,500))
lines(bd_intra_phast$Coord,(bd_intra_phast$Avg),lwd=3,col="grey75")
abline(v=1000,lty=2,lwd=2)
dev.off()

# Supplementary Figure 9e ------------------------------------------------------

get_freqs <- function(x) {
  list(
    A = apply(x, 2, function(x) sum(x=="A")) / nrow(x),
    C = apply(x, 2, function(x) sum(x=="C")) / nrow(x),
    G = apply(x, 2, function(x) sum(x=="G")) / nrow(x),
    T = apply(x, 2, function(x) sum(x=="T")) / nrow(x)
  )
}

bd <- bd_all[bd_all$dist<=500 & bd_all$TypeV2=="nc-nc", ]
bdgr <- GRanges(bd$seq1, IRanges(pmin(bd$pos1,bd$pos2),pmax(bd$pos1,bd$pos2)),strand="+")
bd_u <- unique(c(bd$TSS1, bd$TSS2))

bdfa <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, resize(bdgr, 1000, "center"))
bdfaM <- as.matrix(t(as.data.frame(strsplit(as.character(bdfa), "", fixed=T))))

f_bd <- get_freqs(bdfaM)

pdf("Figure7/Supp9e_bidir_base_freq.pdf", width = 3.8, height = 3.5)
plot(c((-499):0, 1:500), f_bd$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-500, 0, 500))
axis(2, c(0, .5, 1))
lines(c((-499):0, 1:500), f_bd$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-499):0, 1:500), f_bd$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-499):0, 1:500), f_bd$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 9f ------------------------------------------------------

tss1fa2kb <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, promoters(tss1[tss1$name %in% bd_u], 500, 500))
names(tss1fa2kb) <- tss1$name[tss1$name %in% bd_u]
tss1fa2kbM <- as.matrix(t(as.data.frame(strsplit(as.character(tss1fa2kb), "", fixed=T))))

f_bd2 <- get_freqs(tss1fa2kbM[bd_u,,drop=F])

pdf("Figure7/Supp9f_bidir_tss_base_freq.pdf", width = 3.8, height = 3.5)
plot(c((-499):0, 1:500), f_bd2$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-500, 0, 500))
axis(2, c(0, .5, 1))
lines(c((-499):0, 1:500), f_bd2$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-499):0, 1:500), f_bd2$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-499):0, 1:500), f_bd2$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 9g ------------------------------------------------------

bd_intra_h3k4 <- get_unstranded_signal(H3K4ME3_BW,
  resize(bdMid[bdMid$Intragenic],2000,"center"),row.norm=F)
bd_inter_h3k4 <- get_unstranded_signal(H3K4ME3_BW,
  resize(bdMid[!bdMid$Intragenic],2000,"center"),row.norm=F)

pdf("Figure7/Supp9g_h3k4me3.pdf", width = 4, height = 3.5)
plot(bd_inter_h3k4$Coord,bd_inter_h3k4$Avg,type="l",lwd=3,
  xaxt = "n", xlab = "Distance from promoter midpoint (bp)", ylab = "H3K4me3 (RPM)",
  col="black", ylim = c(5,14)
)
axis(1,c(1,500,1000,1500,2000),c(-1000,-500,0,500,1000))
lines(bd_intra_h3k4$Coord,bd_intra_h3k4$Avg,lwd=3,col="grey75")
abline(v=1000,lty=2,lwd=2)
dev.off()

# Supplementary Figure 9h ------------------------------------------------------

bd_intra_h3k9 <- get_unstranded_signal(H3K9AC_BW,
  resize(bdMid[bdMid$Intragenic],2000,"center"),row.norm=F)
bd_inter_h3k9 <- get_unstranded_signal(H3K9AC_BW,
  resize(bdMid[!bdMid$Intragenic],2000,"center"),row.norm=F)

pdf("Figure7/Supp9h_h3k9ac.pdf", width = 4, height = 3.5)
plot(bd_inter_h3k9$Coord,bd_inter_h3k9$Avg,type="l",lwd=3,
  xaxt = "n", xlab = "Distance from promoter midpoint (bp)", ylab = "H3K9ac (RPM)",
  col="black", ylim = c(5,26)
)
axis(1,c(1,500,1000,1500,2000),c(-1000,-500,0,500,1000))
lines(bd_intra_h3k9$Coord,bd_intra_h3k9$Avg,lwd=3,col="grey75")
abline(v=1000,lty=2,lwd=2)
dev.off()

# Supplementary Figure 9i ------------------------------------------------------

ncnc500 <- filter(ncnc, dist <= 500)
ncnc500$MaxCPM1 <- rowMaxs(tssMm[ncnc500$TSS1, ])
ncnc500$MaxCPM2 <- rowMaxs(tssMm[ncnc500$TSS2, ])
ncnc500$MaxRatio <-
  (1 + ifelse(ncnc500$MaxCPM1 > ncnc500$MaxCPM2, ncnc500$MaxCPM1, ncnc500$MaxCPM2)) /
  (1 + ifelse(ncnc500$MaxCPM1 > ncnc500$MaxCPM2, ncnc500$MaxCPM2, ncnc500$MaxCPM1))

ncnc500$CorCat <- "Neutral"
ncnc500$CorCat[ncnc500$PCC >= .25] <- "Positive"
ncnc500$CorCat[ncnc500$PCC <= -.25] <- "Negative"
ncnc500$Log2MaxRatio <- log2(ncnc500$MaxRatio)

# Two-sided Mann-Whitney tests
rstatix::wilcox_test(ncnc500, MaxRatio~CorCat)

# Medians used for the dashed lines
quantile(ncnc500$MaxRatio[ncnc500$PCC>=.25], p = c(.25, .5, .75))
quantile(ncnc500$MaxRatio[ncnc500$PCC<.25 & ncnc500$PCC> -.25], p = c(.25, .5, .75))
quantile(ncnc500$MaxRatio[ncnc500$PCC<= -.25], p = c(.25, .5, .75))

ggplot(ncnc500, aes(log2(MaxRatio), colour = CorCat)) +
  ylab("Density") +
  geom_vline(xintercept = log2(c(2.081)), colour = "red", lty = 2) +
  geom_vline(xintercept = log2(c(3.379)), colour = "black", lty = 2) +
  geom_vline(xintercept = log2(c(4.252)), colour = "blue", lty = 2) +
  geom_density(bw = 1) +
  coord_cartesian(xlim = c(0, 10)) +
  scale_x_continuous(name = "Ratio of max noncoding/noncoding expression",
    breaks = c(-6, -3, 0, 3, 6, 9),
    labels = c("1 / 64", "1 / 8",  "1", "8", "64", "512")
  ) +
  scale_colour_manual(values = c("blue", "black", "red")) +
  theme_bw() +
  theme(
    legend.position = "none",
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
ggsave("Figure7/Supp9i_noncoding_noncoding_ratio.pdf", height = 1.33, width = 2)

# Supplementary Figure 9k ------------------------------------------------------

qpcr <- readxl::read_xlsx(QPCR_SLP2)
qpcr <- tidyr::gather(qpcr, "Primer", "FC_RBP45B", 5:7)
qpcr_s <- qpcr %>%
  group_by(Strand, Timepoint, Genotype, Primer) %>%
  summarise(Mean = mean(FC_RBP45B), Sd = sd(FC_RBP45B))
qpcr$Group <- paste0(qpcr$Timepoint, "_", qpcr$Genotype, "_", qpcr$Strand)

qpcr_s$Strand <- factor(qpcr_s$Strand, levels = c("Sense", "Antisense"))

# One-way ANOVA with Tukey's HSD (P1 primer pair)
rstatix::tukey_hsd(aov(FC_RBP45B~Group,data=filter(qpcr,Primer=="SLP2_R"&Strand=="Sense"&Timepoint=="DS")))
rstatix::tukey_hsd(aov(FC_RBP45B~Group,data=filter(qpcr,Primer=="SLP2_R"&Strand=="Antisense"&Timepoint=="DS")))

ggplot(filter(qpcr_s, Timepoint=="DS"&Primer=="SLP2_R"), aes(Strand, Mean*100, fill = Genotype)) +
  geom_errorbar(aes(ymin = (Mean*100) / 2, ymax = (Mean + Sd)*100),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.3), name = "FC to RBP45B x 102") +
  scale_fill_manual(values = c("black", "grey25", "grey75"), name = element_blank()) +
  theme_bw() +
  theme(
    legend.text = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
ggsave("Figure7/Supp9k_qpcr_slp2_p1.pdf", width = 2.5, height = 1.2)

# Supplementary Figure 9l ------------------------------------------------------

# One-way ANOVA with Tukey's HSD (P2 primer pair)
rstatix::tukey_hsd(aov(FC_RBP45B~Group,data=filter(qpcr,Primer=="SLP2_F"&Strand=="Sense"&Timepoint=="DS")))
rstatix::tukey_hsd(aov(FC_RBP45B~Group,data=filter(qpcr,Primer=="SLP2_F"&Strand=="Antisense"&Timepoint=="DS")))

ggplot(filter(qpcr_s, Timepoint=="DS"&Primer=="SLP2_F"), aes(Strand, Mean*100, fill = Genotype)) +
  geom_errorbar(aes(ymin = (Mean*100) / 2, ymax = (Mean + Sd)*100),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.85), name = "FC to RBP45B x 102") +
  scale_fill_manual(values = c("black", "grey25", "grey75"), name = element_blank()) +
  theme_bw() +
  theme(
    legend.text = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
ggsave("Figure7/Supp9l_qpcr_slp2_p2.pdf", width = 2.5, height = 1.2)

# Supplementary Figure 10b -----------------------------------------------------

atac1 <- resize(import(ACR_500BP_BED), 1, "center")
atacA <- readr::read_tsv(ATAC_CLUSTERS)

csrna_p <- import(CSRNA_ALL_PLUS_BG)
csrna_m <- import(CSRNA_ALL_MINUS_BG)
strand(csrna_p) <- "+"
strand(csrna_m) <- "-"
csrna <- sort(c(csrna_p, csrna_m))
rm(csrna_p, csrna_m)

# Distal intergenic ACRs
atac1inter <- atac1[atac1$name %in% atacA$PeakID[atacA$Type=="Intergenic"]]
atac1inter$Dist <- structure(atacA$NearestTSSDist,names=atacA$PeakID)[atac1inter$name]
strand(atac1inter) <- "+"

atacInter_csrna <- get_stranded_signal(csrna, resize(atac1inter[atac1inter$Dist>=1000],2000,"center"), row.norm=F, keep=c(0,.995))

pdf("Figure7/Supp10b_csrna_intergenic_acrs.pdf", width = 4, height = 3.2)
plot(seq(-1000,999), atacInter_csrna$PosAvg,
  ylim = c(0, max(c(atacInter_csrna$PosAvg,atacInter_csrna$NegAvg))*1.1),
  type = "l", lwd = 2, xlab = "Distance from ACR midpoint (bp)", ylab = "csRNA-seq (RPM)",
  bty = "n", main = "Distal intergenic ACRs (n = 3,093)", col = "blue")
lines(seq(-1000, 999), atacInter_csrna$NegAvg, col = "forestgreen", lwd = 2)
legend("topright", bty = "n", col = c("blue", "forestgreen"), lty = 1, lwd = 2,
  legend = c("Sense", "Antisense"))
dev.off()

# Supplementary Figure 10d -----------------------------------------------------

luc <- readxl::read_xlsx(LUC_ASSAY)
summary(luc$FoldChange[luc$Construct=="SPT"])

ggplot(filter(luc, Construct=="SPT"), aes(Construct, FoldChange)) +
  geom_jitter(width = 0.12,col="grey50") +
  geom_errorbar(aes(ymin=2.47,ymax=4.62), width = 0.2) +
  geom_errorbar(aes(ymin=3.19,ymax=3.19),width=.3) +
  scale_x_discrete(name = element_blank(),labels="eSPT") +
  scale_y_continuous(limits = c(0, 8.5), name = "Fold change negative control") +
  geom_hline(yintercept=1,col="black",lty=2) +
  theme_bw() +
  theme(
    legend.text = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
ggsave("Figure7/Supp10d_luc_spt.pdf",width=1.02,height=1.62)

# Supplementary Figure 10f -----------------------------------------------------

table(ACR=!is.na(enhA$ACR),Bidir=enhA$BiDir)

# Supplementary Figure 10g -----------------------------------------------------

table(enhA$MaxActivitiy)

ggplot(data.frame(Timepoint=factor(colnames(enhActivity),levels=colnames(enhActivity)),
    Count=c(133,261,285,216,450,546)),
  aes(Timepoint,Count))+
  geom_bar(fill="black",stat="identity",width=.75)+
  scale_y_continuous(expand=c(0,0),limits=c(0,600),name="# of enhancers at max activity") +
  xlab(element_blank())+
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(6, "pt")),
    legend.position = "none",
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    legend.key.size = unit(6, "pt"),
    legend.text = element_text(colour = "black", size = unit(6, "pt")),
    panel.grid = element_blank(),
    panel.border = element_blank()
  )
ggsave("Figure7/Supp10g_enhancer_max_activity_counts.pdf",width=1.8,height=1.3)

# Supplementary Figure 10j -----------------------------------------------------

# pcTSSs within 5 kbp of an enhancer
tss1pc <- tss1[tss1$name %in% tss_pc]
tss1pcEnh <- tss1pc[overlapsAny(tss1pc,resize(enhAll500,10000,"center"))]

enhActivityL2 <- log2(1+enhActivity)
tssMmL2 <- log2(1+tssMm)

enhCor <- matrix(0,ncol=length(enhAll500),nrow=length(tss1pcEnh),dimnames=list(tss1pcEnh$name,enhAll500$name))
for (i in 1:nrow(enhCor)) {
  enhCor[i,] <- cor(tssMmL2[rownames(enhCor)[i],],t(enhActivityL2[colnames(enhCor),]))[1,]
}
enhCor2 <- reshape2::melt(enhCor)
enhCor2$Var1 <- as.character(enhCor2$Var1)
enhCor2$Var2 <- as.character(enhCor2$Var2)
colnames(enhCor2) <- c("PcTSS", "Enhancer", "PCC")
enhCor3 <- enhCor2[substr(enhCor2$PcTSS,4,4)!=substr(enhCor2$Enhancer,4,4),]
enhCor2 <- enhCor2[substr(enhCor2$PcTSS,4,4)==substr(enhCor2$Enhancer,4,4),]
enhCor2 <- enhCor2[abs(
    as.integer(gsub("TSS\\d_","",enhCor2$PcTSS)) - as.integer(gsub("ENH\\d_","",enhCor2$Enhancer))
  ) <= 5000, ]
enhCor2 <- enhCor2[!is.na(enhCor2$PCC), ]
enhCor3 <- enhCor3[!is.na(enhCor3$PCC), ]

pdf("Figure7/Supp10j_enhancer_nearby_pctss_pcc.pdf",width=4.5,height=3.5)
plot(density(enhCor2$PCC,bw=.1),lwd=3,xlab="PCC",ylab="Density")
lines(density(enhCor3$PCC,bw=.1),col="red",lwd=3)
legend("topleft",lwd=3,col=c("black","red"),lty=1,bty="n",
  legend=c("Within 5 kb", "Different chromosome"))
dev.off()

table(enhCor2$PCC>=.5)

# Supplementary Figure 10k -----------------------------------------------------

enhLinks <- readr::read_tsv(ENHANCER_GENE_LINKS)

enhGO <- gprofiler2::gost(unique(enhLinks$Gene), "athaliana")$result
enhGOtop <- enhGO[order(enhGO$p_value),]
enhGOtop <- enhGOtop[1:10,]

ggplot(enhGOtop,aes(reorder(term_name,p_value),-log10(p_value),fill=source))+
  geom_bar(stat="identity",width=.75) +
  scale_y_continuous(expand=c(0,0),limits=c(0,14),name="-log10(P-value)") +
  scale_fill_manual(values=c("black","grey"),name=element_blank()) +
  xlab(element_blank()) +
  theme_bw() +
  theme(
    axis.text.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.x = element_text(colour = "black", size = unit(6, "pt")),
    axis.line.x.bottom = element_line(colour = "black", size = .5),
    axis.line.y.left = element_line(colour = "black", size = .5),
    legend.key.size = unit(6, "pt"),
    legend.text = element_text(colour = "black", size = unit(6, "pt")),
    panel.grid = element_blank(),
    panel.border = element_blank()
  ) +
  coord_flip()
ggsave("Figure7/Supp10k_enhancer_genes_go.pdf",width=4.44,height=1.72)
