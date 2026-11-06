# Benjamin J.M. Tremblay
#
# Figure 6 and Supplementary Figures 7, 8
#
# Package requirements: rtracklayer, matrixStats, ComplexHeatmap, circlize,
# viridisLite, magick, ggplot2, dplyr, tidyr, readxl, rstatix, eulerr, readr,
# zoo, genomation, BSgenome.Athaliana.TAIR.TAIR9, universalmotif

library(rtracklayer)
library(matrixStats)
library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(eulerr)
library(BSgenome.Athaliana.TAIR.TAIR9)
seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

source("Functions.R")

dir.create("Figure6", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_SIZE1_BED <- "tss.size1.bed"
BIDIR_PAIRS <- "Bidir_all_pairs.tsv"
ATAC_CLUSTERS <- "ATACseq_clustering3/Clusters_anno.tsv"
THIEFFRY2020_PROMPTS <- "Thieffry2020_prompts.bed"

# External tables
THIEFFRY2020_TABLE <- "TPC2019-LSB-00815R2_Supplemental_Data_Set_1.xlsx"
THIEFFRY2022_TABLE <- "tpc.22.00124Supplemental Dataset S8  Enhancer candidates.xlsx"
KINDGREN2020_TABLE <- "Supplementary Tables.xlsx"

# RT-qPCR
QPCR_AT1G04170 <- "Cristina_qPCR_foldchange_RBP45B.xlsx"
QPCR_AT3G26650 <- "Cristina_qPCR_foldchange_RBP45B_c2.xlsx"

# Read density
CSRNA_L57_PLUS_BW <- "L57.merged.csRNA.plus.bw"
CSRNA_L57_MINUS_BW <- "L57.merged.csRNA.minus.bw"
HEN2_RNA_PLUS_BW <- "hen2.merged.RNA.plus.bw"
HEN2_RNA_MINUS_BW <- "hen2.merged.RNA.minus.bw"
ATAC_DS_BW <- "atacseq_bw/DS.bw"
ATAC_L6_BW <- "atacseq_bw/L6.bw"
ATAC_L26_BW <- "atacseq_bw/L26.bw"
ATAC_L57_BW <- "atacseq_bw/L57.bw"
RNAPII_BW <- "external_bigWigs/RNAPII_14D.bw"
MNASE_BW <- "Ath_leaf_NPS.bw"
PHYLOP_BW <- "PhyloP.bw"
PHASTCONS_BW <- "PhastCons.bw"

# Shared data ------------------------------------------------------------------

colorBlindBlack8  <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tss1 <- import(TSS_SIZE1_BED)

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

atac <- readr::read_tsv(ATAC_CLUSTERS)
atacM <- as.matrix(atac[, 13:20])
rownames(atacM) <- atac$PeakID
atacMm <- as.matrix(data.frame(row.names = rownames(atacM),
  DS = rowMeans(atacM[, 1:2]),
  L6 = rowMeans(atacM[, 3:4]),
  L26 = rowMeans(atacM[, 5:6]),
  L57 = rowMeans(atacM[, 7:8])
))

# Correlation over the Col-0 samples
calc_cor <- function(x, y) cor(tssM[x,],tssM[y,])

bd_all <- readr::read_tsv(BIDIR_PAIRS)
bd_all$cor2 <- mapply(calc_cor, bd_all$TSS1, bd_all$TSS2)

# Divergent (ncTSS-pcTSS) and bidirectional non-coding (ncTSS-ncTSS) pairs within 1 kbp
div <- bd_all[bd_all$TypeV2=="nc-pc" & bd_all$dist <= 1000, ]
ncnc <- bd_all[bd_all$TypeV2=="nc-nc" & bd_all$dist <= 1000, ]

ncncBedMid <- with(ncnc, GRanges(seq1,
  IRanges(floor((pos2+pos1)/2), floor((pos2+pos1)/2)),
  len = pos2 - pos1,
  TSS1 = TSS1, TSS2 = TSS2,
))
ncncBedMid$len <- abs(ncncBedMid$len)
ncncBedMid <- ncncBedMid[order(ncncBedMid$len)]

# Midpoint between the two TSSs
divBedMid <- with(div, GRanges(seq1,
  IRanges(floor((pos2+pos1)/2), floor((pos2+pos1)/2)),
  strand = ifelse(Type2=="protein_coding", "+", "-"),
  len = pos2 - pos1,
  DivTSS = ifelse(Type2=="protein_coding", TSS1, TSS2),
  PcTSS = ifelse(Type2=="protein_coding", TSS2, TSS1)
))
strand(divBedMid)[divBedMid$len < 0] <- c("+" = "-", "-" = "+")[as.character(strand(divBedMid))[divBedMid$len < 0]]
divBedMid$len <- abs(divBedMid$len)
divBedMid <- divBedMid[order(divBedMid$len)]

# Centred on the pcTSS, sorted by inter-TSS distance
divBed <- with(div, GRanges(seq1,
  IRanges(ifelse(Type2=="protein_coding",pos2,pos1),ifelse(Type2=="protein_coding",pos2,pos1)),
  strand = ifelse(Type2=="protein_coding", "+", "-"),
  len = pos2 - pos1, PCC = PCC,
  DivTSS = ifelse(Type2=="protein_coding", TSS1, TSS2),
  PcTSS = ifelse(Type2=="protein_coding", TSS2, TSS1)
))
strand(divBed)[divBed$len < 0] <- c("+" = "-", "-" = "+")[as.character(strand(divBed))[divBed$len < 0]]
divBed$len <- abs(divBed$len)
divBed <- divBed[order(divBed$len)]
divBed500 <- divBed[divBed$len<=500]

# Figure 6a --------------------------------------------------------------------

# Counts for the pie chart
table(bd_all$TypeV2[bd_all$dist<=500])

ggplot(data.frame(Type=c("nc-nc","nc-pc","pc-pc"),N=c(381,1643,1127)),
  aes("",N,fill=Type)) +
  geom_bar(width=1,stat="identity")+
  scale_fill_manual(values=c("grey30","grey60","black")) +
  coord_polar("y", start=0)+
  theme_void()
ggsave("Figure6/Fig6a_pie1.pdf", width=2,height=1)

# Figure 6b --------------------------------------------------------------------

ggplot() +
  geom_smooth(data = filter(bd_all, TypeV2=="nc-pc"&dist<=4000),
    aes(dist, cor2), colour = "black", method = "loess", span = .65) +
  geom_smooth(data = filter(bd_all, TypeV2=="nc-nc"&dist<=4000),
    aes(dist, cor2), colour = "red", method = "loess", span = .65)+
  geom_smooth(data = filter(bd_all, TypeV2=="pc-pc"&dist<=4000),
    aes(dist, cor2), colour = "blue", method = "loess", span = .65) +
  theme_bw() +
  xlab("Distance between TSSs") +
  ylab("Correlation") +
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
ggsave("Figure6/Fig6b_corr.pdf", width=1.91,height=1.41)

# Figure 6c --------------------------------------------------------------------

div500 <- filter(div, dist <= 500)
div500$MaxCPM1 <- rowMaxs(tssMm[div500$TSS1, ])
div500$MaxCPM2 <- rowMaxs(tssMm[div500$TSS2, ])
div500$MaxRatio <-
  (1 + ifelse(div500$Type1=="protein_coding", div500$MaxCPM1, div500$MaxCPM2)) /
  (1 + ifelse(div500$Type1=="protein_coding", div500$MaxCPM2, div500$MaxCPM1))

div500$CorCat <- "Neutral"
div500$CorCat[div500$PCC >= .25] <- "Positive"
div500$CorCat[div500$PCC <= -.25] <- "Negative"
div500$Log2MaxRatio <- log2(div500$MaxRatio)

# Two-sided Mann-Whitney tests
rstatix::wilcox_test(div500, MaxRatio~CorCat)

# Medians used for the dashed lines
quantile(div500$MaxRatio[div500$PCC>=.25], p = c(.25, .5, .75))
quantile(div500$MaxRatio[div500$PCC<.25 & div500$PCC> -.25], p = c(.25, .5, .75))
quantile(div500$MaxRatio[div500$PCC<= -.25], p = c(.25, .5, .75))

ggplot(div500, aes(log2(MaxRatio), colour = CorCat)) +
  ylab("Density") +
  geom_vline(xintercept = log2(c(7.065)), colour = "red", lty = 2) +
  geom_vline(xintercept = log2(c(8.534)), colour = "black", lty = 2) +
  geom_vline(xintercept = log2(c(12.290)), colour = "blue", lty = 2) +
  geom_density(bw = 1) +
  coord_cartesian(xlim = c(-8, 11)) +
  scale_x_continuous(name = "Ratio of max coding/non-coding expression",
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
ggsave("Figure6/Fig6c_coding_noncoding_ratio.pdf", height = 1.33, width = 2)

# Figure 6d --------------------------------------------------------------------

# csRNA-seq (L57)

csrna_57 <- read_signal(
  pos = CSRNA_L57_PLUS_BW,
  neg = CSRNA_L57_MINUS_BW
)

csrna_div <- get_stranded_signal(
  csrna_57,
  windows = promoters(divBed, 2000, 1000),
  row.norm=T,
  return.matrix = T,
  merge = T
)

pdf("Figure6/Fig6d_div_csrna.pdf", width = 2, height = 1.72)
Heatmap(t(apply(csrna_div, 1, nn0)),
  show_row_names = FALSE, show_column_names = FALSE,
  use_raster=T,raster_quality=4,
  col = circlize::colorRamp2(c(-.1, 0, .1),
      c("red", "white", "blue")),
  raster_by_magick = TRUE,
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# RNA-seq (hen2-4)

hen2 <- read_signal(
  pos = HEN2_RNA_PLUS_BW,
  neg = HEN2_RNA_MINUS_BW)

hen2_div <- get_stranded_signal(hen2, keep = c(0., .90),
  windows = promoters(divBed, 2000, 1000), return.matrix = TRUE, merge = TRUE)

pdf("Figure6/Fig6d_div_rna_hen2.pdf", width = 2, height = 1.72)
Heatmap(hen2_div,
  show_row_names = FALSE, show_column_names = FALSE,
  use_raster=T,raster_quality=4,
  col = circlize::colorRamp2(c(-.1, 0, .1),
      c("red", "white", "blue")),
  raster_by_magick = TRUE,
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# RNAPII ChIP-seq

rnap <- import(RNAPII_BW)

rnap_div <- get_unstranded_signal(rnap, promoters(divBed, 2000, 1000), keep=c(0.,.9),
  return.matrix=TRUE)

pdf("Figure6/Fig6d_div_rnapii.pdf", width = 2, height = 1.72)
Heatmap(rnap_div,
  show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
  col=viridisLite::plasma(100,begin=0),
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# MNase-seq

mnase <- import(MNASE_BW)
seqlevels(mnase) <- as.character(1:5)
mnase_div <- get_unstranded_signal(mnase, keep = c(0., .90),
  windows = promoters(divBed, 2000, 1000), return.matrix = TRUE)

pdf("Figure6/Fig6d_div_mnase.pdf", width = 2, height = 1.72)
Heatmap(mnase_div,
  show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
  col=viridisLite::plasma(100,begin=0),
  cluster_rows = FALSE, cluster_columns= FALSE
)
dev.off()

# Figure 6g --------------------------------------------------------------------

qpcr <- readxl::read_xlsx(QPCR_AT1G04170)
qpcr <- tidyr::gather(qpcr, "Sample", "FC_RBP45B", 5:ncol(qpcr))
qpcr$Timepoint <- gsub("[.][123]$", "", qpcr$Sample)
qpcr_s <- qpcr %>%
  group_by(Plate, Genotype, Gene, RNA, Timepoint) %>%
  summarise(Mean = mean(FC_RBP45B), Sd = sd(FC_RBP45B))

qpcr$Group <- paste0(qpcr$Genotype, "_", qpcr$Timepoint)

# Bonferroni correction
adjp <- function(x) { x$p.adj2 <- x$p * 3 ; x }

# Col-0 vs SALK_201027C
adjp(rstatix::t_test(filter(qpcr, Plate==3 & Gene=="PcGene"), FC_RBP45B~Group))[, c(2,3,8,11)]
adjp(rstatix::t_test(filter(qpcr, Plate==3 & Gene=="NcRNA"), FC_RBP45B~Group))[, c(2,3,8,11)]

ggplot(filter(qpcr_s, RNA=="Stable" & Plate == 3 & Gene=="PcGene"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "S72", "L168"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[1:2], name = element_blank()) +
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
ggsave("Figure6/Fig6g_qpcr_pcgene.pdf", width = 2.6, height = 1.01)

ggplot(filter(qpcr_s, RNA=="Stable" & Plate == 3 & Gene=="NcRNA"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "S72", "L168"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.005), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[1:2], name = element_blank()) +
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
ggsave("Figure6/Fig6g_qpcr_ncrna.pdf", width = 2.6, height = 1.01)

# Figure 6h --------------------------------------------------------------------

# Col-0 vs SALK_073206
adjp(rstatix::t_test(filter(qpcr, Plate==4 & Gene=="PcGene"), FC_RBP45B~Group))[, c(2,3,8,11)]
adjp(rstatix::t_test(filter(qpcr, Plate==4 & Gene=="NcRNA"), FC_RBP45B~Group))[, c(2,3,8,11)]

ggplot(filter(qpcr_s, RNA=="Stable" & Plate == 4 & Gene=="PcGene"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "S72", "L168"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 2), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[c(1,3)], name = element_blank()) +
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
ggsave("Figure6/Fig6h_qpcr_pcgene.pdf", width = 2.6, height = 1.01)

ggplot(filter(qpcr_s, RNA=="Stable" & Plate == 4 & Gene=="NcRNA"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "S72", "L168"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.006), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[c(1,3)], name = element_blank()) +
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
ggsave("Figure6/Fig6h_qpcr_ncrna.pdf", width = 2.6, height = 1.01)

# Supplementary Figure 7a ------------------------------------------------------

t20_div <- import(THIEFFRY2020_PROMPTS)

k20 <- readxl::read_xlsx(KINDGREN2020_TABLE, sheet = "Table S1 Novel transcripts", skip = 2)
k20_tx <- GRanges(as.character(k20$Chromosome), IRanges(k20$start, k20$end),
  strand = k20$strand)

# Full extent of the divergent promoters
dStart <- dEnd <- start(divBed)
dStart[as.character(strand(divBed)) == "+"] <- dStart[as.character(strand(divBed)) == "+"] -
  divBed$len[as.character(strand(divBedMid)) == "+"]
dEnd[as.character(strand(divBed)) == "-"] <- dEnd[as.character(strand(divBed)) == "-"] +
  divBed$len[as.character(strand(divBed)) == "-"]
divBedFull <- divBed
start(divBedFull) <- pmin(dStart, dEnd)
end(divBedFull) <- pmax(dStart, dEnd)

ncncBedFull <- resize(ncncBedMid, ncncBedMid$len, "center")

divBedFull <- sort(divBedFull[divBedFull$len <= 500])
ncncBedFull <- sort(ncncBedFull[ncncBedFull$len <= 500])

table(overlapsAny(divBedFull, promoters(t20_div, 400, 100), ignore.strand=T))
table(overlapsAny(divBedFull, promoters(k20_tx, 400, 100), ignore.strand=T))

divBedFull$T2020 <- overlapsAny(divBedFull, promoters(t20_div, 400, 100), ignore.strand=T)
divBedFull$K2020 <- overlapsAny(divBedFull, promoters(k20_tx, 400, 100), ignore.strand=T)

t20_k20_div <- t20_div[overlapsAny(promoters(t20_div, 400, 100),promoters(k20_tx,400,100),ignore.strand=T)]

# Total = 2141
pdf("Figure6/Supp7a_in_common_div.pdf", width = 1, height = 1)
plot(euler(c(
  "Tremblay&Thieffry2020&Kindgren2020" = sum(divBedFull$T2020 & divBedFull$K2020),
  "Thieffry2020&Kindgren2020" = length(t20_k20_div[!overlapsAny(promoters(t20_k20_div,400,100), divBedFull, ignore.strand=T)]),
  "Tremblay&Thieffry2020" = sum(divBedFull$T2020 & !divBedFull$K2020),
  "Tremblay&Kindgren2020" = sum(!divBedFull$T2020 & divBedFull$K2020),
  "Tremblay" = sum(!divBedFull$T2020 & !divBedFull$K2020),
  "Thieffry2020" = sum(!overlapsAny(promoters(t20_div,400,100)[!overlapsAny(promoters(t20_div,400,100),promoters(k20_tx,400,100),ignore.strand=T)],divBedFull,ignore.strand=T)),
  "Kindgren2020" = 917 - 35 - 20 - 401
)), edges = FALSE, fills = list(fill = c("grey25", "grey50", "grey75"))
)
dev.off()

# Supplementary Figure 7b ------------------------------------------------------

t20 <- readxl::read_xlsx(THIEFFRY2020_TABLE, sheet = "Bidirectional TCs")
t20_ncnc <- GRanges(sub("Chr", "", t20$chr), IRanges(t20$start, t20$end))

t22 <- readxl::read_xlsx(THIEFFRY2022_TABLE, sheet = "enhancer_candidates", skip = 1)
t22_ncnc <- GRanges(sub("Chr", "", t22$chr), IRanges(t22$start, t22$end))

table(overlapsAny(resize(ncncBedFull, 500, "center"), t20_ncnc, ignore.strand = TRUE))
table(overlapsAny(resize(ncncBedFull, 500, "center"), t22_ncnc, ignore.strand = TRUE))
table(overlapsAny(t20_ncnc, t22_ncnc, ignore.strand = TRUE))

t20_22_ncnc <- t20_ncnc[overlapsAny(t20_ncnc, t22_ncnc, ignore.strand=T)]

ncncBedFull$T2020 <- overlapsAny(resize(ncncBedFull, 500, "center"), t20_ncnc, ignore.strand = TRUE)
ncncBedFull$T2022 <- overlapsAny(resize(ncncBedFull, 500, "center"), t22_ncnc, ignore.strand = TRUE)

# Total = 581
pdf("Figure6/Supp7b_in_common_ncnc.pdf", width = 1, height = 1)
plot(euler(c(
  "Tremblay&Thieffry2020&Thieffry2022" = sum(ncncBedFull$T2020 & ncncBedFull$T2022),
  "Thieffry2020&Thieffry2022" = length(t20_22_ncnc[!overlapsAny(t20_22_ncnc, resize(ncncBedFull, 500, "center"), ignore.strand=T)]),
  "Tremblay&Thieffry2020" = sum(ncncBedFull$T2020 & !ncncBedFull$T2022),
  "Tremblay&Thieffry2022" = sum(!ncncBedFull$T2020 & ncncBedFull$T2022),
  "Tremblay" = sum(!ncncBedFull$T2020 & !ncncBedFull$T2022),
  "Thieffry2020" = sum(!overlapsAny(t20_ncnc[!overlapsAny(t20_ncnc,t22_ncnc,ignore.strand=T)],resize(ncncBedFull, 500, "center"),ignore.strand=T)),
  "Thieffry2022" = sum(!overlapsAny(t22_ncnc[!overlapsAny(t22_ncnc,t20_ncnc,ignore.strand=T)],resize(ncncBedFull, 500, "center"),ignore.strand=T))
)), edges = FALSE, fills = list(fill = c("grey25", "grey50", "grey75"))
)
dev.off()

# Supplementary Figure 7c ------------------------------------------------------

pdf("Figure6/Supp7c_div_size_density.pdf", width = 4, height = 3.5)
plot(density(divBed500$len, bw=25), lwd = 3, xlim = c(0, 500),
  ylim = c(0, 0.0035), main = NA, xlab = "Distance from pcTSS (bp)",
  ylab = "Density"
)
abline(v=206.7, lty = 2, lwd = 2, col = "grey50")
dev.off()

# Supplementary Figure 7d ------------------------------------------------------

div_atac_ds <- get_unstranded_signal(ATAC_DS_BW, promoters(divBed500,1000,0))
div_atac_l6 <- get_unstranded_signal(ATAC_L6_BW, promoters(divBed500,1000,0))
div_atac_l26 <- get_unstranded_signal(ATAC_L26_BW, promoters(divBed500,1000,0))
div_atac_l57 <- get_unstranded_signal(ATAC_L57_BW, promoters(divBed500,1000,0))

which.max(div_atac_ds$Avg)
which.max(div_atac_l6$Avg)
which.max(div_atac_l26$Avg)
which.max(div_atac_l57$Avg)

pdf("Figure6/Supp7d_div_atac_prob.pdf", width = 4, height = 3.5)
plot(-div_atac_ds$Coord[501:1000],div_atac_ds$Avg[501:1000],type="l",lwd=3,
  xaxt = "n", xlab = "Distance from pcTSS (bp)", ylab = "ATAC-seq signal probability",
  col=colorBlindBlack8[1], ylim=0:1)
axis(1,-c(500,600,700,800,900,1000),c(-500,-400,-300,-200,-100,0))
lines(-div_atac_l6$Coord[501:1000], div_atac_l6$Avg[501:1000],lwd=3,col=colorBlindBlack8[4])
lines(-div_atac_l26$Coord[501:1000], div_atac_l26$Avg[501:1000],lwd=3,col=colorBlindBlack8[5])
lines(-div_atac_l57$Coord[501:1000], div_atac_l57$Avg[501:1000],lwd=3,col=colorBlindBlack8[6])
abline(v=-894, lty=2,lwd=2,col="black")
dev.off()

# Supplementary Figure 7e ------------------------------------------------------

# Unidirectional pcTSSs
tss1_pc_uni <- tss1[tss1$name %in% tss_pc]
tss1_pc_uni <- tss1_pc_uni[!tss1_pc_uni$name %in% divBed$PcTSS]

uni_atac_ds <- get_unstranded_signal(ATAC_DS_BW, promoters(tss1_pc_uni,1000,0))
uni_atac_l6 <- get_unstranded_signal(ATAC_L6_BW, promoters(tss1_pc_uni,1000,0))
uni_atac_l26 <- get_unstranded_signal(ATAC_L26_BW, promoters(tss1_pc_uni,1000,0))
uni_atac_l57 <- get_unstranded_signal(ATAC_L57_BW, promoters(tss1_pc_uni,1000,0))

which.max(uni_atac_ds$Avg)
which.max(uni_atac_l6$Avg)
which.max(uni_atac_l26$Avg)
which.max(uni_atac_l57$Avg)

pdf("Figure6/Supp7e_uni_atac_prob.pdf", width = 4, height = 3.5)
plot(-uni_atac_ds$Coord[501:1000],uni_atac_ds$Avg[501:1000],type="l",lwd=3,
  xaxt = "n", xlab = "Distance from pcTSS (bp)", ylab = "ATAC-seq signal probability",
  col=colorBlindBlack8[1], ylim=0:1)
axis(1,-c(500,600,700,800,900,1000),c(-500,-400,-300,-200,-100,0))
lines(-uni_atac_l6$Coord[501:1000], uni_atac_l6$Avg[501:1000],lwd=3,col=colorBlindBlack8[4])
lines(-uni_atac_l26$Coord[501:1000], uni_atac_l26$Avg[501:1000],lwd=3,col=colorBlindBlack8[5])
lines(-uni_atac_l57$Coord[501:1000], uni_atac_l57$Avg[501:1000],lwd=3,col=colorBlindBlack8[6])
abline(v=-950, lty=2,lwd=2,col="black")
dev.off()

# Supplementary Figure 7f ------------------------------------------------------

div_phylop <- get_unstranded_signal(PHYLOP_BW, promoters(divBed500,1000,0), row.norm=F)
uni_phylop <- get_unstranded_signal(PHYLOP_BW, promoters(tss1_pc_uni,1000,0), row.norm=F)

pdf("Figure6/Supp7f_phylop.pdf", width = 4, height = 3.5)
plot(-uni_phylop$Coord[501:1000],uni_phylop$Avg[501:1000],type="l",lwd=3,
  xaxt = "n", xlab = "Distance from pcTSS (bp)", ylab = "PhyloP score",
  col="black")
axis(1,-c(500,600,700,800,900,1000),c(-500,-400,-300,-200,-100,0))
lines(-div_phylop$Coord[501:1000], div_phylop$Avg[501:1000],lwd=3,col="grey75")
dev.off()

# Supplementary Figure 7g ------------------------------------------------------

div_phast <- get_unstranded_signal(PHASTCONS_BW, promoters(divBed500,1000,0), row.norm=F)
uni_phast <- get_unstranded_signal(PHASTCONS_BW, promoters(tss1_pc_uni,1000,0), row.norm=F)

pdf("Figure6/Supp7g_phastcons.pdf", width = 4, height = 3.5)
plot(-uni_phast$Coord[501:1000],uni_phast$Avg[501:1000],type="l",lwd=3,
  xaxt = "n", xlab = "Distance from pcTSS (bp)", ylab = "PhastCons score",
  col="black")
axis(1,-c(500,600,700,800,900,1000),c(-500,-400,-300,-200,-100,0))
lines(-div_phast$Coord[501:1000], div_phast$Avg[501:1000],lwd=3,col="grey75")
dev.off()

# Supplementary Figure 7h ------------------------------------------------------

get_freqs <- function(x) {
  list(
    A = apply(x, 2, function(x) sum(x=="A")) / nrow(x),
    C = apply(x, 2, function(x) sum(x=="C")) / nrow(x),
    G = apply(x, 2, function(x) sum(x=="G")) / nrow(x),
    T = apply(x, 2, function(x) sum(x=="T")) / nrow(x)
  )
}

bidir <- readr::read_tsv(BIDIR_PAIRS)
div2 <- bidir[bidir$dist<=500 & bidir$TypeV2 == "nc-pc", ]
div_pc <- ifelse(div2$Type1=="protein_coding",div2$TSS1,div2$TSS2)
div_nc <- ifelse(div2$Type1=="protein_coding",div2$TSS2,div2$TSS1)

tss1fadiv <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, promoters(tss1[tss1$name %in% div_pc], 1000, 500))
names(tss1fadiv) <- tss1$name[tss1$name %in% div_pc]
tss1fadivM <- as.matrix(t(as.data.frame(strsplit(as.character(tss1fadiv), "", fixed=T))))

f_div <- get_freqs(tss1fadivM[div_pc,,drop=F])

pdf("Figure6/Supp7h_div_pctss_base_freq.pdf", width = 3.8, height = 3.5)
plot(c((-999):0, 1:500), f_div$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-1000, -500, 0, 500))
axis(2, c(0, .5, 1))
lines(c((-999):0, 1:500), f_div$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-999):0, 1:500), f_div$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-999):0, 1:500), f_div$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 7i ------------------------------------------------------

tss1fadiv2 <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, promoters(tss1[tss1$name %in% div_nc], 1000, 500))
names(tss1fadiv2) <- tss1$name[tss1$name %in% div_nc]
tss1fadiv2M <- as.matrix(t(as.data.frame(strsplit(as.character(tss1fadiv2), "", fixed=T))))

f_div3 <- get_freqs(tss1fadiv2M[div_nc,,drop=F])

pdf("Figure6/Supp7i_div_nctss_base_freq.pdf", width = 3.8, height = 3.5)
plot(c((-999):0, 1:500), f_div3$A, type = "l", lwd = 2, ylim = 0:1, xaxt = "n",
  ylab = "Base frequency", xlab = "Position relative to TSS (bp)", yaxt = "n",
  col = universalmotif:::DNA_COLOURS["A"], bty = "n")
axis(1, c(-1000, -500, 0, 500))
axis(2, c(0, .5, 1))
lines(c((-999):0, 1:500), f_div3$C, lwd = 2,
  col = universalmotif:::DNA_COLOURS["C"])
lines(c((-999):0, 1:500), f_div3$G, lwd = 2,
  col = universalmotif:::DNA_COLOURS["G"])
lines(c((-999):0, 1:500), f_div3$T, lwd = 2,
  col = universalmotif:::DNA_COLOURS["T"])
legend("topright", lwd = 2, lty = 1, col = universalmotif:::DNA_COLOURS,
  legend = c("A", "C", "G", "T"), bty="n")
dev.off()

# Supplementary Figure 8b ------------------------------------------------------

xa <- c(1.25, 1.75, 3.25, 3.75, 5.25, 5.75, 7.25, 7.75, 9.25, 9.75, 11.25, 11.75)
xb <- c(1.5, 3.5, 5.5, 7.5, 9.5, 11.5)

# LHY pcTSS
pdf("Figure6/Supp8b_lhy_pc_csrna.pdf", width = 3.9, height = 3.5)
plot(xa,
  tssM["TSS1_37853", ], ylim = c(0, 70),
  pch = 19, bty = "n", xaxt = "n", xlab = NA, ylab = "csRNA-seq (CPM)"
)
lines(xb, tssMm["TSS1_37853", ])
axis(1, xb, c("DS", "S24", "S72", "L6", "L26", "L57"))
dev.off()

# LHY divergent ncTSS
pdf("Figure6/Supp8b_lhy_div_csrna.pdf", width = 3.9, height = 3.5)
plot(xa,
  tssM["TSS1_38250", ], ylim = c(0, 5),
  pch = 19, bty = "n", xaxt = "n", xlab = NA, ylab = "csRNA-seq (CPM)"
)
lines(xb, tssMm["TSS1_38250", ])
axis(1, xb, c("DS", "S24", "S72", "L6", "L26", "L57"))
dev.off()

# LHY promoter ACR
pdf("Figure6/Supp8b_lhy_atac.pdf", width = 3.9, height = 3.5)
plot(xa[1:8],
  atacM["ATAC1_38105", ], ylim = c(0, 200),
  pch = 19, bty = "n", xaxt = "n", xlab = NA, ylab = "ATAC-seq (RPM)"
)
lines(xb[1:4], atacMm["ATAC1_38105", ])
axis(1, xb[1:4], c("DS", "L6", "L26", "L57"))
dev.off()

# Supplementary Figure 8e ------------------------------------------------------

qpcr2 <- readxl::read_xlsx(QPCR_AT3G26650)
qpcr2 <- tidyr::gather(qpcr2, "Sample", "FC_RBP45B", 4:ncol(qpcr2))
qpcr2$Timepoint <- gsub("[.][123]$", "", qpcr2$Sample)
qpcr2_s <- qpcr2 %>%
  group_by(Plate, Genotype, Gene, Timepoint) %>%
  summarise(Mean = mean(FC_RBP45B), Sd = sd(FC_RBP45B))

qpcr2$Group <- paste0(gsub("_", "", gsub("-", "", qpcr2$Genotype, fixed=T), fixed=T), "_", gsub("_", "", qpcr2$Timepoint,fixed=T))

# Col-0 vs SAIL_1250_D04
adjp(rstatix::t_test(filter(qpcr2, Plate==1 & Gene=="PcGene"), FC_RBP45B~Group))[, c(2,3,8,11)]
adjp(rstatix::t_test(filter(qpcr2, Plate==1 & Gene=="NcRNA"), FC_RBP45B~Group))[, c(2,3,8,11)]

ggplot(filter(qpcr2_s, Plate == 1 & Gene=="PcGene"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "L12", "L48"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[1:2], name = element_blank()) +
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
ggsave("Figure6/Supp8e_qpcr_pcgene.pdf", width = 2.6, height = 1.01)

ggplot(filter(qpcr2_s, Plate==1 & Gene=="NcRNA"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "L12", "L48"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.025), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[1:2], name = element_blank()) +
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
ggsave("Figure6/Supp8e_qpcr_ncrna.pdf", width = 2.6, height = 1.01)

# Supplementary Figure 8f ------------------------------------------------------

# Col-0 vs SALK_138567
adjp(rstatix::t_test(filter(qpcr2, Plate==2 & Gene=="PcGene"), FC_RBP45B~Group))[, c(2,3,8,11)]
adjp(rstatix::t_test(filter(qpcr2, Plate==2 & Gene=="NcRNA"), FC_RBP45B~Group))[, c(2,3,8,11)]

ggplot(filter(qpcr2_s, Plate == 2 & Gene=="PcGene"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "L12", "L48"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 25), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[c(1,3)], name = element_blank()) +
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
ggsave("Figure6/Supp8f_qpcr_pcgene.pdf", width = 2.6, height = 1.01)

ggplot(filter(qpcr2_s, Plate == 2 & Gene=="NcRNA"),
  aes(Timepoint, Mean, fill = Genotype)) +
  geom_errorbar(aes(ymin = Mean / 2, ymax = Mean + Sd),
    position = position_dodge(0.47), width = 0.3) +
  geom_bar(stat = "identity", position = "dodge", width = 0.5) +
  scale_x_discrete(limits = c("DS", "L12", "L48"), name = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.012), name = "FC to RBP45B") +
  scale_fill_manual(values = colorBlindBlack8[c(1,3)], name = element_blank()) +
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
ggsave("Figure6/Supp8f_qpcr_ncrna.pdf", width = 2.6, height = 1.01)
