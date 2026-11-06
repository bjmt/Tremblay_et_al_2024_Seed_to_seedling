# Benjamin J.M. Tremblay, June 2023
#
# Antisense TSSs

library(rtracklayer)
library(GenomicFeatures)

# Load in TSS data

TSSs1 <- import("tss.size1.bed")
TSS_anno <- readr::read_tsv("TSS.annotations.tsv")
DetTx <- import("detected.transcripts.bed")

TSS_anno$NearestGeneType[is.na(TSS_anno$NearestGeneType)] <- "none"

TSS_anno_pc <- TSS_anno[TSS_anno$NearestGeneType == "protein_coding", ]
TSS_anno_nc <- TSS_anno[TSS_anno$NearestGeneType != "protein_coding", ]

PcTx <- DetTx[DetTx$name %in% TSS_anno_pc$NearestTranscript]
TSSs1_nc <- TSSs1[TSSs1$name %in% TSS_anno_nc$TSS]

PcTx_rc <- PcTx
strand(PcTx_rc) <- c("+" = "-", "-" = "+")[as.character(strand(PcTx))]

# Find non-coding TSSs which overlap the antisense strand of a detected protein
# coding gene, or 200 bp past the TTS

Nc_ov <- findOverlaps(TSSs1_nc, c(promoters(PcTx_rc, 200, 0), PcTx_rc))

Antisenses <- data.frame(row.names = NULL,
  asTSS = TSSs1_nc$name[queryHits(Nc_ov)],
  pcGene = gsub("[.]\\d+", "", rep(PcTx_rc$name, 2)[subjectHits(Nc_ov)])
)
Antisenses <- Antisenses[!duplicated(Antisenses), ]

# Now start annotating

tssM <- as.matrix(TSS_anno[, 10:21])
rownames(tssM) <- TSS_anno$TSS
tssM <- tssM[order(rowSums(tssM), decreasing = TRUE), ]

tss2gene <- structure(TSS_anno$NearestGene, names = TSS_anno$TSS)
tss2gene <- tss2gene[rownames(tssM)]
gene2tss <- structure(names(tss2gene)[!is.na(tss2gene)], names = unname(tss2gene)[!is.na(tss2gene)])
gene2tss <- gene2tss[!duplicated(names(gene2tss))]

Antisenses$pcTSS <- gene2tss[Antisenses$pcGene]
Antisenses <- Antisenses[order(rowSums(tssM[Antisenses$asTSS, ]), decreasing = TRUE), ]

Antisenses$PCC <- mapply(function(x, y) cor(tssM[x, ], tssM[y, ]),
  Antisenses$asTSS, Antisenses$pcTSS)
Antisenses <- Antisenses[!is.na(Antisenses$PCC), ]

A11 <- rtracklayer::import("Araport11_novel.gff3")
gene2symbol <- structure(A11$symbol, names = A11$ID)
gene2symbol[is.na(gene2symbol)] <- names(gene2symbol)[is.na(gene2symbol)]

Antisenses$pcSymbol <- gene2symbol[Antisenses$pcGene]
Antisenses$asGene <- structure(TSS_anno$NearestGene, names = TSS_anno$TSS)[Antisenses$asTSS]

Genes <- genes(makeTxDbFromGFF("Araport11_novel.gff3"))
gene2size <- structure(width(Genes), names = Genes$gene_id)

Antisenses$SenseSize <- gene2size[Antisenses$pcGene]
Antisenses$TSS_dist <- NA_integer_

pair_dist <- function(x, y) {
  distance(TSSs1[TSSs1$name==x], TSSs1[TSSs1$name==y], ignore.strand=TRUE)
}

for (i in 1:nrow(Antisenses)) {
  Antisenses$TSS_dist[i] <- pair_dist(Antisenses$pcTSS[i], Antisenses$asTSS[i])
}
Antisenses$TTS_dist <- Antisenses$SenseSize - Antisenses$TSS_dist

Antisenses$Dist_ratio <- with(Antisenses, pmin(1, (SenseSize - TTS_dist) / SenseSize))

# Proximal TSS: with first half of gene and at least 1 kbp away from gene TTS

Antisenses$Proximal <- Antisenses$Dist_ratio <= 0.5 & Antisenses$TTS_dist > 1000

# Save all antisense TSSs

readr::write_tsv(Antisenses, "All_antisenses.tsv")

# For the manuscript, we only analyzed one antisense TSS per protein coding gene
# to prevent duplicate values from influencing distributions. For genes with multiple
# antisenses, the highest expressed one across all samples was kept.

Antisenses <- Antisenses[!duplicated(Antisenses$pcGene), ]

readr::write_tsv(Antisenses, "Filtered_antisenses.tsv")

# Antisense tables used in the figures -----------------------------------------

FigPairs <- data.frame(row.names = NULL,
  asTSS = TSSs1_nc$name[queryHits(Nc_ov)],
  pcGene = gsub("[.]\\d+", "", rep(PcTx_rc$name, 2)[subjectHits(Nc_ov)])
)
FigPairs <- FigPairs[!duplicated(FigPairs), ]
FigPairs$pcTSS <- gene2tss[FigPairs$pcGene]
FigPairs <- FigPairs[order(rowSums(tssM[FigPairs$asTSS, ]), decreasing = TRUE), ]

readr::write_tsv(FigPairs, "all_antisense_pairs.txt")

FigPairs <- FigPairs[!duplicated(FigPairs$pcGene), ]
FigPairs$PCC <- mapply(function(x, y) cor(tssM[x, ], tssM[y, ]),
  FigPairs$asTSS, FigPairs$pcTSS)
FigPairs <- FigPairs[!is.na(FigPairs$PCC), ]

# Mean CPM per time point

tssMm <- cbind(
  DS = rowMeans(tssM[, 1:2]),
  S24 = rowMeans(tssM[, 3:4]),
  S72 = rowMeans(tssM[, 5:6]),
  L6 = rowMeans(tssM[, 7:8]),
  L26 = rowMeans(tssM[, 9:10]),
  L57 = rowMeans(tssM[, 11:12])
)

FigAntisenses <- data.frame(row.names = NULL,
  PcTSS = FigPairs$pcTSS,
  PcGene = FigPairs$pcGene,
  PcSymbol = gene2symbol[FigPairs$pcGene],
  AntisenseTSS = FigPairs$asTSS,
  AntisenseGene = structure(TSS_anno$NearestGene, names = TSS_anno$TSS)[FigPairs$asTSS],
  Correlation = FigPairs$PCC,
  AveragePcAntisenseRatio = NA_real_,
  MaxPcCPM = matrixStats::rowMaxs(tssMm[FigPairs$pcTSS, ], na.rm = TRUE),
  MinPcCPM = matrixStats::rowMins(tssMm[FigPairs$pcTSS, ], na.rm = TRUE),
  MaxAntisenseCPM = matrixStats::rowMaxs(tssMm[FigPairs$asTSS, ], na.rm = TRUE),
  MinAntisenseCPM = matrixStats::rowMins(tssMm[FigPairs$asTSS, ], na.rm = TRUE)
)
FigAntisenses$AntisenseGene[is.na(FigAntisenses$AntisenseGene)] <- ""

# Mean sense/antisense ratio (samples with at least 1 CPM)

for (i in 1:nrow(FigAntisenses)) {
  i_1 <- tssM[FigAntisenses$PcTSS[i], ]
  i_2 <- tssM[FigAntisenses$AntisenseTSS[i], ]
  i_1 <- i_1[i_2 >= 1]
  i_2 <- i_2[i_2 >= 1]
  i_2 <- i_2[i_1 >= 1]
  i_1 <- i_1[i_1 >= 1]
  FigAntisenses$AveragePcAntisenseRatio[i] <- mean(i_1 / i_2, na.rm = TRUE)
}

FigAntisenses$MaxPcTimepoint <- colnames(tssMm)[apply(tssMm[FigAntisenses$PcTSS, ], 1, which.max)]
FigAntisenses$MinPcTimepoint <- colnames(tssMm)[apply(tssMm[FigAntisenses$PcTSS, ], 1, which.min)]
FigAntisenses$MaxAntisenseTimepoint <-
  colnames(tssMm)[apply(tssMm[FigAntisenses$AntisenseTSS, ], 1, which.max)]
FigAntisenses$MinAntisenseTimepoint <-
  colnames(tssMm)[apply(tssMm[FigAntisenses$AntisenseTSS, ], 1, which.min)]

for (i in colnames(tssMm)) {
  FigAntisenses[[paste0("PcCPM_", i)]] <- tssMm[FigAntisenses$PcTSS, i]
}
for (i in colnames(tssMm)) {
  FigAntisenses[[paste0("AntisenseCPM_", i)]] <- tssMm[FigAntisenses$AntisenseTSS, i]
}

readr::write_tsv(FigAntisenses, "Antisenses.txt")

