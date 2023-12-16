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

Antisenses$PcTSS <- gene2tss[Antisenses$pcGene]
Antisenses <- Antisenses[order(rowSums(tssM[Antisenses$asTSS, ]), decreasing = TRUE), ]

Antisenses$PCC <- mapply(function(x, y) cor(tssM[x, ], tssM[y, ]),
  Antisenses$asTSS, Antisenses$pcTSS)
Antisenses <- Antisenses[!is.na(Antisenses$PCC), ]

A11 <- rtracklayer::import("Araport11_novel.gff3")
gene2symbol <- structure(A11$symbol, names = A11$ID)
gene2symbol[is.na(gene2symbol)] <- names(gene2symbol)[is.na(gene2symbol)]

Antisenses$pcSymbol <- gene2symbol[Antisenses$pcGene]
Antisenses$asGene <- structure(TSS_anno$NearestGene, names = TSS_anno$TSS)[Antisenses$asTSS]

Genes <- genes(A11)
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

Antisenses$Dist_ratio <- with(Antisenses, pmin(1, SenseSize - TSS_dist) / SenseSize)

# Proximal TSS: with first half of gene and at least 1 kbp away from gene TTS

Antisenses$Proximal <- Antisenses$Dist_ratio <= 0.5 & Antisenses$TTS_dist > 1000

# Save all antisense TSSs

readr::write_tsv(Antisenses, "All_antisenses.tsv")

# For the manuscript, we only analyzed one antisense TSS per protein coding gene
# to prevent duplicate values from influencing distributions. For genes with multiple
# antisenses, the highest expressed one across all samples was kept.

Antisenses <- Antisenses[!duplicated(Antisenses$pcGene), ]

readr::write_tsv(Antisenses, "Filtered_antisenses.tsv")

