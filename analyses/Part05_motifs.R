# Benjamin J.M. Tremblay, June 2023
#
# Motif analyses

library(rtracklayer)
library(GenomicFeatures)
library(universalmotif)
library(Biostrings)
library(matrixStats)
library(BSgenome.Athaliana.TAIR.TAIR9)

seqlevels(BSgenome.Athaliana.TAIR.TAIR9) <- c(1:5, "Mt", "Pt")

# Get ACR and promoter sequences per cluster, and from everyone for background

ACRs500 <- import("peaks_final_500bp.bed")
TSSs500 <- promoters(import("tss.size1.bed"), 400, 100)

TSS_anno <- readr::read_tsv("TSS.annotations.clusters.tsv")
ACR_anno <- readr::read_tsv("ACR.annotations.clusters.tsv")

TSSs500$Cluster <- structure(TSS_anno$Cluster, names = TSS_anno$TSS)[TSSs500$name]
ACRs500$Cluster <- structure(ACR_anno$Cluster, names = ACR_anno$PeakID)[ACRs500$name]

TSS_bkg <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, TSSs500)
ACR_bkg <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, ACRs500)

writeXStringSet(TSS_bkg, "Promoters/Promoters_bkg.fa")
writeXStringSet(ACR_bkg, "ACRs/ACRs_bkg.fa")

for (i in paste0("C", 1:6)) {
  writeXStringSet(TSS_bkg[TSSs500$Cluster == i], paste0("Promoters/Promoters_", i, ".fa"))
}

for (i in paste0("A", 1:5)) {
  writeXStringSet(ACR_bkg[ACRs500$Cluster == i], paste0("ACRs/ACRs_", i, ".fa"))
}

# Find motifs with streme

system("bash Part05_motifs_helper1.sh")

# Load motifs and merge those with high overlap coefficients

Motifs <- c(
  read_meme("Promoters_streme/C1/streme.txt"),
  read_meme("Promoters_streme/C2/streme.txt"),
  read_meme("Promoters_streme/C3/streme.txt"),
  read_meme("Promoters_streme/C4/streme.txt"),
  read_meme("Promoters_streme/C5/streme.txt"),
  read_meme("Promoters_streme/C6/streme.txt"),
  read_meme("ACRs_streme/A1/streme.txt"),
  read_meme("ACRs_streme/A2/streme.txt"),
  read_meme("ACRs_streme/A3/streme.txt"),
  read_meme("ACRs_streme/A4/streme.txt"),
  read_meme("ACRs_streme/A5/streme.txt")
)

names(Motifs) <- sapply(Motifs, function(x) x["name"])
Motif_hits <- scan_sequences(Motifs, c(TSS_bkg[TSSs500$Cluster!="C0"],
  ACR_bkg[ACRs500$Cluster!="A0"]), threshold = 0.0001, RC = TRUE,
  return.granges = TRUE, no.overlaps = TRUE)

Motif_ov <- matrix(0, ncol = length(Motifs), nrow = length(Motifs),
  dimnames = list(names(Motifs), names(Motifs)))

for (i in 1:ncol(Motif_ov)) {
  for (j in 1:ncol(Motif_ov)) {
    Motif_ov[j, i] <- sum(overlapsAny(
      Motif_hits[Motif_hits$motif == rownames(Motif_ov)[j]],
      Motif_hits[Motif_hits$motif == rownames(Motif_ov)[i]]
    ))
    Motif_ov[j, i] <- Motif_ov[j, i] / min(
      sum(Motif_hits$motif == rownames(Motif_ov)[j]),
      sum(Motif_hits$motif == rownames(Motif_ov)[i])
    )
  }
}

Motif_t <- cutree(hclust(as.dist(1 - Motif_ov)), h = 0.8)

FinalMotifs <- tapply(Motifs[names(Motif_t)], Motif_t, merge_motifs)
FinalMotifs <- trim_motifs(FinalMotifs, 0.25)
names(FinalMotifs) <- paste0("M", seq_along(FinalMotifs))
for (i in seq_along(FinalMotifs)) {
  FinalMotifs[[i]]["name"] <- paste0("M", i)
  FinalMotifs[[i]]["altname"] <- FinalMotifs[[i]]["consensus"]
}

# Export final motifs

write_meme(FinalMotifs, "Final_motifs.meme")

# Calculate enrichment of final motifs and similarity with known TFBSs

system("bash Part05_motifs_helper2.sh")

ACR_sea <- list(
  readr::read_tsv("sea/sea_A1/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_A2/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_A3/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_A4/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_A5/sea.tsv", comment = "#")
)
names(ACR_sea) <- paste0("A", 1:5)

TSS_sea <- list(
  readr::read_tsv("sea/sea_C1/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_C2/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_C3/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_C4/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_C5/sea.tsv", comment = "#"),
  readr::read_tsv("sea/sea_C6/sea.tsv", comment = "#")
)
names(TSS_sea) <- paste0("C", 1:6)

na0 <- function(x) { x[is.na(x)] <- 0 ; x }
na1 <- function(x) { x[is.na(x)] <- 1 ; x }

# Create cluster enrichment tables

ACR_enr <- as.matrix(data.frame(row.names = names(FinalMotifs),
  A1 = -log10(na1(structure(ACR_sea$A1$QVALUE, names = ACR_sea$A1$ID)[names(FinalMotifs)])),
  A2 = -log10(na1(structure(ACR_sea$A2$QVALUE, names = ACR_sea$A2$ID)[names(FinalMotifs)])),
  A3 = -log10(na1(structure(ACR_sea$A3$QVALUE, names = ACR_sea$A3$ID)[names(FinalMotifs)])),
  A4 = -log10(na1(structure(ACR_sea$A4$QVALUE, names = ACR_sea$A4$ID)[names(FinalMotifs)])),
  A5 = -log10(na1(structure(ACR_sea$A5$QVALUE, names = ACR_sea$A5$ID)[names(FinalMotifs)]))
))

TSS_enr <- as.matrix(data.frame(row.names = names(FinalMotifs),
  C1 = -log10(na1(structure(TSS_sea$C1$QVALUE, names = TSS_sea$C1$ID)[names(FinalMotifs)])),
  C2 = -log10(na1(structure(TSS_sea$C2$QVALUE, names = TSS_sea$C2$ID)[names(FinalMotifs)])),
  C3 = -log10(na1(structure(TSS_sea$C3$QVALUE, names = TSS_sea$C3$ID)[names(FinalMotifs)])),
  C4 = -log10(na1(structure(TSS_sea$C4$QVALUE, names = TSS_sea$C4$ID)[names(FinalMotifs)])),
  C5 = -log10(na1(structure(TSS_sea$C5$QVALUE, names = TSS_sea$C5$ID)[names(FinalMotifs)])),
  C6 = -log10(na1(structure(TSS_sea$C6$QVALUE, names = TSS_sea$C6$ID)[names(FinalMotifs)]))
))

# Do a final filtering of the motifs, keeping those enriched in at least one dataset

ToKeep <- rownames(ACR_enr)[rowMaxs(ACR_enr) >= 2 | rowMaxs(TSS_enr) >= 2]
ACR_enr <- ACR_enr[ToKeep, ]
TSS_enr <- TSS_enr[ToKeep, ]

write_meme(FinalMotifs[ToKeep], "Final_motifs_filtered.meme")
write.table(ACR_enr, "Motif_enrichment_ACRs.tsv", quote = FALSE, sep = "\t")
write.table(TSS_enr, "Motif_enrichment_TSSs.tsv", quote = FALSE, sep = "\t")

TomTom <- as.data.frame(readr::read_tsv("tomtom/tomtom.tsv", comment = "#"))
TomTom <- TomTom[TomTom$Query_ID %in% ToKeep, ]

TomTom <- TomTom[order(TomTom$`E-value`), ]

readr::write_tsv(TomTom, "Motif_comparisons_known.tsv")

# Motif tables for the figures -------------------------------------------------

# Promoters as -400/+100 bp around the TSS, as made by bedtools slop

FigTSSs501 <- promoters(import("tss.size1.bed"), 400, 101)
write.table(data.frame(as.character(seqnames(FigTSSs501)), start(FigTSSs501) - 1L,
    end(FigTSSs501), FigTSSs501$name, FigTSSs501$score, as.character(strand(FigTSSs501))),
  "tss.size1_u400d100.bed", quote = FALSE, sep = "\t",
  row.names = FALSE, col.names = FALSE)

# Top 2000 TSSs and ACRs per cluster (by their highest mean CPM in any sample)

FigTSSA <- readr::read_tsv("csRNAseq_clustering/Clusters_anno.tsv")
FigACRA <- readr::read_tsv("ATACseq_clustering3/Clusters_anno.tsv")

FigTSSM <- as.matrix(FigTSSA[, 11:22])
rownames(FigTSSM) <- FigTSSA$TSS
FigTSSMm <- as.matrix(data.frame(row.names = rownames(FigTSSM),
  DS = rowMeans(FigTSSM[, 1:2]),
  S24 = rowMeans(FigTSSM[, 3:4]),
  S72 = rowMeans(FigTSSM[, 5:6]),
  L6 = rowMeans(FigTSSM[, 7:8]),
  L26 = rowMeans(FigTSSM[, 9:10]),
  L57 = rowMeans(FigTSSM[, 11:12])
))

FigACRM <- as.matrix(FigACRA[, 13:20])
rownames(FigACRM) <- FigACRA$PeakID
FigACRMm <- as.matrix(data.frame(row.names = rownames(FigACRM),
  DS = rowMeans(FigACRM[, 1:2]),
  L6 = rowMeans(FigACRM[, 3:4]),
  L26 = rowMeans(FigACRM[, 5:6]),
  L57 = rowMeans(FigACRM[, 7:8])
))

FigTSSTop <- lapply(1:6, function(x) {
  names(sort(rowMaxs(FigTSSMm[FigTSSA$TSS[FigTSSA$Cluster == x], ], useNames = TRUE),
    decreasing = TRUE)[1:2000])
})
FigACRTop <- lapply(1:5, function(x) {
  names(sort(rowMaxs(FigACRMm[FigACRA$PeakID[FigACRA$Cluster == x], ], useNames = TRUE),
    decreasing = TRUE)[1:2000])
})

FigACRs500 <- import("peaks_final_500bp.bed")
seqlengths(FigACRs500) <- seqlengths(BSgenome.Athaliana.TAIR.TAIR9)[1:5]
FigACRs500 <- trim(FigACRs500)

FigSeqs <- list()
FigSeqs$Cbkg <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, FigTSSs501)
names(FigSeqs$Cbkg) <- FigTSSs501$name
FigSeqs$Abkg <- getSeq(BSgenome.Athaliana.TAIR.TAIR9, FigACRs500)
names(FigSeqs$Abkg) <- FigACRs500$name
for (i in 1:6) {
  FigSeqs[[paste0("C", i)]] <- FigSeqs$Cbkg[FigTSSs501$name %in% FigTSSTop[[i]]]
}
for (i in 1:5) {
  FigSeqs[[paste0("A", i)]] <- FigSeqs$Abkg[FigACRs500$name %in% FigACRTop[[i]]]
}

dir.create("fasta2", showWarnings = FALSE)
writeXStringSet(FigSeqs$Cbkg, "fasta2/tss_bkg.fa")
writeXStringSet(FigSeqs$Abkg, "fasta2/atac_bkg.fa")
for (i in 1:6) {
  writeXStringSet(FigSeqs[[paste0("C", i)]], paste0("fasta2/tss_C", i, ".fa"))
}
for (i in 1:5) {
  writeXStringSet(FigSeqs[[paste0("A", i)]], paste0("fasta2/atac_A", i, ".fa"))
}

# Find motifs with streme and check their centrality with centrimo

system("bash Part05_motifs_helper3.sh streme")

# Keep central motifs with enough hits, removing redundant ones

fig_proc_streme <- function(x) {
  mots <- read_meme(paste0("streme_top2k/", x, "/streme.txt"))
  names(mots) <- sapply(mots, function(y) y["name"])
  cen <- readr::read_tsv(paste0("streme_top2k/", x, "_centrimo/centrimo.tsv"), comment = "#")
  mots_good <- cen$motif_id[cen$`E-value` <= 1e-10]
  if (!length(mots_good)) return(NULL)
  mots2 <- mots[mots_good]
  hits2 <- scan_sequences(mots2, FigSeqs[[x]], threshold = 1e-4,
    RC = TRUE, return.granges = TRUE, no.overlaps = TRUE)
  nhits <- table(hits2$motif)
  mots2 <- mots2[nhits[names(mots2)] >= 200]
  ovs <- matrix(0, ncol = length(mots2), nrow = length(mots2),
    dimnames = list(names(mots2), names(mots2)))
  for (i in 1:ncol(ovs)) {
    for (j in 1:nrow(ovs)) {
      ovs[j, i] <- sum(overlapsAny(
        hits2[hits2$motif == rownames(ovs)[j]],
        hits2[hits2$motif == colnames(ovs)[i]]
      ))
      ovs[j, i] <- ovs[j, i] / min(
        sum(hits2$motif == rownames(ovs)[j]),
        sum(hits2$motif == colnames(ovs)[i])
      )
    }
  }
  ovs[upper.tri(ovs)] <- NA
  ovs <- reshape2::melt(ovs)
  ovs[[1]] <- as.character(ovs[[1]])
  ovs[[2]] <- as.character(ovs[[2]])
  ovs <- ovs[!is.na(ovs[[3]]), ]
  ovs <- ovs[ovs[[1]] != ovs[[2]], ]
  ovs <- ovs[order(ovs[[3]], decreasing = TRUE), ]
  ovs$Cen1 <- structure(cen$`E-value`, names = cen$motif_id)[ovs$Var1]
  ovs$Cen2 <- structure(cen$`E-value`, names = cen$motif_id)[ovs$Var2]
  ovs$Enr1 <- structure(sapply(mots, function(y) y["eval"]), names = names(mots))[ovs$Var1]
  ovs$Enr2 <- structure(sapply(mots, function(y) y["eval"]), names = names(mots))[ovs$Var2]
  ovs <- ovs[ovs$value >= 0.1, ]
  for (i in seq_len(nrow(ovs))) {
    if (ovs$Cen1[i] < ovs$Cen2[i] && ovs$Enr1[i] < ovs$Enr2[i]) {
      mots2 <- mots2[names(mots2) != ovs$Var2[i]]
    } else if (ovs$Cen2[i] < ovs$Cen1[i] && ovs$Enr2[i] < ovs$Enr1[i]) {
      mots2 <- mots2[names(mots2) != ovs$Var1[i]]
    } else {
      if (all(ovs$Cen1[i] < c(ovs$Cen2[i], ovs$Enr1[i], ovs$Enr2[i])) ||
        all(ovs$Enr1[i] < c(ovs$Cen1[i], ovs$Cen2[i], ovs$Enr2[i]))) {
        mots2 <- mots2[names(mots2) != ovs$Var2[i]]
      } else if (all(ovs$Cen2[i] < c(ovs$Cen1[i], ovs$Enr1[i], ovs$Enr2[i])) ||
        all(ovs$Enr2[i] < c(ovs$Cen1[i], ovs$Cen2[i], ovs$Enr1[i]))) {
        mots2 <- mots2[names(mots2) != ovs$Var1[i]]
      }
    }
  }
  mots2
}

fig_motif_ovs <- function(mots, seqs) {
  names(mots) <- sapply(mots, function(x) x["name"])
  hits <- scan_sequences(mots, seqs, threshold = 0.0001, RC = TRUE,
    return.granges = TRUE, no.overlaps = TRUE)
  ovs <- matrix(0, ncol = length(mots), nrow = length(mots),
    dimnames = list(names(mots), names(mots)))
  for (i in 1:ncol(ovs)) {
    for (j in 1:nrow(ovs)) {
      ovs[j, i] <- sum(overlapsAny(
        hits[hits$motif == rownames(ovs)[j]],
        hits[hits$motif == colnames(ovs)[i]]
      ))
      ovs[j, i] <- ovs[j, i] / min(
        sum(hits$motif == rownames(ovs)[j]),
        sum(hits$motif == colnames(ovs)[i])
      )
    }
  }
  ovs
}

FigMotsA <- do.call(c, lapply(paste0("A", 1:5), fig_proc_streme))
FigMotsC <- do.call(c, lapply(paste0("C", 1:6), fig_proc_streme))
names(FigMotsA) <- sapply(FigMotsA, function(x) x["name"])
names(FigMotsC) <- sapply(FigMotsC, function(x) x["name"])

FigSeqsA <- c(FigSeqs$A1, FigSeqs$A2, FigSeqs$A3, FigSeqs$A4, FigSeqs$A5)
FigSeqsC <- c(FigSeqs$C1, FigSeqs$C2, FigSeqs$C3, FigSeqs$C4, FigSeqs$C5, FigSeqs$C6)

# Merge the ACR and promoter motifs separately, then together

FigOvsA <- fig_motif_ovs(FigMotsA, FigSeqsA)
FigHA <- cutree(hclust(as.dist(1 - FigOvsA)), h = 0.9)
FigMotsAFinal <- tapply(FigMotsA[names(FigHA)], FigHA, merge_motifs)

FigOvsC <- fig_motif_ovs(FigMotsC, FigSeqsC)
FigHC <- cutree(hclust(as.dist(1 - FigOvsC)), h = 0.9)
FigMotsCFinal <- tapply(FigMotsC[names(FigHC)], FigHC, merge_motifs)

for (i in seq_along(FigMotsAFinal)) {
  FigMotsAFinal[[i]]["name"] <- FigMotsAFinal[[i]]["consensus"]
}
for (i in seq_along(FigMotsCFinal)) {
  FigMotsCFinal[[i]]["name"] <- FigMotsCFinal[[i]]["consensus"]
}
names(FigMotsAFinal) <- sapply(FigMotsAFinal, function(x) x["name"])
names(FigMotsCFinal) <- sapply(FigMotsCFinal, function(x) x["name"])

FigOvsAll <- fig_motif_ovs(c(FigMotsAFinal, FigMotsCFinal), c(FigSeqsA, FigSeqsC))
FigHAll <- cutree(hclust(as.dist(1 - FigOvsAll)), h = 0.8)
FigMotsFinal <- tapply(c(FigMotsAFinal, FigMotsCFinal)[names(FigHAll)], FigHAll, merge_motifs)
FigMotsFinal <- trim_motifs(FigMotsFinal, 0.25)
names(FigMotsFinal) <- sapply(FigMotsFinal, function(x) x["consensus"])
for (i in seq_along(FigMotsFinal)) {
  FigMotsFinal[[i]]["name"] <- FigMotsFinal[[i]]["consensus"]
  FigMotsFinal[[i]]["altname"] <- paste0("M", i)
}

write_meme(FigMotsFinal, "all_final.meme", overwrite = TRUE)

# Motif enrichment (sea), similarity with known TFBSs (tomtom) and genome scan

system("bash Part05_motifs_helper3.sh sea")
system("bash Part05_motifs_helper3.sh scan")

FigSeaA <- lapply(paste0("streme_final/all_A", 1:5, "_sea/sea.tsv"),
  function(x) readr::read_tsv(x, comment = "#"))
FigSeaC <- lapply(paste0("streme_final/all_C", 1:6, "_sea/sea.tsv"),
  function(x) readr::read_tsv(x, comment = "#"))
names(FigSeaA) <- paste0("A", 1:5)
names(FigSeaC) <- paste0("C", 1:6)

FigEnrA <- as.matrix(data.frame(row.names = names(FigMotsFinal),
  lapply(FigSeaA, function(x) -log10(na1(structure(x$QVALUE, names = x$ID)[names(FigMotsFinal)])))))
FigEnrC <- as.matrix(data.frame(row.names = names(FigMotsFinal),
  lapply(FigSeaC, function(x) -log10(na1(structure(x$QVALUE, names = x$ID)[names(FigMotsFinal)])))))

FigToKeep <- rownames(FigEnrA)[rowMaxs(FigEnrA) >= 2 | rowMaxs(FigEnrC) >= 2]
FigEnrA <- FigEnrA[FigToKeep, ]
FigEnrC <- FigEnrC[FigToKeep, ]

saveRDS(list(ATAC = FigEnrA, csRNA = FigEnrC), "enrichment_matrix_final.RDS")

# Best match per TF, and correlation of TF expression with motif enrichment

FigA11 <- import("Araport11_novel.gff3")
FigGene2Sym <- structure(FigA11$symbol, names = FigA11$ID)
FigGene2Sym[is.na(FigGene2Sym)] <- names(FigGene2Sym)[is.na(FigGene2Sym)]
rm(FigA11)

FigTFs <- sapply(read_meme("PlantTFDB_Ath_TF_binding_motifs.meme"), function(x) x["name"])

FigTom <- as.data.frame(readr::read_tsv("streme_final/all_final_tomtom/tomtom.tsv", comment = "#"))
FigTom <- FigTom[FigTom$Query_ID %in% FigToKeep, ]
FigTom <- FigTom[order(FigTom$`E-value`), ]
FigTom <- FigTom[!duplicated(FigTom$Target_ID), ]
FigTom$Target_Symbol <- FigGene2Sym[FigTom$Target_ID]

FigTSSAtf <- FigTSSA[FigTSSA$NearestGene %in% FigTFs, ]
FigTSSAtf <- FigTSSAtf[order(rowSums(as.matrix(FigTSSAtf[, 11:22])), decreasing = TRUE), ]
FigTSSAtf <- FigTSSAtf[!duplicated(FigTSSAtf$NearestGene), ]
FigTF2TSS <- structure(FigTSSAtf$TSS, names = FigTSSAtf$NearestGene)

FigTSSMmZ <- t(scale(t(FigTSSMm)))
FigACRMmZ <- t(scale(t(FigACRMm)))
FigTFMmZ <- FigTSSMmZ[FigTF2TSS, ]
rownames(FigTFMmZ) <- names(FigTF2TSS)

FigTSSMcls <- matrix(0, ncol = 6, nrow = 6,
  dimnames = list(paste0("C", 1:6), colnames(FigTSSMm)))
FigACRMcls <- matrix(0, ncol = 4, nrow = 5,
  dimnames = list(paste0("A", 1:5), colnames(FigACRMm)))
for (i in 1:5) {
  FigACRMcls[i, ] <- colMeans(FigACRMmZ[FigACRA$PeakID[FigACRA$Cluster == i], ], na.rm = TRUE)
}
for (i in 1:6) {
  FigTSSMcls[i, ] <- colMeans(FigTSSMmZ[FigTSSA$TSS[FigTSSA$Cluster == i], ], na.rm = TRUE)
}

FigTFTSSMem <- cor(t(FigTFMmZ), t(FigTSSMcls))
FigTFACRMem <- cor(t(FigTFMmZ[, c("DS", "L6", "L26", "L57")]), t(FigACRMcls))

FigTom$AtacCor <- NA_real_
FigTom$TssCor <- NA_real_
FigTom <- FigTom[FigTom$Target_ID %in% rownames(FigTFMmZ), ]
for (i in 1:nrow(FigTom)) {
  FigTom$AtacCor[i] <- cor(FigTFACRMem[FigTom$Target_ID[i], ], FigEnrA[FigTom$Query_ID[i], ])
  FigTom$TssCor[i] <- cor(FigTFTSSMem[FigTom$Target_ID[i], ], FigEnrC[FigTom$Query_ID[i], ])
}
FigTom$Target_Cluster <- structure(FigTSSAtf$Cluster, names = FigTSSAtf$NearestGene)[FigTom$Target_ID]
FigTom$Target_MaxExpr <- sapply(FigTom$Target_ID, function(x) max(FigTSSMm[FigTF2TSS[x], ]))
FigTom$Target_PeakExpr <- sapply(FigTom$Target_ID, function(x) which.max(FigTSSMm[FigTF2TSS[x], ]))
FigTom$Target_MaxZ <- sapply(FigTom$Target_ID, function(x) max(FigTFMmZ[x, ]))

readr::write_tsv(FigTom, "tomtom_final_motifs.tsv")

# ABI5 peaks from seedlings and leaves (ReMap 2022 non-redundant peaks)

FigReMap <- import("remap2022_nr_macs2_TAIR10_v1_0.bed.gz")
FigReMap <- FigReMap[grepl("seedling", FigReMap$name) | grepl("leaves", FigReMap$name)]
FigReMap <- FigReMap[!grepl("inflo|shoot|etiol", FigReMap$name)]
FigReMap$name <- vapply(strsplit(FigReMap$name, ":", TRUE), function(x) x[1], "")
FigABI5 <- sort(FigReMap[FigReMap$name == "ABI5"])
write.table(data.frame(as.character(seqnames(FigABI5)), start(FigABI5) - 1L,
    end(FigABI5), FigABI5$name, 1L, "."),
  "ABI5.bed", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)

