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

