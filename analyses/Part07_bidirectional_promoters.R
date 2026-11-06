# Benjamin J.M. Tremblay, June 2023
#
# Bidirectional promoters

library(rtracklayer)
library(GenomicFeatures)

# Load in TSS data

TSSs1 <- import("tss.size1.bed")
TSSs <- import("tss.shrunk.bed")

TSS_anno <- readr::read_tsv("TSS.annotations.tsv")
TSS_anno$NearestGeneType[is.na(TSS_anno$NearestGeneType)] <- "none"

tssM <- as.matrix(TSS_anno[, 10:25])
rownames(tssM) <- TSS_anno$TSS

TSSs1_rc <- TSSs1
strand(TSSs1_rc) <- c("+" = "-", "-" = "+")[as.character(strand(TSSs1))]

# Find all TSS pairs

BiDi_pairs <- follow(TSSs1, TSSs1_rc)

BiDir <- data.frame(
  TSS1 = TSSs1$name,
  seq1 = as.character(seqnames(TSSs1)),
  pos1 = start(TSSs1),
  strand1 = as.character(strand(TSSs1)),
  TSS2 = TSSs1$name[BiDi_pairs],
  seq2 = as.character(seqnames(TSSs1))[BiDi_pairs],
  pos2 = start(TSSs1)[BiDi_pairs],
  strand2 = as.character(strand(TSSs1))[BiDi_pairs]
)
BiDir <- BiDir[complete.cases(BiDir), ]

Dups <- mapply(function(x, y) paste0(sort(c(x, y)), collapse = "-"), BiDir$TSS1, BiDir$TSS2)
Dups <- duplicated(Dups)
BiDir <- BiDir[!Dups, ]

# Start annotating

BiDir$Distance <- mapply(function(x, y) distance(TSSs[TSSs$name == x], TSSs[TSSs$name == y],
  ignore.strand = TRUE), BiDir$TSS1, BiDir$TSS2)

BiDir$PCC <- NA_real_
for (i in 1:nrow(BiDir)) {
  BiDir$PCC[i] <- cor(tssM[BiDir$TSS1[i], ], tssM[BiDir$TSS2[i], ])
}

BiDir$Type1 <- structure(TSS_anno$NearestGeneType, names = TSS_anno$TSS)[BiDir$TSS1]
BiDir$Type2 <- structure(TSS_anno$NearestGeneType, names = TSS_anno$TSS)[BiDir$TSS2]

BiDir$Type <- mapply(function(x, y) {
  if (x == "protein_coding" && y == "protein_coding") "pc-pc"
  else if (x == "protein_coding" || y == "protein_coding") "nc-pc"
  else "nc-nc"
}, BiDir$Type1, BiDir$Type2)

# Export all pairs

readr::write_tsv(BiDir, "All_bidirectional_TSS.tsv")

# Filter to those within 500 bp

BiDir <- BiDir[BiDir$Distance <= 500, ]
BiDir_ncnc <- BiDir[BiDir$Type == "nc-nc", ]
BiDir_ncpc <- BiDir[BiDir$Type == "nc-pc", ]

readr::write_tsv(BiDir_ncpc, "Divergent_TSS.tsv")

# Extract bidirectional non-coding pairs

GR_ncnc <- with(BiDir_ncnc, GRanges(seq1,
  IRanges(
    pmin(
      structure(start(TSSs), names = TSSs$name)[BiDir_ncnc$TSS1],
      structure(start(TSSs), names = TSSs$name)[BiDir_ncnc$TSS2]
    ),
  pmax(
    structure(end(TSSs), names = TSSs$name)[BiDir_ncnc$TSS1],
    structure(end(TSSs), names = TSSs$name)[BiDir_ncnc$TSS2]
  )
), TSS1 = TSS1, TSS2 = TSS2, PCC = PCC, Distance = Distance))

# Find out if there is an overlapping ACR

ACRs <- import("peaks_final.bed")

GR_ncnc$ACR <- NA_character_
GR_ncnc$ACR[queryHits(findOverlaps(GR_ncnc, ACRs, ignore.strand = TRUE))] <- 
  ACRs$name[subjectHits(findOverlaps(GR_ncnc, ACRs, ignore.strand = TRUE))]

DetTx <- import("detected.transcripts.bed")

GR_ncnc$Intragenic <- overlapsAny(resize(GR_ncnc, 1, "center"), DetTx)

# Export bidirectional non-coding TSSs

readr::write_tsv(as.data.frame(GR_ncnc), "Bidirectional_noncoding_TSS.tsv")

GR_ncnc$name <- paste0("BIDIR", as.character(seqnames(GR_ncnc)), "_",
  round((end(GR_ncnc) + start(GR_ncnc)) / 2))
export(GR_ncnc, "Bidirectional_noncoding_TSS.bed")

# Bidirectional TSS table used in the figures -----------------------------------

FigBidir <- data.frame(
  TSS1 = TSSs1$name,
  seq1 = as.character(seqnames(TSSs1)),
  pos1 = start(TSSs1),
  strand1 = as.character(strand(TSSs1)),
  TSS2 = TSSs1$name[BiDi_pairs],
  seq2 = as.character(seqnames(TSSs1))[BiDi_pairs],
  pos2 = start(TSSs1)[BiDi_pairs],
  strand2 = as.character(strand(TSSs1))[BiDi_pairs]
)
FigBidir <- FigBidir[complete.cases(FigBidir), ]
FigBidir <- FigBidir[!duplicated(mapply(function(x, y) paste0(sort(c(x, y)), collapse = "-"),
  FigBidir$TSS1, FigBidir$TSS2)), ]

FigBidir$dist <- abs(FigBidir$pos1 - FigBidir$pos2)

FigBidir$PCC <- NA_real_
for (i in 1:nrow(FigBidir)) {
  FigBidir$PCC[i] <- cor(tssM[FigBidir$TSS1[i], ], tssM[FigBidir$TSS2[i], ])
}

FigBidir$Type1 <- structure(TSS_anno$NearestGeneType, names = TSS_anno$TSS)[FigBidir$TSS1]
FigBidir$Type2 <- structure(TSS_anno$NearestGeneType, names = TSS_anno$TSS)[FigBidir$TSS2]

FigBidir$Type <- mapply(function(x, y) paste0(sort(c(x, y)), collapse = "-"),
  FigBidir$Type1, FigBidir$Type2)
FigBidir$TypeV2 <- mapply(function(x, y) {
  if (x == "protein_coding" && y == "protein_coding") "pc-pc"
  else if (x == "protein_coding" || y == "protein_coding") "nc-pc"
  else "nc-nc"
}, FigBidir$Type1, FigBidir$Type2)

# Median ratio of the higher to the lower expressed TSS

median_ratio <- function(x, y, cols) {
  x_m <- tssM[x, cols]
  y_m <- tssM[y, cols]
  xy <- (x_m >= 0.1) & (y_m >= 0.1)
  median(pmax(x_m[xy], y_m[xy]) / pmin(x_m[xy], y_m[xy]))
}
FigBidir$MedianRatio <- mapply(median_ratio, FigBidir$TSS1, FigBidir$TSS2,
  MoreArgs = list(cols = 1:16))
FigBidir$MedianRatioCol0 <- mapply(median_ratio, FigBidir$TSS1, FigBidir$TSS2,
  MoreArgs = list(cols = 1:12))
FigBidir$MedianRatioExo <- mapply(median_ratio, FigBidir$TSS1, FigBidir$TSS2,
  MoreArgs = list(cols = 13:16))

FigBidir$distWide <- distance(TSSs[match(FigBidir$TSS1, TSSs$name)],
  TSSs[match(FigBidir$TSS2, TSSs$name)], ignore.strand = TRUE)

readr::write_tsv(FigBidir, "Bidir_all_pairs.tsv")

# Divergent transcripts from Thieffry et al. 2020 (Supplemental Data Set 1)

Prompts <- readxl::read_excel("TPC2019-LSB-00815R2_Supplemental_Data_Set_1.xlsx",
  sheet = "PROMPTs")
Prompts <- Prompts[Prompts$chr != "ChrC", ]
write.table(data.frame(
  chr = sub("Chr", "", Prompts$chr),
  start = as.integer(Prompts$start - 1),
  end = as.integer(Prompts$end),
  name = Prompts$geneID_anti,
  score = Prompts$TPMsupport,
  strand = Prompts$strand
), "Thieffry2020_prompts.bed",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
