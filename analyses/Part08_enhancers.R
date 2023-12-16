# Benjamin J.M. Tremblay, June 2023
#
# Putative unidirectional and bidirectional enhancers

library(rtracklayer)
library(GenomicFeatures)

# Load in TSS/ACR data

TSSs1 <- import("tss.size1.bed")
TSSs <- import("tss.shrunk.bed")

ACRs <- import("peaks_final.bed")

ACR_anno <- readr::read_tsv("ACR.annotations.tsv")

DetTx <- import("detected.transcripts.bed")

Antisenses <- readr::read_tsv("All_antisenses.tsv")

Divergent <- readr::read_tsv("Divergent_TSS.tsv")

TSS_anno <- readr::read_tsv("TSS.annotations.tsv")
TSS_anno$NearestGeneType[is.na(TSS_anno$NearestGeneType)] <- "none"

BiDir_ncnc <- readr::read_tsv("Bidirectional_noncoding_TSS.tsv")
GR_ncnc <- import("Bidirectional_noncoding_TSS.bed")

mcols(GR_ncnc) <- as(as.data.frame(BiDir_ncnc), "DataFrame")

# Start looking for enhancers: intragenic/intergenic bidirectional non-coding TSSs
# and intergenic non-antisense non-divergent unidirectional non-coding TSSs

Enhancers <- TSS_anno$TSS[is.na(TSS_anno$NearestTranscriptType) |
  TSS_anno$NearestTranscriptType == "putative_rna"]
Enhancers <- Enhancers[structure(TSS_anno$TSS_type, names = TSS_anno$TSS)[Enhancers] %in%
  c("unknown_intergenic", "known_tss", "unknown_antisense")]
Enhancers <- Enhancers[!Enhancers %in% c(Antisenses$pcTSS, Antisenses$asTSS)]
Enhancers <- Enahncers[!Enhancers %in% c(Divergent_TSS$TSS1, Divergent_TSS$TSS2)]
Enhancers1 <- TSSs1[TSSs1$name %in% Enhancers]

EnhancersACR <- ACRs[overlapsAny(ACRs, Enhancers1)]
ACR_anno_enh <- ACR_anno[ACR_anno$PeakID %in% EnhancersACR$name, ]
ACR_anno_enh <- ACR_anno_enh[ACR_anno_enh$Type != "Intragenic", ]
EnhancersACR <- EnhancersACR[EnhancersACR %in% ACR_anno_enh$PeakID]

EnhancersACR <- EnhancersACR[!overlapsAny(EnhancersACR, TSSs1[!TSSs1$name %in% Enhancers1$name])]
EnhancersACR <- EnhancersACR[!overlapsAny(EnhancersACR, DetTx[grepl("^AT", DetTx$name)])]

# Resize bidirectional non-coding TSS enhancers if their ACR is bigger

GR_ncnc$ACR_start <- structure(start(ACRs), names = ACRs$name)[GR_ncnc$ACR]
GR_ncnc$ACR_end <- structure(end(ACRs), names = ACRs$name)[GR_ncnc$ACR]
start(GR_ncnc) <- pmin(start(GR_ncnc), GR_ncnc$ACR_start, na.rm = TRUE)
end(GR_ncnc) <- pmax(end(GR_ncnc), GR_ncnc$ACR_end, na.rm = TRUE)

EnhancersACR <- EnhancersACR[!EnhancersACR$name %in% GR_ncnc$ACR]
EnhancersACR$TSS <- as.character(tapply(Enhancers1$name[
  subjectHits(findOverlaps(EnhancersACR, Enhancers1))],
  queryHits(findOverlaps(EnhancersACR, Enhancers1)),
  function(x) paste0(x, collapse = ",")))

EnhancersACR$Intragenic <- FALSE
EnhancersACR$ACR <- EnhancersACR$name
GR_ncnc$TSS <- paste0(GR_ncnc$TSS1, ",", GR_ncnc$TSS2)
GR_ncnc$BiDir <- TRUE
EnhancersACR$BiDir <- FALSE

FinalEnhancers <- GRanges(c(seqnames(GR_ncnc), seqnames(EnhancersACR)),
  IRanges(c(start(GR_ncnc), start(EnhancersACR)), c(end(GR_ncnc), end(EnhancersACR))),
  TSS = c(GR_ncnc$TSS, EnhancersACR$TSS), ACR = c(GR_ncnc$ACR, EnhancersACR$ACR),
  Intragenic = c(GR_ncnc$Intragenic, EnhancersACR$Intragenic),
  BiDir = c(GR_ncnc$BiDir, EnhancersACR$BiDir))
FinalEnhancers <- sort(FinalEnhancers)

DetTx_pc <- DetTx[DetTx$name %in% TSS_anno$NearestTranscript[!is.na(TSS_anno$NearestGeneType) &
  TSS_anno$NearestGeneType == "protein_coding"]]
ACRs_inter <- ACRs[ACRs$name %in% ACR_anno$PeakID[!ACR_anno$Type %in% c("Intragenic", "Promoter")]]

TxDb <- makeTxDbFromGFF("Araport11_novel.gff3")
A11 <- import("Araport11_novel.gff3")
pcGenes <- A11$ID[!is.na(A11$locus_type) & A11$locus_type == "protein_coding"]
DetTx_pc <- DetTx_pc[gsub("[.]\\d+", "", DetTx_pc$name) %in% pcGenes]

# Get unidirectional enhancers

Enhancers1_old <- Enhancers1
Enhancers1 <- Enhancers1[!Enhancers1$name %in% c(Divergent$TSS1, Divergent$TSS2,
  unlist(strsplit(FinalEnhancers$TSS, ",", fixed = TRUE)))]
Enhancers1 <- reduce(resize(Enhancers1, 500, "center"))
Enhancers1 <- Enhancers1[!overlapsAny(Enhancers1, DetTx[grepl("^AT", DetTx$name)],
  ignore.strand = TRUE)]
Enhancers1$TSS <- as.character(tapply(Enhancers1$name[subjectHits(findOverlaps(,
  Enhancers1, Enhancers1_old))], queryHits(findOverlaps(Enhancers1, Enhancers1_old)),
  function(x) paste0(x, collapse = ",")))
Enhancers1 <- Enhancers1[!overlapsAny(Enhancers1, ACRs)]
Enhancers1 <- Enhancers1[!overlapsAny(Enhancers1, transcripts(TxDb)[
  grepl("^AT", transcripts(TxDb)$tx_name)], ignore.strand = TRUE)]

Enhancers1$ACR <- NA_character_
Enhancers1$Intragenic <- FALSE
Enhancers1$BiDir <- FALSE

# Final list of enhancers

FinalEnhancers <- sort(c(FinalEnhancers, Enhancers1))
strand(FinalEnhancers) <- "*"

FinalEnhancers$name <- paste0("ENH", as.character(seqnames(FinalEnhancers)),
  "_", round((start(FinalEnhancers) + end(FinalEnhancers)) / 2))
FinalEnhancers <- FinalEnhancers[!duplicated(FinalEnhancers$name)]

# Perform manual adjustment of enhancers that overlap promoters

TSSs_pc <- TSSs[TSSs$name %in% TSS_anno$TSS[!is.na(TSS_anno$NearestGeneType) &
  TSS_anno$NearestGeneType == "protein_coding"]]
export(FinalEnhancers[overlapsAny(FinalEnhancers, TSSs_pc, ignore.strand = TRUE)],
  "Enhancers_to_fix.bed")

# Read adjusted enhancers back into R

ToFix <- import("Enhancers_to_fix.bed")
start(FinalEnhancers)[FinalEnhancers$name %in% ToFix$name] <- start(ToFix)
end(FinalEnhancers)[FinalEnhancers$name %in% ToFix$name] <- end(ToFix)

# Fix names again

FinalEnhancers$name <- paste0("ENH", as.character(seqnames(FinalEnhancers)),
  "_", round((start(FinalEnhancers) + end(FinalEnhancers)) / 2))
FinalEnhancers <- sort(FinalEnhancers)

# Export final enhancers

export(FinalEnhancers, "Final_enhancers.bed")
readr::write_tsv(as.data.frame(FinalEnhancers), "Final_enhancers.tsv")

# Also a version where the enhancers are at least 500 bp

FinalEnhancers[width(FinalEnhancers) < 500] <- resize(
  FinalEnhancers[width(FinalEnhancers) < 500], 500, "center")

export(FinalEnhancers, "Final_enhancers_500bp.bed")
readr::write_tsv(as.data.frame(FinalEnhancers), "Final_enhancers_500bp.tsv")

# Quantify enhancer activity

system("bash Part08_enhancers_helper.sh")

DS_p <- readr::read_tsv("Enhancer_activities/DS_plus.txt", col_names = FALSE)
DS_m <- readr::read_tsv("Enhancer_activities/DS_minus.txt", col_names = FALSE)

S24_p <- readr::read_tsv("Enhancer_activities/S24_plus.txt", col_names = FALSE)
S24_m <- readr::read_tsv("Enhancer_activities/S24_minus.txt", col_names = FALSE)

S72_p <- readr::read_tsv("Enhancer_activities/S72_plus.txt", col_names = FALSE)
S72_m <- readr::read_tsv("Enhancer_activities/S72_minus.txt", col_names = FALSE)

L6_p <- readr::read_tsv("Enhancer_activities/L6_plus.txt", col_names = FALSE)
L6_m <- readr::read_tsv("Enhancer_activities/L6_minus.txt", col_names = FALSE)

L26_p <- readr::read_tsv("Enhancer_activities/L26_plus.txt", col_names = FALSE)
L26_m <- readr::read_tsv("Enhancer_activities/L26_minus.txt", col_names = FALSE)

L57_p <- readr::read_tsv("Enhancer_activities/L57_plus.txt", col_names = FALSE)
L57_m <- readr::read_tsv("Enhancer_activities/L57_minus.txt", col_names = FALSE)

EnhancerActivity <- as.matrix(data.frame(row.names = DS_p$X1,
  DS = DS_p$X4 + DS_m$X4,
  S24 = S24_p$X4 + S24_m$X4,
  S72 = S72_p$X4 + S72_m$X4,
  L6 = L6_p$X4 + L6_m$X4,
  L26 = L26_p$X4 + L26_m$X4,
  L57 = L57_p$X4 + L57_m$X4
))

mcols(FinalEnhancers) <- cbind(mcols(FinalEnhancers), EnhancerActivity[FinalEnhancers$name, ])

FinalEnhancers$MaxActivitiy <- apply(mcols(FinalEnhancers)[,
  c("DS", "S24", "S72", "L6", "L26", "L57")], 1,
  function(x) c("DS", "S24", "S72", "L6", "L26", "L57")[which.max(x)])

readr::write_tsv(as.data.frame(FinalEnhancers), "Final_enhancer_activity.tsv")

# Look for nearby (within 5 Kbp in either direction) correlating protein coding TSSs

TSSs1_pc <- TSSs1[TSSs1$name %in% TSS_anno$TSS[!is.na(TSS_anno$NearestGeneType) &
  TSS_anno$NearestGeneType == "protein_coding"]]

tssM <- as.matrix(TSS_anno[, 10:21])
rownames(tssM) <- TSS_anno$TSS
tssMm <- as.matrix(data.frame(row.names = rownames(tssM),
  DS = rowMeans(tssM[, 1:2]),
  S24 = rowMeans(tssM[, 3:4]),
  S72 = rowMeans(tssM[, 5:6]),
  L6 = rowMeans(tssM[, 7:8]),
  L26 = rowMeans(tssM[, 9:10]),
  L57 = rowMeans(tssM[, 11:12])
))

TSSs1_pc_ov <- TSSs1_pc[overlapsAny(TSSs1_pc, resize(FinalEnhancers, 10000, "center"))]

EnhancerActivityL2 <- log2(1 + EnhancerActivity)
tssMmL2 <- log2(1 + tssMm)

EnhancerPCC <- matrix(0, ncol = length(FinalEnhancers), nrow = length(TSSs1_pc_ov),
  dimnames = list(TSSs1_pc_ov$name, FinalEnhancers$name))

for (i in 1:nrow(EnhancerPCC)) {
  EnhancerPCC[i, ] <- cor(tssMmL2[rownames(EnhancerPCC)[i], ],
    t(EnhancerActivityL2[colnames(EnhancerPCC), ]))[1, ]
}

EnhancerPCC <- reshape2::melt(EnhancerPCC)
EnhancerPCC$Var1 <- as.character(EnhancerPCC$Var1)
EnhancerPCC$Var2 <- as.character(EnhancerPCC$Var2)
colnames(EnhancersPCC) <- c("pcTSS", "Enhancer", "PCC")

EnhancerPCC <- EnhancerPCC[substr(EnhancerPCC$pcTSS, 4, 4) == substr(EnhancerPCC$Enhancer, 4, 4), ]
EnhancerPCC <- EnhancerPCC[abs(as.integer(gsub("TSS\\d_", EnhancerPCC$pcTSS)) - 
  as.integer(gsub("ENH\\d_", "", EnhancerPCC$Enhancer))) <= 5000, ]
EnhancerPCC <- EnhancerPCC[!is.na(EnhancerPCC$PCC) & EnhancerPCC$PCC >= 0.5, ]
EnhancerPCC$Transcript <- structure(TSS_anno$NearestTranscript, names = TSS_anno$TSS)[
  EnhancerPCC$pcTSS]
EnhancerPCC$Gene <- structure(TSS_anno$NearestGene, names = TSS_anno$TSS)[
  EnhancerPCC$pcTSS]
EnhancerPCC$Symbol <- structure(A11$symbol, names = A11$ID)[EnhancerPCC$Gene]
EnhancerPCC$Symbol[is.na(EnhancerPCC$Symbol)] <- ""

# Export the results

readr::write_tsv(EnhancerPCC, "Enhancer_gene_links.tsv")

