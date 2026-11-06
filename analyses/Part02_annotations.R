# Benjamin J.M. Tremblay, June 2023
#
# Annotating TSSs and ACRs

library(rtracklayer)
library(GenomicFeatures)

TxDb <- makeTxDbFromGFF("Araport11_novel.gff3")
A11 <- import("Araport11_novel.gff3")

# csRNA-seq --------------------------------------------------------------------

TSSs <- import("tss.shrunk.bed")

Tx <- transcripts(TxDb)
Tx <- Tx[, "tx_name"]

# Clean up transcript/gene annotations

MIRs <- A11[grepl("mir", A11$Name, ignore.case = TRUE)]
MIRsV <- vapply(MIRs$Parent,
  function(x) if (!length(x)) NA_character_ else unique(gsub("[.]\\d+", "", x)), character(1))
MIRsV[is.na(MIRsV)] <- gsub("[.]\\d+", "", MIRs$Derives_from[is.na(MIRsV)])
MIRs <- structure(MIRsV, names = MIRs$Name)

Tx2Type <- structure(as.character(A11$type), names = A11$ID)

Tx2LocType <- structure(as.character(A11$locus_type), names = A11$ID)
Tx2LocType[is.na(Tx2LocType)] <- Tx2Type[is.na(Tx2LocType)]
Tx2LocType <- Tx2LocType[names(Tx2LocType) %in% Tx$tx_name]

Tx2GeneType <- structure(as.character(A11$locus_type), names = A11$ID)
Tx2GeneType[is.na(Tx2GeneType)] <- Tx2Type[is.na(Tx2GeneType)]
Tx2GeneType <- Tx2GeneType[names(Tx2GeneType) %in% gsub("[.]\\d+", "", Tx$tx_name)]

# Assign nearest annotations

Tx_TSS <- resize(Tx, 1)

FiveP <- unlist(fiveUTRsByTranscript(TxDb, use.names = TRUE))
FiveP$tx <- names(FiveP)

Anno <- data.frame(
  TSS = TSSs$name,
  NearestTranscript = NA_character_,
  NearestTranscriptType = NA_character_,
  NearestGene = NA_character_,
  NearestGeneType = NA_character_,
  AbsoluteDistance = NA_integer_,
  RelativeDistance = NA_real_
)

TSS_dist <- distanceToNearest(TSSs, Tx_TSS)

Anno$NearestTranscript <- Tx_TSS$tx_name[subjectHits(TSS_dist)]
Anno$NearestTranscriptType <- Tx2LocType[gsub("[.]\\d+", "", Anno$NearestTranscript)]
Anno$NearestGeneType <- Tx2GeneType[gsub("[.]\\d+", "", Anno$NearestTranscript)]
Anno$AbsoluteDistance <- mcols(TSS_dist)$distance
Anno$RelativeDistance <- Anno$AbsoluteDistance / width(Tx)[subjectHits(TSS_dist)]

# Refine distances

for (i in 1:nrow(Anno)) {
  if (Anno$AbsoluteDistance[i] == 0) next
  TSS_strand <- as.character(strand(TSSs[i]))
  Tx_strand <- as.character(strand(Tx[Tx$tx_name == Anno$NearestTranscript[i]]))
  if (Tx_strand != TSS_strand) stop(i) # This shouldn't happen
  TSS_pos <- c(start(TSSs[i]), end(TSSs[i]))
  Tx_pos <- if (TSS_strand == "+")
    start(Tx[Tx$tx_name == Anno$NearestTranscript[i]])
  else
    end(Tx[Tx$tx_name == Anno$NearestTranscript[i]])
  if (TSS_strand == "+" && Tx_pos > TSS_pos[2]) {
    Anno$RelativeDistance[i] <- Anno$RelativeDistance[i] * -1
  } else if (TSS_strand == "-" && Tx_pos < TSS_pos[2]) {
    Anno$RelativeDistance[i] <- Anno$RelativeDistance[i] * -1
  }
}

Anno$AbsoluteDistance[Anno$RelativeDistance < 0] <-
  Anno$AbsoluteDistance[Anno$RelativeDistance < 0] * -1

Anno$NearestTranscriptType <- Tx2LocType[Anno$NearestTranscript]
Anno$NearestTranscriptType[grepl("mir", Anno$NearestTranscript, ignore.case = TRUE)] <- "miRNA"
Anno$NearestTranscriptType[grepl("HOMER", Anno$NearestTranscript) | grepl("STRINGTIE", Anno$NearestTranscript)] <- 
  "putative_rna"

Anno$NearestTranscriptSize <- width(Tx)[subjectHits(TSS_dist)]

# Filter, refine assignments

filter_anno <- function(x, tofilter) {
  x$NearestTranscript[tofilter] <- NA_character_
  x$NearestTranscriptType[tofilter] <- NA_character_
  x$NearestGeneType[tofilter] <- NA_character_
  x$AbsoluteDistance[tofilter] <- NA
  x$RelativeDistance[tofilter] <- NA
  x
}

Anno <- filter_anno(Anno, Anno$AbsoluteDistance < -500)
Anno <- filter_anno(Anno, Anno$RelativeDistance > 1)
Anno <- filter_anno(Anno, Anno$AbsoluteDistance > 500 &
  (!is.na(Anno$NearestTranscriptType) & Anno$NearestTranscriptType != "mRNA"))
for (i in 1:nrow(Anno)) {
  if (!is.na(Anno$NearestGeneType[i]) && Anno$NearestGeneType[i] == "protein_coding") {
    TSS_i <- TSSs[TSSs$name == Anno$TSS[i]]
    FiveP_i <- FiveP[FiveP$tx == Anno$NearestTranscript[i]]
    if (!overlapsAny(TSS_i, FiveP_i) &&
      (Anno$AbsoluteDistance[i] < -500 || Anno$AbsoluteDistance[i] > 200)) {
      Anno <- filter_anno(Anno, i)
    }
  }
}

# Peform manual adjustments

if (file.exists("TSS_manual_assignment.txt")) {
  ma <- readr::read_tsv("TSS_manual_assignment.txt", col_names = FALSE)
  for (i in 1:nrow(ma)) {
    if (ma[[2]][i] == "none") {
      Anno <- filter_anno(Anno, Anno$TSS == ma[[1]][i])
    } else {
      if (is.na(Anno$NearestTranscript[Anno$TSS == ma[[1]][i]]) || 
          Anno$NearestTranscript[Anno$TSS == ma[[1]][i]] != ma[[2]][i]) {
        Anno$NearestTranscript[Anno$TSS == ma[[1]][i]] <- ma[[2]][i]
        Anno$NearestTranscriptType[Anno$TSS == ma[[1]][i]] <- Tx2LocType[ma[[2]][i]]
        Anno$NearestGeneType[Anno$TSS == ma[[1]][i]] <- Tx2GeneType[gsub("[.]\\d+", "", ma[[2]][i])]
        if (grepl("mir", Anno$NearestTranscript[Anno$TSS == ma[[1]][i]], ignore.case = TRUE)) {
          Anno$NearestTranscriptType[Anno$TSS == ma[[1]][i]] <- "mirna"
        } else if (grepl("HOMER", Anno$NearestTranscriptType[Anno$TSS == ma[[1]][i]])) {
          Anno$NearestTranscriptType[Anno$TSS == ma[[1]][i]] <- "putative_rna"
        } else if (grepl("STRINGTIE", Anno$NearestTranscriptType[Anno$TSS == ma[[1]][i]])) {
          Anno$NearestTranscriptType[Anno$TSS == ma[[1]][i]] <- "putative_rna"
        }
      }
      Tx_i <- Tx[Tx$tx_name == ma[[2]][i]]
      Anno$AbsoluteDistance[Anno$TSS == ma[[1]][i]] <- distance(Tx_i, TSSs[TSSs$name == ma[[1]][i]])
      Anno$NearestTranscriptSize[Anno$TSS == ma[[1]][i]] <- width(Tx_i)
      Anno$RelativeDistance[Anno$TSS == ma[[1]][i]] <- Anno$AbsoluteDistance[Anno$TSS == ma[[1]][i]] /
        Anno$NearestTranscriptSize[Anno$TSS == ma[[1]][i]]
      if (!overlapsAny(Tx_i, TSSs[TSSs$name == ma[[1]][i]])) {
        Anno$AbsoluteDistance[Anno$TSS == ma[[1]][i]] <- -Anno$AbsoluteDistance[Anno$TSS == ma[[1]][i]]
        Anno$RelativeDistance[Anno$TSS == ma[[1]][i]] <- -Anno$RelativeDistance[Anno$TSS == ma[[1]][i]]
      }
    }
  }
}

Anno <- Anno[order(abs(Anno$RelativeDistance)), ]

Anno$NearestTranscriptSize[is.na(Anno$AbsoluteDistance)] <- NA

Anno$TSS_type <- "known_tss"
Anno$TSS_type[is.na(Anno$NearestTranscript)] <- NA

names(TSSs) <- TSSs$name

# Annotate TSSs with no transcript

for (i in 1:nrow(Anno)) {
  if (!is.na(Anno$NearestTranscript[i])) next
  if (!countOverlaps(TSSs[Anno$TSS[i]], Tx, ignore.strand = TRUE)) {
    Anno$TSS_type[i] <- "unknown_intergenic"
  }
}

for (i in 1:nrow(Anno)) {
  if (!is.na(Anno$TSS_type[i])) next
  if (countOverlaps(TSSs[Anno$TSS[i]], Tx)) {
    Anno$TSS_type[i] <- "unknown_genic"
  }
}

TxRC <- Tx
strand(TxRC) <- c("+" = "-", "-" = "+", "*" = "*")[as.character(strand(Tx))]

for (i in 1:nrow(Anno)) {
  if (!is.na(Anno$TSS_type[i])) next
  if (countOverlaps(TSSs[Anno$TSS[i]], TxRC)) {
    Anno$TSS_type[i] <- "unknown_antisense"
  }
}

# Final cleanup

Anno$NearestGeneType[is.na(Anno$NearestGeneType)] <- Anno$NearestTranscriptType[is.na(Anno$NearestGeneType)]
Anno$NearestTranscriptType[Anno$NearestTranscriptType == "miRNA_primary_transcript"] <- "miRNA"
Anno$NearestGeneType[Anno$NearestGeneType == "mirna" | Anno$NearestTranscriptType == "mirna"] <- "miRNA"
Anno$NearestTranscriptType[Anno$NearestTranscriptType %in% "mirna"] <- "miRNA"
# Novel transcripts from the RNA-seq are labelled as putative RNAs
Anno$NearestGeneType[Anno$NearestGeneType %in% "transcript"] <- "putative_rna"

Anno$NearestGene <- gsub("[.]\\d+", "", Anno$NearestTranscript)

Anno$NearestGene[grepl("mir", Anno$NearestTranscript, ignore.case = TRUE)] <- 
  MIRs[Anno$NearestTranscript[grepl("mir", Anno$NearestTranscript, ignore.case = TRUE)]]

# Add csRNA-seq CPMs and RNA-seq TPMs

AnnoCPMs <- as.data.frame(readr::read_tsv("csRNAseq.quantification.cpm.txt"))
AnnoTPMs <- as.data.frame(readr::read_tsv("RNAseq.quantifications.tpm.tsv"))
colnames(AnnoTPMs)[-1] <- paste0("RNAseq_", colnames(AnnoTPMs)[-1])

readr::write_tsv(cbind(Anno,
    AnnoCPMs[match(Anno$TSS, AnnoCPMs$TSS), -1],
    AnnoTPMs[match(Anno$NearestTranscript, AnnoTPMs$Transcript), -1]),
  "TSS.annotations.tsv")

Tx_detected <- Tx[Tx$tx_name %in% Anno$NearestTranscript]
names(Tx_detected) <- Tx_detected$tx_name
Tx_detected <- Tx_detected[as.character(strand(Tx_detected)) %in% c("+", "-")]
export(Tx_detected, "detected.transcripts.bed")

# TSS table used in the figures -------------------------------------------------

OldNames <- c(DS1 = "DS1", DS2 = "DS2", S241 = "I241", S242 = "I242",
  S721 = "I725", S722 = "I727", L61 = "6H1", L62 = "6H2", L261 = "26H1",
  L262 = "26H2", L571 = "57H1", L572 = "57H2", hen21 = "H1", hen22 = "H2",
  rrp41 = "R1", rrp42 = "R2")

# CAGE and GRO-cap at TSSs, GRO-seq over detected transcripts

system("bed2pos.pl tss.shrunk.bed > tmp.tss.pos")
system("bed2pos.pl detected.transcripts.bed > tmp.tx.pos")
system(paste("annotatePeaks.pl tmp.tss.pos tair10 -strand + -fragLength 1 -norm 1e6",
  "-d external_tagDirs/CAGE_14D_merged-tagDir external_tagDirs/GROcap-tagDir > grocap.norm.txt"))
system(paste("analyzeRepeats.pl tmp.tx.pos tair10 -strand + -tpm",
  "-d external_tagDirs/GROseq1-tagDir external_tagDirs/GROseq2-tagDir > groseq.tpm.txt"))
invisible(file.remove("tmp.tss.pos", "tmp.tx.pos"))

GROcap <- readr::read_tsv("grocap.norm.txt")[, c(1, 20:21)]
colnames(GROcap) <- c("TSS", "CAGE", "GROcap")
GROseq <- readr::read_tsv("groseq.tpm.txt")[, c(1, 9:10)]
colnames(GROseq) <- c("Transcript", "GROseq1", "GROseq2")

CPMs <- readr::read_tsv("csRNAseq.quantification.cpm.txt")
CPMs <- CPMs[, c("TSS", names(OldNames))]
colnames(CPMs)[-1] <- OldNames
TPMs <- readr::read_tsv("RNAseq.quantifications.tpm.tsv")
TPMs <- TPMs[, c("Transcript", names(OldNames))]
colnames(TPMs)[-1] <- paste0("RNAseq_", OldNames)

TSS_gro <- Anno
colnames(TSS_gro)[colnames(TSS_gro) == "NearestTranscript"] <- "NearestFeature"
colnames(TSS_gro)[colnames(TSS_gro) == "NearestTranscriptType"] <- "NearestFeatureType"
colnames(TSS_gro)[colnames(TSS_gro) == "NearestTranscriptSize"] <- "NearestFeatureSize"

TSS_gro <- dplyr::left_join(TSS_gro, CPMs, by = "TSS")
TSS_gro <- dplyr::left_join(TSS_gro, TPMs, by = c("NearestFeature" = "Transcript"))

readr::write_tsv(TSS_gro, "tss.annotated.tsv")
TSS_gro <- readr::read_tsv("tss.annotated.tsv")

TSS_gro <- dplyr::left_join(TSS_gro, GROcap, by = "TSS")
TSS_gro <- dplyr::left_join(TSS_gro, GROseq, by = c("NearestFeature" = "Transcript"))

TSS_gro <- TSS_gro[, c(1:25, 43, 42, 26:41, 44:45)]
TSS_gro <- TSS_gro[order(as.integer(substr(TSS_gro$TSS, 4, 4)),
  as.integer(gsub("TSS\\d_", "", TSS_gro$TSS))), ]

readr::write_tsv(TSS_gro, "tss.annotated.gro.tsv")

# ATAC-seq ---------------------------------------------------------------------

ACRs500 <- import("peaks_final.bed")
ACRs1 <- import("peaks_final_500bp.bed")

# ACRs are annotated using their centre position

ACRs1 <- resize(ACRs1, 1, "center")
start(ACRs1) <- start(ACRs1) + 1
end(ACRs1) <- start(ACRs1)

TSSs1 <- import("tss.size1.bed")

TEs <- A11[as.character(A11$type) %in% c("transposable_element", "transposable_element_gene")]

AnnoTSS <- Anno
rm(Anno)

Genes <- genes(TxDb)

TSSs1_prom <- promoters(TSSs1, 400, 100)
AllTx <- transcripts(TxDb)
MissingTSS <- resize(AllTx[!AllTx$tx_name %in% AnnoTSS$NearestTranscript], 1, "start")
colnames(mcols(MissingTSS)) <- c("id", "name")
MissingProm <- promoters(MissingTSS, 400, 100)

AllProms <- sort(c(TSSs1_prom, MissingProm[!overlapsAny(MissingProm, TSSs1_prom)]))

Genes_detected <- Genes[Genes$gene_id %in% AnnoTSS$NearestGene]

ACRsTSSsDistCS <- distanceToNearest(ACRs1, TSSs1)
ACRsTSSsDistAll <- distanceToNearest(ACRs1, c(TSSs1, MissingTSS))

# Start by assigning all peaks as intergenic, then re-assign as needed

Anno <- data.frame(row.names = NULL,
  PeakID = ACRs1$name,
  Type = "Intergenic",
  NearestTSS = NA_character_,
  NearestTSSDist = NA_integer_,
  NearestTSS_any = NA_character_,
  NearestTSSDist_any = NA_integer_,
  OverlappingProm = NA_character_,
  OverlappingProm_any = NA_character_,
  OverlappingGene = NA_character_,
  OverlappingTE = NA_character_
)

Anno$Type[overlapsAny(ACRs1, Genes)] <- "Intragenic"
Anno$Type[overlapsAny(ACRs1, TEs)] <- "TE"
Anno$Type[overlapsAny(ACRs1, c(TSSs1_prom, MissingProm))] <- "Promoter"

Anno$NearestTSS <- TSSs1$name[subjectHits(ACRsTSSsDistCS)]
Anno$NearestTSSDist <- mcols(ACRsTSSsDistCS)$distance
Anno$NearestTSS_any <- c(TSSs1, MissingTSS)$name[subjectHits(ACRsTSSsDistAll)]
Anno$NearestTSSDist_any <- mcols(ACRsTSSsDistAll)$distance

ACRsPromOvs <- findOverlaps(ACRs1, TSSs1_prom)
ACRsPromOvs <- tapply(TSSs1_prom$name[subjectHits(ACRsPromOvs)], queryHits(ACRsPromOvs),
  function(x) paste0(unique(x), collapse = "/"))
Anno$OverlappingProm[as.integer(names(ACRsPromOvs))] <- unname(ACRsPromOvs)

ACRsPromAllOvs <- findOverlaps(ACRs1, AllProms)
ACRsPromAllOvs <- tapply(AllProms$name[subjectHits(ACRsPromAllOvs)], queryHits(ACRsPromAllOvs),
  function(x) paste0(unique(x), collapse = "/"))
Anno$OverlappingProm_any[as.integer(names(ACRsPromAllOvs))] <- unname(ACRsPromAllOvs)

ACRsGeneOvs <- findOverlaps(ACRs1, Genes)
ACRsGeneOvs <- tapply(Genes$gene_id[subjectHits(ACRsGeneOvs)], queryHits(ACRsGeneOvs),
  function(x) paste0(unique(x), collapse = "/"))
Anno$OverlappingGene[as.integer(names(ACRsGeneOvs))] <- unname(ACRsGeneOvs)

ACRsTEOvs <- findOverlaps(ACRs1, TEs)
ACRsTEOvs <- tapply(TEs$ID[subjectHits(ACRsTEOvs)], queryHits(ACRsTEOvs),
  function(x) paste0(unique(x), collapse = "/"))
Anno$OverlappingTE[as.integer(names(ACRsTEOvs))] <- unname(ACRsTEOvs)

readr::write_tsv(Anno, "ACR.annotations.tsv")
