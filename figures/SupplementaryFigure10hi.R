# Note: the original code for this figure was lost, this is a recreation to reproduce the figure near exactly
#
# Supplementary Figure 10h, i
#
# Package requirements: GenomicRanges, readr

library(GenomicRanges)

dir.create("SupplementaryFigure10", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

# Tested candidate enhancers (Zhu et al. 2015, Yan et al. 2019)
ZHU2015_BED <- "Zhu2015_validated_enhancers.bed"
YAN2019_BED <- "Yan2019_validated_enhancers.bed"

ENHANCER_ACTIVITY <- "Final_enhancer_activity.tsv"
TSS_TABLE <- "tss.annotated.gro.tsv"

# Setup ------------------------------------------------------------------------

read_candidates <- function(bed) {
  x <- as.data.frame(readr::read_tsv(bed, col_names = FALSE,
    col_types = readr::cols(X2 = "i", X3 = "i", .default = "c")))
  x <- x[!duplicated(paste(x$X1, x$X2, x$X3)), ]
  GRanges(x$X1, IRanges(x$X2 + 1, x$X3), name = gsub("_", " ", x$X4))
}

enh <- as.data.frame(readr::read_tsv(ENHANCER_ACTIVITY))
enhGR <- GRanges(enh$seqnames, IRanges(enh$start, enh$end), name = enh$name)
enhActivity <- as.matrix(enh[, c("DS", "S24", "S72", "L6", "L26", "L57")])
rownames(enhActivity) <- enh$name

tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tssA <- tssA[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding", ]
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

detected_enhancers <- function(candidates) {
  hits <- findOverlaps(candidates, enhGR)
  as.character(tapply(enhGR$name[subjectHits(hits)],
    factor(queryHits(hits), levels = seq_along(candidates)),
    function(x) paste0(x, collapse = ",")))
}

enhancer_gene_pcc <- function(enhancers, gene) {
  tss <- tssA$TSS[tssA$NearestGene == gene]
  if (is.na(enhancers) || !length(tss)) return(NA_real_)
  enhancers <- strsplit(enhancers, ",", fixed = TRUE)[[1]]
  max(cor(t(log2(1 + tssMm[tss, , drop = FALSE])),
    t(log2(1 + enhActivity[enhancers, , drop = FALSE]))))
}

# Supplementary Figure 10h -----------------------------------------------------

zhu <- read_candidates(ZHU2015_BED)

zhuTable <- data.frame(
  Name = paste0("(", seq_along(zhu), ") ", zhu$name),
  GUS_positive = !zhu$name %in% c("Common 3", "Common 6", "Common 7", "Flower 3"),
  Enhancer_detected = overlapsAny(zhu, enhGR),
  Enhancer = detected_enhancers(zhu)
)
zhuTable

readr::write_tsv(zhuTable, "SupplementaryFigure10/SuppFig10h_zhu2015_enhancers.tsv", na = "")

# Supplementary Figure 10i -----------------------------------------------------

yan <- read_candidates(YAN2019_BED)

yanTargets <- c(
  "AP1" = "AT1G69120", "CRC" = "AT1G69180", "JAG" = "AT1G68480",
  "RAV1" = "AT1G13260", "RGA2" = "AT1G14920", "RGL1" = "AT1G66350",
  "RSB8/FT" = "AT1G65480", "SWEET1" = "AT1G21460", "UFO" = "AT1G30950",
  "ATC" = "AT2G27550", "BOP2" = "AT2G41370", "ETT/ARF3" = "AT2G33860",
  "MYB91/AS1" = "AT2G37630", "PIL1" = "AT2G46970", "RSB1/AGL6" = "AT2G45650",
  "SOC1/AGL20" = "AT2G45660", "SPL9" = "AT2G42200", "WOX3" = "AT2G28610",
  "BOP1" = "AT3G57130", "SHP1" = "AT3G58780", "TCP18" = "AT3G18550",
  "HXK1" = "AT4G29130", "MED12" = "AT4G00450", "SPT" = "AT4G36930"
)

yanTable <- data.frame(
  Name = paste0("(", seq_along(yan), ") ", yan$name),
  Gene = unname(yanTargets[yan$name]),
  GUS_positive = !yan$name %in% c("ATC", "MED12"),
  Enhancer_detected = overlapsAny(yan, enhGR),
  Enhancer = detected_enhancers(yan)
)
yanTable$PCC <- unname(mapply(enhancer_gene_pcc, yanTable$Enhancer, yanTable$Gene))
yanTable$PCC_above_0.5 <- !is.na(yanTable$PCC) & yanTable$PCC > 0.5
yanTable

readr::write_tsv(yanTable, "SupplementaryFigure10/SuppFig10i_yan2019_enhancers.tsv", na = "")
