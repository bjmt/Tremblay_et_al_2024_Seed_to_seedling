# Note: the original code for this figure was lost, this is a recreation to reproduce the figure near exactly
#
# Supplementary Figure 14
#
# Package requirements: rtracklayer, readr

library(rtracklayer)

dir.create("SupplementaryFigure14", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_SHRUNK_BED <- "tss.shrunk.bed"
DETECTED_TX_BED <- "detected.transcripts.bed"

# csRNA-seq and sRNA-seq (input) tag directories
CSRNA_TAGDIRS <- "csrnaseq_tagDirs"
SRNA_TAGDIRS <- "srnaseq_tagDirs"

# Setup ------------------------------------------------------------------------

samples <- c("DS1", "DS2", "S241", "S242", "S721", "S722", "L61", "L62",
  "L261", "L262", "L571", "L572", "rrp41", "rrp42", "hen21", "hen22")

tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tssA <- tssA[!is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding", ]

tssW <- import(TSS_SHRUNK_BED)
detTx <- import(DETECTED_TX_BED)

# One TSS region and one gene region per protein coding transcript (n = 19,688)
pcTx <- unique(tssA$NearestFeature)
pcTSS <- structure(tssA$TSS, names = tssA$NearestFeature)[pcTx]
geneRegions <- detTx[match(pcTx, detTx$name)]
tssRegions <- tssW[match(pcTSS, tssW$name)]

read_tags <- function(tagDir) {
  tags <- do.call(rbind, lapply(list.files(tagDir, pattern = "[.]tags[.]tsv$",
    full.names = TRUE), readr::read_tsv, col_names = FALSE, col_types = "-cicd-"))
  GRanges(tags$X2, IRanges(tags$X3, width = 1),
    strand = ifelse(tags$X4 == 0, "+", "-"), count = tags$X5)
}

total_tags <- function(tagDir) {
  as.numeric(strsplit(readLines(file.path(tagDir, "tagInfo.txt"))[2], "\t")[[1]][3])
}

count_tags <- function(regions, tags) {
  hits <- findOverlaps(regions, tags)
  out <- numeric(length(regions))
  counts <- tapply(tags$count[subjectHits(hits)], queryHits(hits), sum)
  out[as.integer(names(counts))] <- counts
  out
}

count_library <- function(tagDirs) {
  geneCounts <- matrix(0, nrow = length(pcTx), ncol = length(samples),
    dimnames = list(pcTx, samples))
  tssCounts <- geneCounts
  totals <- structure(numeric(length(samples)), names = samples)
  for (i in samples) {
    tagDir <- file.path(tagDirs, paste0(i, "-tagDir"))
    tags <- read_tags(tagDir)
    geneCounts[, i] <- count_tags(geneRegions, tags)
    tssCounts[, i] <- count_tags(tssRegions, tags)
    totals[i] <- total_tags(tagDir)
  }
  list(gene = geneCounts, tss = tssCounts, total = totals)
}

csrna <- count_library(CSRNA_TAGDIRS)
srna <- count_library(SRNA_TAGDIRS)

csrna_ratio <- (csrna$tss + 1) / (csrna$gene + 1)
srna_ratio <- (srna$tss + 1) / (srna$gene + 1)

csrna_frac <- rbind(Gene = colSums(csrna$gene), TSS = colSums(csrna$tss)) /
  rbind(csrna$total, csrna$total)
srna_frac <- rbind(Gene = colSums(srna$gene), TSS = colSums(srna$tss)) /
  rbind(srna$total, srna$total)

round(csrna_frac, 4)
round(srna_frac, 4)

# Replicate pairs, with labels alternating between two lines
sample_axis <- function(x) {
  groups <- c("DS", "S24", "S72", "L6", "L26", "L57", "rrp4-2", "hen2-4")
  x1 <- x[seq(1, length(x), 2)]
  x2 <- x[seq(2, length(x), 2)]
  line <- ifelse(seq_along(groups) %% 2, 3, 2)
  y_line <- function(l) grconvertY(grconvertY(par("usr")[3], "user", "inches") -
    l * par("csi"), "inches", "user")
  rect(x1, y_line(0.6), x2, y_line(0), xpd = NA)
  segments((x1 + x2) / 2, y_line(0.6), (x1 + x2) / 2, y_line(line + 0.1), xpd = NA)
  mtext(groups, side = 1, at = (x1 + x2) / 2, line = line + 0.2,
    font = ifelse(groups %in% c("rrp4-2", "hen2-4"), 3, 1))
}

# Supplementary Figure 14a -----------------------------------------------------

pdf("SupplementaryFigure14/SuppFig14a_csrna_tss_gene_ratio.pdf", width = 3.4, height = 4)
par(mar = c(5, 4, 4, 1) + 0.1)
boxplot(csrna_ratio, outline = FALSE, col = "lightgray", ylim = c(0, 2),
  xaxt = "n", main = "csRNA libraries", ylab = "TSS/Gene read ratio")
sample_axis(seq_along(samples))
dev.off()

# Supplementary Figure 14b -----------------------------------------------------

pdf("SupplementaryFigure14/SuppFig14b_input_tss_gene_ratio.pdf", width = 3.4, height = 4)
par(mar = c(5, 4, 4, 1) + 0.1)
boxplot(srna_ratio, outline = FALSE, col = "lightgray", ylim = c(0, 2),
  xaxt = "n", main = "Input libraries", ylab = "TSS/Gene read ratio")
sample_axis(seq_along(samples))
dev.off()

# Supplementary Figure 14c -----------------------------------------------------

pdf("SupplementaryFigure14/SuppFig14c_csrna_read_fractions.pdf", width = 3.8, height = 4)
par(mar = c(5, 4, 4, 1) + 0.1)
bp <- barplot(csrna_frac, beside = TRUE, space = c(0, 0.2), col = c("darkgray", "black"),
  border = NA, ylim = c(0, 1), yaxs = "i", axes = FALSE,
  names.arg = rep("", length(samples)), main = "csRNA", ylab = "Fraction of total reads")
axis(2, seq(0, 1, 0.25), las = 1)
box(lwd = 2)
sample_axis(colMeans(bp))
dev.off()

# Supplementary Figure 14d -----------------------------------------------------

pdf("SupplementaryFigure14/SuppFig14d_input_read_fractions.pdf", width = 3.8, height = 4)
par(mar = c(5, 4, 4, 1) + 0.1)
bp <- barplot(srna_frac, beside = TRUE, space = c(0, 0.2), col = c("darkgray", "black"),
  border = NA, ylim = c(0, 1), yaxs = "i", axes = FALSE,
  names.arg = rep("", length(samples)), main = "Input", ylab = "Fraction of total reads")
axis(2, seq(0, 1, 0.25), las = 1)
box(lwd = 2)
legend("topright", legend = c("Gene", "TSS"), fill = c("darkgray", "black"),
  border = NA, bty = "n")
sample_axis(colMeans(bp))
dev.off()
