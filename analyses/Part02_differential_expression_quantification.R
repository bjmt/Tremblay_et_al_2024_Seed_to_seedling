# Benjamin J.M. Tremblay, June 2023
#
# Differential expression, normalize quantification data

library(edgeR)
library(GenomicFeatures)
library(csaw)
library(rtracklayer)

# csRNA-seq --------------------------------------------------------------------

# Read in TSS counts per sample and switch TSS cluster names with new names

TSS_counts <- readr::read_tsv("all.tss.merged.counts.txt")
NameTab <- readr::read_tsv("tss.name_table.txt")
NameVec <- structure(NameTab$NewName, names = NameTab$OldName)
colnames(TSS_counts)[1] <- "TSS"
TSS_counts$TSS <- NameVec[TSS_counts$TSS]

m <- as.matrix(TSS_counts[, 20:ncol(TSS_counts)])
rownames(m) <- TSS_counts$TSS

# Note: Make sure the order of samples actually matches first, reorder if not!!!

colnames(m) <- c("L261", "L262", "L571", "L572", "L61", "L62", "DS1", "DS2",
  "hen21", "hen22", "S241", "S242", "S721", "S722", "rrp41", "rrp42")
m <- m[, c("DS1", "DS2", "S241", "S242", "S721", "S722",
  "L61", "L62", "L261", "L262", "L571", "L572", "hen21", "hen22", "rrp41", "rrp42")]
mode(m) <- "integer"

y <- DGEList(counts = m, group = gsub("[12]$", "", colnames(m)))
y <- calcNormFactors(y)

# Save TMM normalization factors for creating normalized bigWigs

write.table(y$samples, "csRNAseq.norm_factors.txt", sep = "\t", quote = FALSE)

m_cpm <- cpm(y)

# Normalized CPM tables

readr::write_tsv(cbind(TSS = rownames(m_cpm), data.frame(m_cpm, check.names = FALSE)),
  "csRNAseq.quantification.cpm.txt")

# Differential expression

d <- model.matrix(~0+y$samples$group)
y <- estimateDisp(y, d)
fit <- glmQLFit(y, d)

DS_S24 <- glmQLFTest(fit, contrast = c(0, 0, 0, -1, 0, 1, 0, 0))
S24_S72 <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, 0, -1, 1, 0))
S72_L6 <- glmQLFTest(fit, contrast = c(0, 0, 1, 0, 0, 0, -1, 0))
L6_L26 <- glmQLFTest(fit, contrast = c(1, 0, -1, 0, 0, 0, 0, 0))
L26_L57 <- glmQLFTest(fit, contrast = c(-1, 1, 0, 0, 0, 0, 0, 0))
L57_H <- glmQLFTest(fit, contrast = c(0, -1, 0, 0, 1, 0, 0, 0))
L57_R <- glmQLFTest(fit, contrast = c(0, -1, 0, 0, 0, 0, 0, 1))
H_R <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, -1, 0, 0, 1))
DS_L57 <- glmQLFTest(fit, contrast = c(0, 1, 0, -1, 0, 0, 0, 0))

DS_S24_res <- as.data.frame(topTags(DS_S24, Inf))
S24_S72_res <- as.data.frame(topTags(S24_S72, Inf))
S72_L6_res <- as.data.frame(topTags(S72_L6, Inf))
L6_L26_res <- as.data.frame(topTags(L6_L26, Inf))
L26_L57_res <- as.data.frame(topTags(L26_L57, Inf))
L57_H_res <- as.data.frame(topTags(L57_H, Inf))
L57_R_res <- as.data.frame(topTags(L57_R, Inf))
H_R_res <- as.data.frame(topTags(H_R, Inf))
DS_L57_res <- as.data.frame(topTags(DS_L57, Inf))

# Export results

write.table(DS_S24_res, "csRNAseq.DE.DS_S24.txt", quote = FALSE, sep = "\t")
write.table(S24_S72_res, "csRNAseq.DE.S24_S72.txt", quote = FALSE, sep = "\t")
write.table(S72_L6_res, "csRNAseq.DE.S72_L6.txt", quote = FALSE, sep = "\t")
write.table(L6_L26_res, "csRNAseq.DE.L6_L26.txt", quote = FALSE, sep = "\t")
write.table(L26_L57_res, "csRNAseq.DE.L26_L57.txt", quote = FALSE, sep = "\t")
write.table(L57_H_res, "csRNAseq.DE.L57_hen2.txt", quote = FALSE, sep = "\t")
write.table(L57_R_res, "csRNAseq.DE.L57_rrp4.txt", quote = FALSE, sep = "\t")
write.table(H_R_res, "csRNAseq.DE.hen2_rrp4.txt", quote = FALSE, sep = "\t")
write.table(DS_L57_res, "csRNAseq.DE.DS_L57.txt", quote = FALSE, sep = "\t")

# RNA-seq ----------------------------------------------------------------------

# Read in transcript level counts

Transcript_counts <- as.matrix(read.table("transcript_counts.tsv", row.names = 1,
  header = TRUE, check.names = FALSE))

# Note: Make sure the order of samples is correct, reorder if not!!

y <- DGEList(counts = Transcript_counts, group = gsub("[12]$", "", colnames(Transcript_counts)))
y <- y_transcripts[filterByExpr(y, group = y$samples$group), , keep.lib.size = FALSE]
y <- calcNormFactors(y)

# Save TMM normalization factors for creating normalized bigWigs

write.table(y$samples, "RNAseq.norm_factors.txt", quote = FALSE, sep = "\t")

# Perform differential expression of the seedling samples, necessary for exosome
# sensitivity analysis

d <- model.matrix(~0+y$samples$group)
y <- estimateDisp(y, d)
fit <- glmQLFit(y, d)

L57_H <- glmQLFTest(fit, contrast = c(0, -1, 0, 0, 1, 0, 0, 0))
L57_R <- glmQLFTest(fit, contrast = c(0, -1, 0, 0, 0, 0, 0, 1))
H_R <- glmQLFTest(fit, contrast = c(0, 0, 0, 0, -1, 0, 0, 1))

L57_H_res <- as.data.frame(topTags(L57_H, Inf))
L57_R_res <- as.data.frame(topTags(L57_R, Inf))
H_R_res <- as.data.frame(topTags(H_R, Inf))

# Export results

write.table(L57_H_res, "RNAseq.DE.L57_hen2.txt", quote = FALSE, sep = "\t")
write.table(L57_R_res, "RNAseq.DE.L57_rrp4.txt", quote = FALSE, sep = "\t")
write.table(H_R_res, "RNAseq.DE.hen2_rrp4.txt", quote = FALSE, sep = "\t")

# Now calculate normalized TPMs

A11 <- makeTxDbFromGFF("Araport11_novel.gff3")

lens_tx <- transcriptLengths(A11)
lens_tx <- lens_tx[!duplicated(lens_tx$tx_name), ]
lens_tx <- data.frame(row.names = lens_tx$tx_name, Length = lens_tx$tx_len)

y$genes <- lens_tx[rownames(y$counts), , drop = FALSE]

rpkm_transcripts <- rpkm(y_transcripts)
tpm_transcripts <- t(t(rpkm_transcripts) / colSums(rpkm_transcripts, na.rm = TRUE)) * 1e6
tpm_transcripts <- tpm_transcripts[complete.cases(tpm_transcripts), ]

readr::write_tsv(cbind("Transcript" = rownames(tpm_transcripts),
    data.frame(tpm_transcripts, check.names = FALSE)),
  "RNAseq.quantifications.tpm.tsv")

# ATAC-seq ---------------------------------------------------------------------

# Read in peaks and counts

ACRs <- import("peaks_final.bed")
names(ACRs) <- ACRs$name

ACR_counts <- as.matrix(read.table("all.counts.txt", header = TRUE,
  row.names = 1, check.names = FALSE))

# Note: Make sure the order of samples is correct, reorder if not!!

ACR_counts <- SummarizedExperiment(assays = SimpleList(counts = ACR_counts),
  rowRanges = ACRs, colData = DataFrame(totals = colSums(ACR_counts)))

y <- normOffsets(ACR_counts, se.out = TRUE)

m_cpm <- calculateCPM(y, use.offsets = TRUE, log = FALSE)

# Now export results

readr::write_tsv(cbind(PeakID = rownames(m_cpm),
  as.data.frame(m_cpm)), "ATACseq.quantification.rpm.txt")


