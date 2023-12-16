# Benjamin J.M. Tremblay, June 2023
#
# Refine TSS peaks to the minimum area containing reads and call summits

library(rtracklayer)

tss <- import("csrnaseq_TSSs_final/all.tss.merged.filtered.bed")

allPlus <- import("csrnaseq_bedGraphs/all.plus.bedGraph.gz")
allMinus <- import("csrnaseq_bedGraphs/all.minus.bedGraph.gz")

strand(allPlus) <- "+"
strand(allMinus) <- "-"
allBG <- sort(c(allPlus, allMinus))
allBG$score <- as.integer(allBG$score * 100)
allBG <- allBG[!is.na(allBG$score)]

quants <- c(0, 0.005, 0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 0.75 0.9, 0.95, 0.99, 0.995, 1)
allMat <- matrix(ncol = length(quants), nrow = length(tss))
rownames(allMat) <- tss$name
colnames(allMat) <- as.character(quants)
allMax <- structure(integer(length(tss)), names = tss$name)

for (i in 1:length(tss)) {
  i_x <- subsetByOverlaps(allBG, tss[i])
  if (length(i_x)) {
    # Get the single position coordinate with the max number of reads (primary TSS)
    allMax[i] <- start(i_x[which.max(i_x$score)])
    allMat[i, ] <- quantile(rep(start(i_x), times = i_x$score), p = quants)
  }
}

# TSS naming scheme:
# TSS[chromosome]_[primary TSS coordinate]
NewNames <- data.frame(
  OldName = tss$name,
  NewName = paste0("TSS", as.character(seqnames(tss)), "_", allMax)
)
readr::write_tsv(NewNames, "csrnaseq_TSSs_final/tss.name_table.txt")

allDF <- as.data.frame(cbind(TSS = NewNames$NewName, MaxTSS = allMax, as.data.frame(allMat)))
readr::write_tsv(allDF, "csrnaseq_TSSs_final/tss.size.stats.txt")

tss <- tss[, "score"]
names(tss) <- NewNames$NewName

# Refine the TSS peaks to only include the area with reads
start(tss) <- allMat[, 1, drop = TRUE]
end(tss) <- allMat[, ncol(allMat), drop = TRUE]
export(tss, "csrnaseq_TSSs_final/tss.shrunk.bed")

start(tss) <- allMax
end(tss) <- allMax
export(tss, "csrnaseq_TSSs_final/tss.size1.bed")


