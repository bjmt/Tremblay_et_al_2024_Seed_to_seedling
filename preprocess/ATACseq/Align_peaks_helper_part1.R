# Benjamin J.M. Tremblay, June 2023
#
# Get ATAC-seq peak summits

library(rtracklayer)

b <- import("atacseq_bedGraphs/all_merged.bedGraph.gz")
p <- import("atacseq_peaks/peaks_unadjusted.bed")
names(p) <- p$name

get_peak_summit <- function(r) {
  pb <- subsetByOverlaps(b, r)
  pb <- pb[which.max(pb$score)]
  as.integer(mean(start(pb), end(pb)))
}

ps <- data.frame(
  Peak = p$name,
  seqnames = as.character(seqnames(p)),
  pos = NA_integer_
)

for (i in 1:nrow(ps)) {
  ps$pos[i] <- get_peak_summit(p[ps$Peak[i]])
}

pos2 <- GRanges(seqnames = ps$seqnames,
  IRanges(start = ps$pos), name = ps$Peak)
names(pos2) <- pos2$name

p <- import("atacseq_peaks/peaks_unadjusted.bed")
p$name <- paste0("ATAC", as.character(seqnames(p)), "_", pos$pos)
names(p) <- p$name

pos2$name <- p$name
names(pos2) <- p$name

export(pos2, "atacseq_peaks/peaks_final_summits.bed")

pos500 <- resize(pos2, 500, "center")

export(pos500, "atacseq_peaks/peaks_final_500bp.bed")

export(p, "atacseq_peaks/peaks_final.bed")

