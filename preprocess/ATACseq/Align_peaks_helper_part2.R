# Benjamin J.M. Tremblay, June 2023
#
# Get normalized bigWig tracks

library(rtracklayer)
library(csaw)

peaks <- import("atacseq_seqs/peaks_final.bed")

bg <- list(
  DS = import("atacseq_bedGraphs/DS.tiled10bp.bedGraph.gz"),
  L6 = import("atacseq_bedGraphs/L6.tiled10bp.bedGraph.gz"),
  L26 = import("atacseq_bedGraphs/L26.tiled10bp.bedGraph.gz"),
  L57 = import("atacseq_bedGraphs/L57.tiled10bp.bedGraph.gz")
)
bins <- import("genome_tiled10bp.bed")

bgm <- as.matrix(cbind(
  DS = bg[["DS"]]$score,
  L6 = bg[["L6"]]$score,
  L26 = bg[["L26"]]$score,
  L57 = bg[["L57"]]$score
))

bgm <- SummarizedExperiment(assays = SimpleList(counts = bgm),
  rowRanges = bins, colData = DataFrame(totals = colSums(bgm)))

bgm <- normOffsets(bgm, se.out = TRUE, weights = as.numeric(overlapsAny(bg[[1]], peaks)))

cpm <- calculateCPM(bgm, use.offsets = TRUE, log = FALSE)

for (i in 1:4) bg[[i]]$score <- cpm[i, ]

chroms <- structure(read.table("chrom.sizes")[[2]][1:5], names = 1:5)
for (i in 1:4) seqlengths(bg[[i]]) <- chroms

for (i in 1:4) export(bg[[i]], paste0("atacseq_bw/", names(bg)[i], ".bw"))

