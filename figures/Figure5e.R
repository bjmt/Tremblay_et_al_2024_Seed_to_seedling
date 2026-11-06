# Note: the original code for this figure was lost, this is a recreation to reproduce the figure near exactly
#
# Figure 5e
#
# Package requirements: rtracklayer, GenomicFeatures, ComplexHeatmap,
# viridisLite, readr, matrixStats, genomation, zoo

library(rtracklayer)
library(GenomicFeatures)
library(ComplexHeatmap)

source("Functions.R")

dir.create("Figure5", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

ARAPORT11_GFF <- "Araport11_novel.gff3"
TSS_TABLE <- "tss.annotated.gro.tsv"
TSS_BED <- "tss.size1.bed"
ANTISENSES <- "Antisenses.txt"
MNASE_BW <- "Ath_leaf_NPS.bw"

# Antisense pairs --------------------------------------------------------------

a11 <- makeTxDbFromGFF(ARAPORT11_GFF)
tssA <- as.data.frame(readr::read_tsv(TSS_TABLE))
tss1 <- import(TSS_BED)

pair_dist <- function(x, y) {
  distance(tss1[tss1$name==x], tss1[tss1$name==y], ignore.strand=TRUE)
}

a11_g <- genes(a11)
g2len <- structure(width(a11_g), names = a11_g$gene_id)
final_as <- readr::read_tsv(ANTISENSES)
final_as$SenseLen <- g2len[final_as$PcGene]
final_as$TSS_dist <- NA_integer_
for (i in 1:nrow(final_as)) {
  final_as$TSS_dist[i] <- pair_dist(final_as$PcTSS[i], final_as$AntisenseTSS[i])
}
final_as$TTS_dist <- final_as$SenseLen - final_as$TSS_dist
final_as$Dist_ratio <- with(final_as, pmin(1,(SenseLen - TTS_dist) / SenseLen))
final_as$Proximal <- final_as$Dist_ratio<=0.5 & final_as$TTS_dist>1000

tss1_pc <- tss1[tss1$name %in% final_as$PcTSS]
tss1_pc$Proximal <- tss1_pc$name %in% final_as$PcTSS[final_as$Proximal]

non_as_pc <- tssA$TSS[!is.na(tssA$NearestGene) & tssA$NearestGeneType == "protein_coding"]
non_as_pc <- non_as_pc[!non_as_pc %in% tss1_pc$name]

tss1_a <- tss1[tss1$name %in% non_as_pc]
tss1_b <- tss1[tss1$name %in% tss1_pc$name[tss1_pc$Proximal]]
tss1_c <- tss1[tss1$name %in% tss1_pc$name[!tss1_pc$Proximal]]

# Figure 5e --------------------------------------------------------------------

mnase3 <- import(MNASE_BW)
seqlevels(mnase3) <- as.character(1:5)

tss2len <- structure(tssA$NearestGene, names = tssA$TSS)
tss2len <- tss2len[!is.na(tss2len)]
tss2len <- tss2len[unname(tss2len) %in% names(g2len)]
tss2len <- structure(unname(g2len[unname(tss2len)]), names = names(tss2len))

tss2dist <- structure(final_as$TSS_dist, names = final_as$PcTSS)

# No antisense (sorted by gene length), proximal and distal antisense (sorted by distance)

tss1_a2 <- tss1_a[sample(1:length(tss1_a), 4000)]
no_as_win <- promoters(tss1_a2[order(tss2len[tss1_a2$name])], 500, 4000)
prox_as_win <- promoters(tss1_b[order(tss2dist[tss1_b$name])], 500, 4000)
dist_as_win <- promoters(tss1_c[order(tss2dist[tss1_c$name])], 500, 4000)

no_as_map <- get_unstranded_signal(mnase3, no_as_win, keep = c(.00, .95), return.matrix = TRUE)
prox_as_map <- get_unstranded_signal(mnase3, prox_as_win, keep = c(.00, .95), return.matrix = TRUE)
dist_as_map <- get_unstranded_signal(mnase3, dist_as_win, keep = c(.00, .95), return.matrix = TRUE)

mnase_maps <- list(
  Fig5e_heatmap_mnase3_no_antisense_4k = list(m = no_as_map, h = 2.5),
  Fig5e_heatmap_mnase3_proxi_antisense = list(m = prox_as_map, h = 0.71),
  Fig5e_heatmap_mnase3_distal_antisense = list(m = dist_as_map, h = 1.72)
)
for (i in names(mnase_maps)) {
  pdf(paste0("Figure5/", i, ".pdf"), width = 2, height = mnase_maps[[i]]$h)
  draw(Heatmap(mnase_maps[[i]]$m,
    col=viridisLite::plasma(100,begin=0),
    show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
    cluster_rows = FALSE, cluster_columns= FALSE
  ))
  dev.off()
}

# Gene lengths

pdf("Figure5/Fig5e_heatmap_mnase3_distal_antisense_lengths.pdf", width = 4, height = 3)
plot(1:length(tss1_c), tss2len[tss1_c$name][order(tss2dist[tss1_c$name])],
  pch = 19, cex = .4, type = "h", ylim = c(0, 20000)
)
dev.off()
