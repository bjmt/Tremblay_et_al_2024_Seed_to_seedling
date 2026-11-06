# Note: the original code for this figure was lost, this is a recreation to reproduce the figure near exactly
#
# Figure 6d (ATAC-seq heatmap)
#
# Package requirements: rtracklayer, ComplexHeatmap, viridisLite, readr,
# matrixStats, genomation, zoo

library(rtracklayer)
library(ComplexHeatmap)

source("Functions.R")

dir.create("Figure6", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

BIDIR_PAIRS <- "Bidir_all_pairs.tsv"
ATAC_L57_BW <- "atacseq_bw/L57.bw"

# Divergent pairs, centred on the pcTSS and sorted by inter-TSS distance -------

bd_all <- readr::read_tsv(BIDIR_PAIRS)
div <- bd_all[bd_all$TypeV2=="nc-pc" & bd_all$dist <= 1000, ]

divBed <- with(div, GRanges(seq1,
  IRanges(ifelse(Type2=="protein_coding",pos2,pos1),ifelse(Type2=="protein_coding",pos2,pos1)),
  strand = ifelse(Type2=="protein_coding", "+", "-"),
  len = pos2 - pos1
))
strand(divBed)[divBed$len < 0] <- c("+" = "-", "-" = "+")[as.character(strand(divBed))[divBed$len < 0]]
divBed$len <- abs(divBed$len)
divBed <- divBed[order(divBed$len)]

# Figure 6d --------------------------------------------------------------------

atac_57 <- import(ATAC_L57_BW)

atac_div <- get_unstranded_signal(atac_57, promoters(divBed, 2000, 1000), keep=c(0.,.9),
  return.matrix=TRUE)

pdf("Figure6/Fig6d_div_atac.pdf", width = 2, height = 1.72)
draw(Heatmap(atac_div,
  show_row_names = FALSE, show_column_names = FALSE, use_raster=T,raster_quality=3,
  col=viridisLite::plasma(100,begin=0),
  cluster_rows = FALSE, cluster_columns= FALSE
))
dev.off()
