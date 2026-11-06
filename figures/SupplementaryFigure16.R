# Note: the original code for this figure was lost, this is a recreation to reproduce the figure near exactly
#
# Supplementary Figure 16
#
# Package requirements: ggplot2, readr

library(ggplot2)

dir.create("SupplementaryFigure16", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

CLUSTERS_TABLE <- "csRNAseq_clustering/Clusters_anno.tsv"

# Supplementary Figure 16 ------------------------------------------------------

tssA <- readr::read_tsv(CLUSTERS_TABLE)
tssA <- tssA[tssA$Cluster %in% 1:6, ]

unstable <- as.data.frame(table(Cluster = paste0("C", tssA$Cluster[is.na(tssA$NearestGeneType)])))
unstable

ggplot(unstable, aes(Cluster, Freq)) +
  geom_col() +
  xlab(NULL) +
  ylab("# of unstable TSS") +
  theme_bw() +
  theme(
    axis.text = element_text(colour = "black"),
    panel.border = element_rect(colour = "black")
  )
ggsave("SupplementaryFigure16/SuppFig16_unstable_tss_per_cluster.pdf", width = 2.6, height = 1.9)
