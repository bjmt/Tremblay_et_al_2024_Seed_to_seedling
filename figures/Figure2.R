# Benjamin J.M. Tremblay
#
# Figure 2 and Supplementary Figure 2
#
# Package requirements: ComplexHeatmap, simplifyEnrichment, ggplot2, circlize,
# colorBlindness, viridisLite, awtools, readr, rtracklayer, dplyr, cowplot,
# gtools

library(ComplexHeatmap)
library(simplifyEnrichment)
library(ggplot2)

dir.create("Figure2", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

TSS_CLUSTERS <- "csRNAseq_clustering/Clusters_anno.tsv"
TSS_CLUSTER_GENES_DIR <- "csRNAseq_clustering/clusters-final-genes"
TSS_CLUSTER_GO_DIR <- "csRNAseq_clustering/clusters-final-genes-go"
ATAC_CLUSTERS <- "ATACseq_clustering3/Clusters_anno.tsv"
ATAC_CLUSTER_GENES_DIR <- "ATACseq_clustering3/clusters-final-genes"
ATAC_CLUSTER_GO_DIR <- "ATACseq_clustering3/clusters-final-genes-go"
RNA_TPM <- "RNAseq.quantifications.tpm.tsv"
RNA_CLUSTER_TRANSCRIPTS_DIR <- "RNAseq_clustering/clusters-final-transcripts"
RNA_CLUSTER_GO_DIR <- "RNAseq_clustering/clusters-final-genes-go"
ARAPORT11_GFF <- "Araport11_novel.gff3"
TF_LIST <- "PlantTFDB_Ath_TF_list.txt"

# Load clusters ----------------------------------------------------------------

avg_cls <- function(m, cls) {
  colMeans(m[cls, ])
}

atacA <- readr::read_tsv(ATAC_CLUSTERS)
atacA <- atacA[atacA$Cluster != 0, ]
atacA <- atacA[sample.int(nrow(atacA)), ]
atacA <- atacA[order(atacA$Cluster), ]
atacM <- as.matrix(atacA[, 13:20])
rownames(atacM) <- atacA$PeakID
atacMz <- t(scale(t(atacM)))

atacMz_m <- t(sapply(1:5, function(x) avg_cls(atacMz, atacA$Cluster == x)))
rownames(atacMz_m) <- paste0("A", 1:5)

atac_g <- lapply(list.files(ATAC_CLUSTER_GENES_DIR, full = TRUE),
  readr::read_lines)
names(atac_g) <- list.files(ATAC_CLUSTER_GENES_DIR)
names(atac_g) <- gsub("[.]txt$", "", gsub("Cluster_", "", names(atac_g), fixed = TRUE))

atac_go <- lapply(list.files(ATAC_CLUSTER_GO_DIR, full = TRUE),
  readr::read_tsv)
names(atac_go) <- list.files(ATAC_CLUSTER_GO_DIR)
names(atac_go) <- gsub("[.]txt$", "", gsub("Cluster_", "", names(atac_go), fixed = TRUE))
for (i in seq_along(atac_go)) atac_go[[i]] <- atac_go[[i]][atac_go[[i]]$source == "GO:BP", ]

tssA <- readr::read_tsv(TSS_CLUSTERS)
tssA <- tssA[tssA$Cluster != 0, ]
tssA <- tssA[sample.int(nrow(tssA)), ]
tssA <- tssA[order(tssA$Cluster), ]
tssM <- as.matrix(tssA[, 11:22])
rownames(tssM) <- tssA$TSS
tssMz <- t(scale(t(tssM)))

tssMz_m <- t(sapply(1:6, function(x) avg_cls(tssMz, tssA$Cluster == x)))
rownames(tssMz_m) <- paste0("C", 1:6)

tss_g <- lapply(list.files(TSS_CLUSTER_GENES_DIR, full = TRUE),
  readr::read_lines)
names(tss_g) <- list.files(TSS_CLUSTER_GENES_DIR)
names(tss_g) <- gsub("[.]txt$", "", gsub("Genes_", "", names(tss_g), fixed = TRUE))

tss_go <- lapply(
  list.files(TSS_CLUSTER_GO_DIR, full = TRUE),
  readr::read_tsv)
names(tss_go) <- list.files(TSS_CLUSTER_GO_DIR)
names(tss_go) <- gsub("[.]txt$", "", gsub("Cluster_", "", names(tss_go), fixed = TRUE))
for (i in seq_along(tss_go)) tss_go[[i]] <- tss_go[[i]][tss_go[[i]]$source == "GO:BP", ]

rnaM <- as.matrix(read.table(RNA_TPM, row.names = 1, header = TRUE,
    check.names = FALSE))[, 1:12]
rnaM <- rnaM[apply(rnaM, 1, function(x) sum(x>=2)>=2), ]
rnaMz <- t(scale(t(rnaM)))
rna_t <- lapply(list.files(RNA_CLUSTER_TRANSCRIPTS_DIR, full = TRUE), readr::read_lines)
names(rna_t) <- paste0("R", 1:6)
rna_g <- lapply(rna_t, function(x) unique(gsub("[.]\\d+", "", x)))

rna2cls <- c(
  structure(rna_t[[1]], names = rep("R1", length(rna_t[[1]]))),
  structure(rna_t[[2]], names = rep("R2", length(rna_t[[2]]))),
  structure(rna_t[[3]], names = rep("R3", length(rna_t[[3]]))),
  structure(rna_t[[4]], names = rep("R4", length(rna_t[[4]]))),
  structure(rna_t[[5]], names = rep("R5", length(rna_t[[5]]))),
  structure(rna_t[[6]], names = rep("R6", length(rna_t[[6]])))
)

rnaMz_m <- t(sapply(rna_t, function(x) colMeans(rnaMz[x, ])))
rownames(rnaMz_m) <- paste0("R", 1:6)

rna_go <- lapply(list.files(RNA_CLUSTER_GO_DIR, full = TRUE), readr::read_tsv)
names(rna_go) <- paste0("R", 1:6)
for (i in seq_along(rna_go)) rna_go[[i]] <- rna_go[[i]][rna_go[[i]]$source == "GO:BP", ]

a11 <- rtracklayer::import(ARAPORT11_GFF)

# Figure 2a --------------------------------------------------------------------

# GO keyword enrichment of each cluster
tss_go_wc <- anno_word_cloud_from_GO(
  paste0("C", tssA$Cluster),
  lapply(tss_go[2:7], function(x) x$term_id[x$p_value < 0.01]),
  fontsize_range = c(2, 8),
  bg_gp = gpar(col = "#000000", fill = "#ffffff"),
  word_cloud_grob_param = list(max_width = unit(30, "mm"), word_space = unit(2, "pt"),
    line_space = unit(2, "pt")),
  stat = "pvalue", max_words = 15, min_stat = 1,
  exclude_words = c("process", "metabolic", "response", "cellular", "cell",
    "compound", "biological", "macromolecule", "levels", "decreased", "organic"))

pdf("Figure2/Fig2a_csrna_heatmap.pdf", height = 2.4, width = 3.9)
Heatmap(tssMz, row_split = paste0("C", tssA$Cluster),
  col = circlize::colorRamp2(c(-1.9, 0, 1.9), c("blue", "white", "red")),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  show_column_names = FALSE, row_title_rot = 0, row_title_gp = gpar(fontsize = 6),
  cluster_column_slices = FALSE, cluster_row_slices = FALSE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  column_gap = unit(0.2, "mm"), row_gap = unit(0.6, "mm"),
  raster_quality = 3,
  right_annotation = rowAnnotation(wc = tss_go_wc))
dev.off()

pdf("Figure2/Fig2a_csrna_heatmap_legend.pdf", height = 0.5, width = 1.5)
draw(
  Legend(col_fun = circlize::colorRamp2(c(-1.9, 0, 1.9), c("blue", "white", "red")),
    title = "Z-score", direction = "horizontal", title_position = "lefttop",
    at = c(-2, 0, 2), labels_gp = gpar(fontsize = 5),
    legend_width = unit(2, "cm"), legend_height = unit(1, "mm"),
    title_gp = gpar(fontsize = 5, fontface = "bold"))
)
dev.off()

# Figure 2b --------------------------------------------------------------------

atac_go_wc <- anno_word_cloud_from_GO(
  paste0("A", atacA$Cluster),
  lapply(atac_go[2:6], function(x) x$term_id[x$p_value < 0.05]),
  fontsize_range = c(2, 8),
  bg_gp = gpar(col = "#000000", fill = "#ffffff"),
  word_cloud_grob_param = list(max_width = unit(30, "mm"), word_space = unit(2, "pt"),
    line_space = unit(2, "pt")),
  stat = "pvalue", max_words = 15, min_stat = 1,
  exclude_words = c("process", "metabolic", "response", "cellular", "cell",
    "compound", "biological", "macromolecule", "levels", "decreased", "organic"))

pdf("Figure2/Fig2b_atac_heatmap.pdf", height = 1.68, width = 3.46)
Heatmap(atacMz, row_split = paste0("A", atacA$Cluster),
  col = circlize::colorRamp2(c(-1.9, 0, 1.9), c("blue", "white", "red")),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  show_column_names = FALSE, row_title_rot = 0, row_title_gp = gpar(fontsize = 6),
  cluster_column_slices = FALSE, cluster_row_slices = FALSE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  column_gap = unit(0.2, "mm"), row_gap = unit(0.6, "mm"),
  raster_quality = 3, use_raster = TRUE,
  right_annotation = rowAnnotation(wc = atac_go_wc))
dev.off()

pdf("Figure2/Fig2b_atac_heatmap_legend.pdf", height = 0.5, width = 1.5)
draw(
  Legend(col_fun = circlize::colorRamp2(c(-1.9, 0, 1.9), c("blue", "white", "red")),
    title = "Z-score", direction = "horizontal", title_position = "lefttop",
    at = c(-2, 0, 2), labels_gp = gpar(fontsize = 5),
    legend_width = unit(1, "cm"), legend_height = unit(1, "mm"),
    title_gp = gpar(fontsize = 5, fontface = "bold"))
)
dev.off()

# Figure 2c --------------------------------------------------------------------

# PCC over the time-points shared with ATAC-seq
cCor <- cor(t(atacMz_m), t(tssMz_m[, -(3:6)]))

# Jaccard coefficient and one-sided Fisher's exact test
ovCoef <- matrix(0, ncol = 6, nrow = 5, dimnames = list(paste0("A", 1:5), paste0("C", 1:6)))
ovFish <- matrix(0, ncol = 6, nrow = 5, dimnames = list(paste0("A", 1:5), paste0("C", 1:6)))
for (i in 1:nrow(ovCoef)) {
  for (j in 1:ncol(ovCoef)) {
    ovCoef[i, j] <- sum(atac_g[[i+1]] %in% tss_g[[j+1]]) /
      length(unique(c(atac_g[[i+1]], tss_g[[j+1]])))
    ovFish[i, j] <- fisher.test(alt = "g",
      matrix(c(
          sum(atac_g[[i+1]] %in% tss_g[[j+1]]),
          sum(!atac_g[[i+1]] %in% tss_g[[j+1]]),
          sum(!tss_g[[j+1]] %in% atac_g[[i+1]]),
          length(unique(c(unlist(tss_g, use.names=FALSE), unlist(atac_g, use.names=FALSE)))) -
            length(unique(c(atac_g[[i+1]], tss_g[[j+1]])))
        ), ncol = 2))$p.value
  }
}
print(ovFish < 1e-6)

pdf("Figure2/Fig2c_comparing_clusters.pdf", width = 3.525, height = 1.3)
Heatmap(cCor,
  cluster_columns = FALSE, cluster_rows = FALSE, name = "PCC",
  col = circlize::colorRamp2(c(-1, 0, 1),
    c(colorBlindness::safeColors[c("vermillion", "skyBlue")], "white")[c(2, 3, 1)]),
  column_names_rot = 0, column_names_centered = TRUE,
  row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), at = c(-1, 0, 1),
    title_gp = gpar(fontsize = 5), legend_height = unit(8, "mm")),
  column_title = "csRNA-seq clusters", row_title = "ATAC-seq clusters") +
Heatmap(ovCoef,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = viridisLite::inferno(100, begin = 0.00, end = 1),
  row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), at = c(0, 0.2),
    title_gp = gpar(fontsize = 5), legend_height = unit(8, "mm")),
  name = "Overlap Coef.", column_names_rot = 0, column_names_centered = TRUE
)
dev.off()

# Figure 2d --------------------------------------------------------------------
tssA$Type2 <- "Protein-coding"
tssA$Type2[is.na(tssA$NearestGeneType) | tssA$NearestGeneType != "protein_coding"] <- "Non-coding"
tssA$Type2[is.na(tssA$NearestGeneType) | tssA$NearestGeneType == "putative_rna"] <- "Unknown"
tssA_types <- table(paste0("C", tssA$Cluster), tssA$Type2)
for (i in 1:nrow(tssA_types)) tssA_types[i, ] <- tssA_types[i, ] / sum(tssA_types[i, ])
tssA_types <- as.data.frame(tssA_types)
tssA_types$Var2 <- factor(as.character(tssA_types$Var2),
  levels = c("Protein-coding", "Non-coding", "Unknown"))

p2 <- ggplot(tssA_types, aes(Var1, Freq, group = Var2, fill = Var2)) +
  geom_bar(stat = "identity", width = .8) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("csRNA-seq cluster") +
  ylab("Proportion") +
  scale_fill_manual(values = rev(awtools::ppalette[-(2:5)])[c(1, 4, 3, 2)],
    name = element_blank()) +
  theme_bw() +
  theme(
    aspect.ratio = 0.6667,
    panel.border = element_rect(colour = "black", size = 1),
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black", size = unit(6, "pt")),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    legend.key.size = unit(2, "mm"),
    axis.text = element_text(colour = "black", size = unit(6, "pt"))
  )

ggsave("Figure2/Fig2d_csRNA_cls_types.pdf", p2, width = 3.5, height = 1.1)

# Figure 2e --------------------------------------------------------------------

atacA_types <- table(paste0("A", atacA$Cluster), atacA$Type)
for (i in 1:nrow(atacA_types)) atacA_types[i, ] <- atacA_types[i, ] / sum(atacA_types[i, ])
atacA_types <- as.data.frame(atacA_types)

p1 <- ggplot(atacA_types, aes(Var1, Freq, group = Var2, fill = Var2)) +
  geom_bar(stat = "identity", width = .8) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("ATAC-seq cluster") +
  ylab("Proportion") +
  scale_fill_manual(values = awtools::ppalette[-1], name = element_blank()) +
  theme_bw() +
  theme(
    aspect.ratio = 0.8,
    panel.border = element_rect(colour = "black", size = 1),
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black", size = unit(6, "pt")),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    legend.key.size = unit(2, "mm"),
    axis.text = element_text(colour = "black", size = unit(6, "pt"))
  )

ggsave("Figure2/Fig2e_ATAC_cls_types.pdf", p1, width = 3.5, height = 1.1)

# Supplementary Figure 2a ------------------------------------------------------
tssA_ispc <- !is.na(tssA$NearestGeneType) & tssA$NearestGeneType == "protein_coding"
tssA_annnc <- !is.na(tssA$NearestGeneType) & !tssA$NearestGeneType %in% c(
  "putative_rna", "protein_coding")
tssA_novnc <- is.na(tssA$NearestGeneType) | tssA$NearestGeneType == "putative_rna"
pct_pcnc <- data.frame(row.names = NULL,
  Sample = rep(colnames(tssA)[11:22], each = 3),
  Timepoint = rep(c("DS", "S24", "S72", "L6", "L26", "L57"), each = 6),
  Type = rep(c("Protein-coding", "Non-coding", "Unknown"), 12),
  Percent = NA_real_
)
for (i in 1:12) {
  i_pct1 <- 100 * (tapply(tssA[[10 + i]], tssA_ispc, sum) / sum(tssA[[10 + i]]))
  i_pct2 <- 100 * (tapply(tssA[[10 + i]], tssA_annnc, sum) / sum(tssA[[10 + i]]))
  i_pct3 <- 100 * (tapply(tssA[[10 + i]], tssA_novnc, sum) / sum(tssA[[10 + i]]))
  pct_pcnc$Percent[pct_pcnc$Sample == colnames(tssA)[10 + i]] <-
    c(i_pct1["TRUE"], i_pct2["TRUE"], i_pct3["TRUE"])
}

p3 <- ggplot(pct_pcnc[pct_pcnc$Type != "Protein-coding", ],
  aes(Sample, Percent, fill = Type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 17)) +
  facet_wrap(~factor(Timepoint, levels = c("DS", "S24", "S72", "L6", "L26", "L57")),
    nrow = 1, scales = "free_x", strip.position = "bottom") +
  scale_fill_manual(values = c("black", "grey50"), name = element_blank()) +
  xlab(element_blank()) +
  ylab("% of transcription") +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(colour = "black", size = unit(6, "pt")),
    panel.grid = element_blank(),
    axis.text.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.title.y = element_text(colour = "black", size = unit(6, "pt")),
    axis.ticks.x = element_blank(),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    legend.position = "top",
    panel.spacing = unit(1, "mm"),
    axis.text.x = element_blank()
  )

ggsave("Figure2/SuppFig2a_Pct_noncoding_transcription.pdf", p3, width = 1.765, height = 1.773)

# Supplementary Figure 2b ------------------------------------------------------

rCor <- cor(t(rnaMz_m), t(tssMz_m))

rovCoef <- matrix(0, ncol = 6, nrow = 6, dimnames = list(paste0("R", 1:6), paste0("C", 1:6)))
rovFish <- matrix(0, ncol = 6, nrow = 6, dimnames = list(paste0("R", 1:6), paste0("C", 1:6)))
for (i in 1:nrow(rovCoef)) {
  for (j in 1:ncol(rovCoef)) {
    rovCoef[i, j] <- sum(rna_g[[i]] %in% tss_g[[j+1]]) /
      length(unique(c(rna_g[[i]], tss_g[[j+1]])))
    rovFish[i, j] <- fisher.test(alt = "g",
      matrix(c(
          sum(rna_g[[i]] %in% tss_g[[j+1]]),
          sum(!rna_g[[i]] %in% tss_g[[j+1]]),
          sum(!tss_g[[j+1]] %in% rna_g[[i]]),
          length(unique(c(unlist(tss_g, use.names=FALSE), unlist(rna_g, use.names=FALSE)))) -
            length(unique(c(rna_g[[i]], tss_g[[j+1]])))
        ), ncol = 2))$p.value
  }
}
print(rovFish < 1e-6)

pdf("Figure2/SuppFig2b_comparing_clusters_rna.pdf", width = 3.525, height = 1.6)
Heatmap(rCor,
  cluster_columns = FALSE, cluster_rows = FALSE, name = "PCC",
  col = circlize::colorRamp2(c(-1, 0, 1),
    c(colorBlindness::safeColors[c("vermillion", "skyBlue")], "white")[c(2, 3, 1)]),
  column_names_rot = 0, column_names_centered = TRUE,
  row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), at = c(-1, 0, 1),
    title_gp = gpar(fontsize = 5), legend_height = unit(8, "mm")),
  column_title = "csRNA-seq clusters", row_title = "RNA-seq clusters") +
Heatmap(rovCoef,
  cluster_rows = FALSE, cluster_columns = FALSE,
  col = viridisLite::inferno(100, begin = 0.00, end = 1),
  row_names_gp = gpar(fontsize = 6), column_names_gp = gpar(fontsize = 6),
  row_title_gp = gpar(fontsize = 6), column_title_gp = gpar(fontsize = 6),
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 5), at = c(0, 0.3),
    title_gp = gpar(fontsize = 5), legend_height = unit(8, "mm")),
  name = "Overlap Coef.", column_names_rot = 0, column_names_centered = TRUE
)
dev.off()

# Supplementary Figure 2c ------------------------------------------------------

gene2type <- structure(as.character(a11$locus_type), names = a11$ID)

rnatypes <- data.frame(row.names = NULL,
  Transcript = unname(rna2cls),
  Cluster = names(rna2cls),
  Type = gene2type[gsub("[.]\\d+", "", unname(rna2cls))]
)
rnatypes$Type2 <- rnatypes$Type
rnatypes$Type2[is.na(rnatypes$Type)] <- "Unknown"
rnatypes$Type2[rnatypes$Type2 == "protein_coding"] <- "Protein-coding"
rnatypes$Type2[!rnatypes$Type2 %in% c("Unknown", "Protein-coding")] <- "Non-coding"

rnatypes2 <- dplyr::summarize(
  dplyr::group_by(rnatypes[, c(2, 4)], Cluster, Type2),
  sum = dplyr::n())
rnatypes2$prop <- rnatypes2$sum / structure(sapply(rna_t, length), names = names(rna_t))[rnatypes2$Cluster]
rnatypes2$Type3 <- factor(rnatypes2$Type2, levels = c("Protein-coding", "Non-coding", "Unknown"))

p4 <- ggplot(rnatypes2, aes(Cluster, prop, group = Type3, fill = Type3)) +
  geom_bar(stat = "identity", width = .8) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("RNA-seq cluster") +
  ylab("Proportion") +
  scale_fill_manual(values = rev(awtools::ppalette[-(2:5)])[c(1, 4, 3, 2)],
    name = element_blank()) +
  theme_bw() +
  theme(
    aspect.ratio = 0.6667,
    panel.border = element_rect(colour = "black", size = 1),
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black", size = unit(6, "pt")),
    legend.text = element_text(colour = "black", size = unit(5, "pt")),
    legend.key.size = unit(2, "mm"),
    axis.text = element_text(colour = "black", size = unit(6, "pt"))
  )

ggsave("Figure2/SuppFig2c_RNA_cls_types.pdf", p4, width = 3.5, height = 1.1)

# Supplementary Figure 2d ------------------------------------------------------

rna_go_wc <- anno_word_cloud_from_GO(
  rep(names(rna_t), times = sapply(rna_t, length)),
  lapply(rna_go, function(x) x$term_id[x$p_value < 0.01]),
  fontsize_range = c(2, 8),
  bg_gp = gpar(col = "#000000", fill = "#ffffff"),
  word_cloud_grob_param = list(max_width = unit(30, "mm"), word_space = unit(2, "pt"),
    line_space = unit(2, "pt")),
  stat = "pvalue", max_words = 15, min_stat = 1,
  exclude_words = c("process", "metabolic", "response", "cellular", "cell",
    "compound", "biological", "macromolecule", "levels", "decreased", "organic"))

pdf("Figure2/SuppFig2d_rna_heatmap.pdf", height = 2.4, width = 3.9)
Heatmap(rnaMz[rna2cls, ], row_split = names(rna2cls),
  col = circlize::colorRamp2(c(-1.9, 0, 1.9), c("blue", "white", "red")),
  cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE,
  show_column_names = FALSE, row_title_rot = 0, row_title_gp = gpar(fontsize = 6),
  cluster_column_slices = FALSE, cluster_row_slices = FALSE,
  show_row_dend = FALSE, show_column_dend = FALSE,
  column_order = colnames(rnaMz), row_order = unname(rna2cls),
  column_gap = unit(0.2, "mm"), row_gap = unit(0.6, "mm"),
  raster_quality = 3,
  right_annotation = rowAnnotation(wc = rna_go_wc))
dev.off()

# Supplementary Figure 2e ------------------------------------------------------

tssCls <- readr::read_tsv(TSS_CLUSTERS)
TFs <- readr::read_tsv(TF_LIST)

tssClsM <- as.matrix(tssCls[, 11:22])
rownames(tssClsM) <- tssCls$TSS
tssClsM <- tssClsM[order(matrixStats::rowMedians(tssClsM), decreasing=T), ]

tssTf <- tssCls[!is.na(tssCls$NearestGene) & tssCls$NearestGene %in% TFs$Gene_ID, ]
tssMtf <- tssClsM[rownames(tssClsM) %in% tssTf$TSS, ]
tssMtf <- tssMtf[!duplicated(unname(structure(tssTf$NearestGene, names=tssTf$TSS)[rownames(tssMtf)])), ]
tssTf <- tssTf[tssTf$TSS %in% rownames(tssMtf), ]

gene2sym <- structure(a11$symbol, names=a11$ID)
gene2sym[is.na(unname(gene2sym))] <- names(gene2sym)[is.na(unname(gene2sym))]
tss2gene <- structure(tssTf$NearestGene, names = tssTf$TSS)

tssTf$MoreThan10 <- apply(tssMtf[tssTf$TSS,], 1, function(x) sum(x>=10) >= 2)
tssTf$MoreThan50 <- apply(tssMtf[tssTf$TSS,], 1, function(x) sum(x>=50) >= 2)
tssTf$MoreThan100 <- apply(tssMtf[tssTf$TSS,], 1, function(x) sum(x>=100) >= 2)

tf_count_plot <- function(x) {
  ggplot(data.frame(cluster = paste0("C", 0:6), size = as.vector(table(x)))[-1,],
    aes(cluster, size)) +
    geom_bar(stat = "identity", fill = "black", width = 0.75) +
    scale_y_continuous(expand = c(0, 0), name = "# of TFs") +
    xlab(element_blank()) +
    theme_bw() +
    theme(
      axis.text.x = element_text(colour = "black", size = unit(5, "pt")),
      legend.position = "none",
      axis.text.y = element_text(colour = "black", size = unit(5, "pt")),
      axis.title.y = element_text(colour = "black", size = unit(5, "pt")),
      axis.title.x = element_text(colour = "black", size = unit(5, "pt")),
      axis.line.x.bottom = element_line(colour = "black", size = .5),
      axis.line.y.left = element_line(colour = "black", size = .5),
      legend.key.size = unit(5, "pt"),
      legend.text = element_text(colour = "black", size = unit(5, "pt")),
      panel.grid = element_blank(),
      panel.border = element_blank()
    )
}

tf1 <- tf_count_plot(tssTf$Cluster)
tf2 <- tf_count_plot(tssTf$Cluster[tssTf$MoreThan10])
tf3 <- tf_count_plot(tssTf$Cluster[tssTf$MoreThan50])
tf4 <- tf_count_plot(tssTf$Cluster[tssTf$MoreThan100])

cowplot::save_plot("Figure2/SuppFig2e_TF_counts.pdf", base_width = 5, base_height = 1,
  cowplot::plot_grid(nrow = 1, tf1, tf2, tf3, tf4))

# Supplementary Figure 2f ------------------------------------------------------

rownames(tssMtf) <- tss2gene[rownames(tssMtf)]
tssMtf <- tssMtf[gtools::mixedorder(rownames(tssMtf)), ]
tssMtf <- tssMtf[apply(tssMtf, 1, function(x) sum(x>=100) >= 2), ]

tssTf <- as.data.frame(tssTf)
rownames(tssTf) <- tssTf$NearestGene
tssTf <- tssTf[rownames(tssMtf), ]
tssTf <- tssTf[tssTf$Cluster != 0, ]
tssMtf <- tssMtf[rownames(tssTf), ]

print(nrow(tssMtf))
print(table(tssTf$Cluster))

pdf("Figure2/SuppFig2f_TF_heatmap.pdf", height = 3.5, width = 3.2)
Heatmap(t(scale(t(tssMtf))),
  col = circlize::colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  name = "Z-score",
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  row_names_gp = gpar(fontsize = 2),
  row_split = paste0("C", tssTf$Cluster),
  show_column_names = FALSE, row_title_rot = 0, row_title_gp = gpar(fontsize = 6),
  column_gap = unit(0.2, "mm"), row_gap = unit(0.6, "mm"),
  row_labels = gene2sym[rownames(tssMtf)]
)
dev.off()
