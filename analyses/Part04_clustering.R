# Benjamin J.M. Tremblay, June 2023
#
# Clustering of the csRNA-seq, RNA-seq, ATAC-seq

library(WGCNA)
library(matrixStats)

# g:Profiler version e106

gprofiler2::set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e106_eg53_p16")

# Functions --------------------------------------------------------------------

calc_rep_repro <- function(x) {
  y <- abs(c(diff(x[1:2]), diff(x[3:4]), diff(x[5:6]), diff(x[7:8])))
  y / sd(x)
}

filter_clus <- function(x, sdMaxRepDiff = .25, minMem = 1000) {
  # This function aims to find large clusters consistent between replicates
  f_count <- table(x$colors)
  f_count <- names(f_count[f_count >= minMem])
  y <- apply(x$newMEs, 2, calc_rep_repro)
  f <- apply(y, 2, function(z) all(z <= sdMaxRepDiff))
  f <- as.numeric(gsub("ME", "", colnames(x$newMEs)[f]))
  f[f %in% f_count]
}

run_go <- function(genes, file) {
  for (attempt in 1:5) {
    go <- try(gprofiler2::gost(query = genes, org = "athaliana"))
    if (!inherits(go, "try-error")) break
    Sys.sleep(30)
  }
  if (is.data.frame(go$result) && nrow(go$result)) {
    go <- go$result
    go <- go[order(go$source, go$p_value), ]
    readr::write_tsv(go, file)
  }
}

TSS_table <- readr::read_tsv("tss.annotated.gro.tsv")

# csRNA-seq --------------------------------------------------------------------

# Read in CPM table (for Col-0 samples only)

CPMs <- readr::read_tsv("csRNAseq.quantification.cpm.txt")
CPMsM <- as.matrix(CPMs[, -1])
rownames(CPMsM) <- CPMs$TSS
CPMsM <- CPMsM[, 1:12]

# Find right power value

powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2));
softhresh <- pickSoftThreshold(t(CPMsM))
plot(softhresh$fitIndices$Power, softhresh$fitIndices$mean.k)

# Cluster

net <- blockwiseModules(t(CPMsM), power = 3,
  TOMType = "signed", minModuleSize = 30,
  numericLabels = TRUE, networkType = "signed",
  saveTOMs = TRUE, saveTOMFileBase = "csRNA_TOM",
  verbose = 3, maxBlockSize = 60000)

# Final clusters

Modules <- c("3" = "1", "7" = "1", "5" = "2", "4" = "3", "6" = "4", "2" = "5", "1" = "6")

Clusters <- net$colors[net$colors %in% 1:7]
Clusters <- data.frame(row.names = NULL,
  TSS = names(Clusters),
  Cluster = Modules[as.character(unname(Clusters))]
)
Clusters <- Clusters[order(Clusters$Cluster), ]

TSS_anno <- readr::read_tsv("TSS.annotations.tsv")
TSS_anno$Cluster <- paste0("C", structure(Clusters$Cluster, names = Clusters$TSS)[TSS_anno$TSS])
TSS_anno$Cluster[TSS_anno$Cluster == "CNA"] <- "C0"

readr::write_tsv(TSS_anno, "TSS.annotations.clusters.tsv")

# Cluster tables, genes and GO enrichment

dir.create("csRNAseq_clustering/clusters-final-genes", recursive = TRUE, showWarnings = FALSE)
dir.create("csRNAseq_clustering/clusters-final-genes-go", recursive = TRUE, showWarnings = FALSE)

Clusters_anno <- rbind(
  data.frame(Cluster = 0L, TSS = TSS_table$TSS[!TSS_table$TSS %in% Clusters$TSS]),
  data.frame(Cluster = as.integer(Clusters$Cluster), TSS = Clusters$TSS)
)
Clusters_anno <- dplyr::left_join(Clusters_anno, TSS_table, by = "TSS")

readr::write_tsv(Clusters_anno, "csRNAseq_clustering/Clusters_anno.tsv")

ClusterGenes <- list()
for (i in unique(Clusters$Cluster)) {
  Genes <- TSS_table$NearestGene[TSS_table$TSS %in% Clusters$TSS[Clusters$Cluster == i]]
  Genes <- unique(gsub("[.]\\d+", "", Genes[!is.na(Genes)]))
  ClusterGenes[[paste0("C", i)]] <- Genes
}
Genes <- TSS_table$NearestGene[!TSS_table$TSS %in% Clusters$TSS]
Genes <- ClusterGenes$C6[!is.na(Genes)]
ClusterGenes$C0 <- unique(gsub("[.]\\d+", "", Genes))

for (i in names(ClusterGenes)) {
  readr::write_lines(ClusterGenes[[i]],
    paste0("csRNAseq_clustering/clusters-final-genes/Genes_", i, ".txt"))
  run_go(ClusterGenes[[i]][grepl("^AT", ClusterGenes[[i]])],
    paste0("csRNAseq_clustering/clusters-final-genes-go/Cluster_", i, ".txt"))
}

# RNA-seq ----------------------------------------------------------------------

TPMs <- readr::read_tsv("RNAseq.quantifications.tpm.tsv")
TPMsM <- as.matrix(TPMs[, -1])
rownames(TPMsM) <- TPMs$Transcript
TPMsM <- TPMsM[, 1:12]

# Only keep higher expression isoforms

TPMsM <- TPMsM[apply(TPMsM, 1, function(x) sum(x >= 2) >= 2), ]

# Find right power value

powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2));
softhresh <- pickSoftThreshold(t(TPMsM))
plot(softhresh$fitIndices$Power, softhresh$fitIndices$mean.k)

# Cluster

net <- blockwiseModules(t(TPMsM), power = 6,
  TOMType = "signed", minModuleSize = 30,
  numericLabels = TRUE, networkType = "signed",
  saveTOMs = TRUE, saveTOMFileBase = "RNA_TOM",
  verbose = 3, maxBlockSize = 60000)
net <- WGCNA::mergeCloseModules(t(TPMsM), net$colors, net$MEs, cutHeight = 0.25)

# Final clusters: keep large clusters

keep <- filter_clus(net, 0.6)

Modules <- c("7" = "1", "3" = "2", "8" = "3", "1" = "4", "2" = "5", "4" = "6")

Clusters <- net$colors[net$colors %in% keep]
Clusters <- data.frame(row.names = NULL,
  Transcript = names(Clusters),
  Cluster = paste0("R", Modules[as.character(unname(Clusters))])
)
Clusters <- Clusters[sample.int(nrow(Clusters)), ]
Clusters <- Clusters[order(Clusters$Cluster), ]

readr::write_tsv(Clusters, "Trancripts.clusters.tsv")

# Cluster transcripts and GO enrichment

dir.create("RNAseq_clustering/clusters-final-transcripts", recursive = TRUE, showWarnings = FALSE)
dir.create("RNAseq_clustering/clusters-final-genes-go", recursive = TRUE, showWarnings = FALSE)

for (i in sort(unique(Clusters$Cluster))) {
  Transcripts <- Clusters$Transcript[Clusters$Cluster == i]
  readr::write_lines(Transcripts,
    paste0("RNAseq_clustering/clusters-final-transcripts/Transcripts_", i, ".txt"))
  run_go(unique(gsub("[.]\\d+", "", Transcripts)),
    paste0("RNAseq_clustering/clusters-final-genes-go/Cluster_", i, ".txt"))
}

# ATAC-seq ---------------------------------------------------------------------

RPMs <- readr::read_tsv("ATACseq.quantification.rpm.txt")
RPMsM <- as.matrix(RPMs[, -1])
rownames(RPMsM) <- RPMs$PeakID
RPMsM <- RPMsM[, 1:8]

# Only keep strong more variable peaks

RPMsM <- RPMsM[log2(rowVars(RPMsM)) >= 4, ]

# Find right power value

powers <- c(seq(4, 10, by = 1), seq(12, 20, by = 2));
softhresh <- pickSoftThreshold(t(RPMsM))
plot(softhresh$fitIndices$Power, softhresh$fitIndices$mean.k)

# Cluster

net <- blockwiseModules(t(RPMsM), power = 5,
  TOMType = "signed", minModuleSize = 30,
  numericLabels = TRUE, networkType = "signed",
  saveTOMs = TRUE, saveTOMFileBase = "ATAC_TOM",
  verbose = 3, maxBlockSize = 60000)
net <- WGCNA::mergeCloseModules(t(RPMsM), net$colors, net$MEs, cutHeight = 0.35)

# Final clusters: keep large clusters

keep <- filter_clus(net, 0.5)

Modules <- c("1" = "1", "3" = "2", "8" = "3", "4" = "4", "2" = "5")

Clusters <- net$colors[net$colors %in% keep]
Clusters <- data.frame(row.names = NULL,
  PeakID = names(Clusters),
  Cluster = Modules[as.character(unname(Clusters))]
)

ACR_anno <- readr::read_tsv("ACR.annotations.tsv")
ACR_anno$Cluster <- paste0("A", structure(Clusters$Cluster, names = Clusters$PeakID)[ACR_anno$PeakID])
ACR_anno$Cluster[ACR_anno$Cluster == "ANA"] <- "A0"

readr::write_tsv(ACR_anno, "ACR.annotations.clusters.tsv")

# Cluster tables, genes and GO enrichment

dir.create("ATACseq_clustering3/clusters-final-genes", recursive = TRUE, showWarnings = FALSE)
dir.create("ATACseq_clustering3/clusters-final-genes-go", recursive = TRUE, showWarnings = FALSE)

ATAC_cpm <- as.matrix(read.table("ATAC_peaks_cpm_fixed.txt", check.names = FALSE))

Clusters_anno <- readr::read_tsv("ACR.annotations.tsv")
Clusters_anno <- cbind(Cluster = structure(Clusters$Cluster,
  names = Clusters$PeakID)[Clusters_anno$PeakID], Clusters_anno)
Clusters_anno$Cluster[is.na(Clusters_anno$Cluster)] <- 0
Clusters_anno$Cluster <- as.integer(Clusters_anno$Cluster)
Clusters_anno$Gene <- structure(TSS_table$NearestGene,
  names = TSS_table$TSS)[Clusters_anno$NearestTSS]
Clusters_anno <- cbind(Clusters_anno, ATAC_cpm[Clusters_anno$PeakID, ])

readr::write_tsv(Clusters_anno, "ATACseq_clustering3/Clusters_anno.tsv")

for (i in 0:5) {
  Genes <- Clusters_anno$Gene[Clusters_anno$Cluster == i]
  Genes <- unique(Genes[!is.na(Genes)])
  readr::write_lines(Genes, paste0("ATACseq_clustering3/clusters-final-genes/Cluster_A", i, ".txt"))
  run_go(Genes, paste0("ATACseq_clustering3/clusters-final-genes-go/Cluster_A", i, ".txt"))
}
