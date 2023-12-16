# Benjamin J.M. Tremblay, June 2023
#
# Clustering of the csRNA-seq, RNA-seq, ATAC-seq

library(WGCNA)
library(matrixStats)

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
net <- WGCNA::mergeCloseModules(t(CPMsM), net$colors, net$MEs, cutHeight = 0.25)

# Reorder/rename clusters as needed, keep large clusters

filter_clus <- function(x, sdMaxRepDiff = .25, minMem = 1000) {
  # This function aims to find large clusters consistent between replicates
  f_count <- table(x$colors)
  f_count <- names(f_count[f_count >= minMem])
  y <- apply(x$newMEs, 2, calc_rep_repro)
  f <- apply(y, 2, function(z) all(z <= sdMaxRepDiff))
  f <- as.numeric(gsub("ME", "", colnames(x$newMEs)[f]))
  f[f %in% f_count]
}

keep <- filter_clus(net)
clusters <- net$colors[net$colors %in% keep]

# Extract final cluster membership information

clusters <- data.frame(row.names = NULL,
  TSS = names(clusters),
  Cluster = paste0("C", names(clusters))
)
clusters <- clusters[sample.int(nrow(clusters)), ]
clusters <- clusters[order(clusters$Cluster), ]

TSS_anno <- readr::read_tsv("TSS.annotations.tsv")
TSS_anno$Cluster <- structure(clusters$Cluster, names = clusters$TSS)[TSS_anno$TSS]
TSS_anno$Cluster[is.na(TSS_anno$Cluster)] <- "C0"

# Export

readr::write_tsv(TSS_anno, "TSS.annotations.clusters.tsv")

# GO enrichment

for (i in unique(clusters$Cluster)) {
  Genes <- TSS_anno$NearestGene[TSS_anno$Cluster == i]
  Genes <- unique(Genes[!is.na(Genes)])
  go <- gprofiler2::gost(query = Genes, org = "athaliana")
  go <- go$result
  go <- go[order(go$source, go$p_value), ]
  readr::write_tsv(go, paste("Cluster.GO.", i, ".tsv"))
}

# RNA-seq ----------------------------------------------------------------------

TPMs <- readr::read_tsv("RNAseq.quantifications.tpm.tsv")
TPMsM <- as.matrix(TPMs[, -1])
rownames(TPMsM) <- TPMs$Transcript
TPMsM <- TPMsM[, 1:12]

# Only keep higher expression isoforms

TPMsM <- TPMsM[apply(TPMsM, 1, function(sum(x >= 2) >= 2)), ]

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

# Reorder/rename clusters as needed, keep large clusters

keep <- filter_clus(net, 0.6)
clusters <- net$colors[net$colors %in% keep]

# Extract final cluster membership information

clusters <- data.frame(row.names = NULL,
  Transcript = names(clusters),
  Cluster = paste0("R", names(clusters))
)
clusters <- clusters[sample.int(nrow(clusters)), ]
clusters <- clusters[order(clusters$Cluster), ]

# Export

readr::write_tsv(clusters, "Trancripts.clusters.tsv")

# GO enrichment

for (i in unique(clusters$Cluster)) {
  Genes <- clusters$Transcript[clusters$Cluster == i]
  Genes <- Genes[!is.na(Genes)]
  Genes <- unique(gsub("[.]\\d+", "", Genes))
  go <- gprofiler2::gost(query = Genes, org = "athaliana")
  go <- go$result
  go <- go[order(go$source, go$p_value), ]
  readr::write_tsv(go, paste("Cluster.GO.", i, ".tsv"))
}

# ATAC-seq ---------------------------------------------------------------------

RPMs <- readr::read_tsv("ATACseq.quantifications.rpm.txt")
RPMsM <- as.matrix(RPMs[, -1])
rownames(RPMsM) <- RPMs$PeakID
RPMsM <- RPMsM[, 1:12]

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

# Reorder/rename clusters as needed, keep large clusters

keep <- filter_clus(net, 0.5)
clusters <- net$colors[net$colors %in% keep]

# Extract final cluster membership information

clusters <- data.frame(row.names = NULL,
  PeakID = names(clusters),
  Cluster = paste0("A", names(clusters))
)
clusters <- clusters[sample.int(nrow(clusters)), ]
clusters <- clusters[order(clusters$Cluster), ]

ACR_anno <- readr::read_tsv("ACR.annotations.tsv")
ACR_anno$Cluster <- structure(clusters$Cluster, names = Clusters$PeakID)[ACR_anno$PeakID]
ACR_anno$Cluster[is.na(ACR_anno$Cluster)] <- "A0"

# Export

readr::write_tsv(ACR_anno, "ACR.annotations.clusters.tsv")

# GO enrichment

for (i in unique(clusters$Cluster)) {
  Genes <- ACR_anno$NearesTSS[ACR_anno$Cluster == i]
  Genes <- TSS_anno$NearestGene[TSS_anno$TSS %in% Genes]
  Genes <- unique(Genes[!is.na(Genes)])
  go <- gprofiler2::gost(query = Genes, org = "athaliana")
  go <- go$result
  go <- go[order(go$source, go$p_value), ]
  readr::write_tsv(go, paste("Cluster.GO.", i, ".tsv"))
}

