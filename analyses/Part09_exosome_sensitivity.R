# Benjamin J.M. Tremblay, June 2023
#
# Exosome sensitive TSSs and transcripts


Anno <- readr::read_tsv("TSS.annotations.tsv")
Tss1 <- import("TSSs_size1.bed")

cs_L57_H <- read.table("csRNAseq.DE.L57_hen2.txt", header = TRUE, sep = "\t")
cs_L57_R <- read.table("csRNAseq.DE.L57_rrp4.txt", header = TRUE, sep = "\t")

r_L57_H <- read.table("RNAseq.DE.L57_hen2.txt", header = TRUE, sep = "\t")
r_L57_R <- read.table("RNAseq.DE.L57_rrp4.txt", header = TRUE, sep = "\t")

# Classify TSS/transcripts as coding/non-coding

tss_unknown <- Tss1$name[Tss1$name %in% Anno$TSS[is.na(Anno$NearestGeneType)]]
tss_putative <- Tss1$name[Tss1$name %in% Anno$TSS[!is.na(Anno$NearestGeneType) &
  Anno$NearestGeneType == "putative_rna"]]
tss_other <- Tss1$name[Tss1$name %in% Anno$TSS[!is.na(Anno$NearestGeneType) &
  Anno$NearestGeneType %in% c("miRNA", "pre_trna", "small_nuclear_rna", "small_nucleolar_rna")]]
tss_lncrna <- Tss1$name[Tss1$name %in% Anno$TSS[!is.na(Anno$NearestGeneType) &
  Anno$NearestGeneType %in% c("antisense_long_noncoding_rna", "antisense_rna", "long_noncoding_rna", "novel_transcribed_region", "other_rna", "pseudogene", "transposable_element_gene")]]
tss_mrna <- Tss1$name[Tss1$name %in% Anno$TSS[!is.na(Anno$NearestGeneType) &
  Anno$NearestGeneType == "protein_coding"]]

tx_putative <- unique(Anno$NearestFeature[Anno$TSS %in% tss_putative])
tx_lncrna <- unique(Anno$NearestFeature[Anno$TSS %in% tss_lncrna])
tx_mrna <- unique(Anno$NearestFeature[Anno$TSS %in% tss_mrna])

# Only look at non-coding RNAs

cs_L57_H <- cs_L57_H[rownames(cs_L57_H %in% c(tss_lncrna, tss_unknown, tss_putative)), ]
cs_L57_R <- cs_L57_R[rownames(cs_L57_R %in% c(tss_lncrna, tss_unknown, tss_putative)), ]

r_L57_H <- r_L57_H[rownames(r_L57_H %in% c(tx_putative, tx_lncrna)), ]
r_L57_R <- r_L57_R[rownames(r_L57_R %in% c(tx_putative, tx_lncrna)), ]

# Identify putatively exosome-(in)sensitive TSSs/transcripts

cs_Ins_H <- rownames(cs_L57_H)[cs_L57_H$logFC < 1 & cs_L57_H$FDR >= 0.05]
cs_Ins_R <- rownames(cs_L57_R)[cs_L57_R$logFC < 1 & cs_L57_R$FDR >= 0.05]
rn_Ins_H <- rownames(cs_L57_H)[r_L57_H$logFC < 1 & r_L57_H$FDR >= 0.05]
rn_Ins_R <- rownames(cs_L57_R)[r_L57_R$logFC < 1 & r_L57_R$FDR >= 0.05]

cs_Sen_H <- rownames(cs_L57_H)[cs_L57_H$logFC >= 1 & cs_L57_H$FDR < 0.05]
cs_Sen_R <- rownames(cs_L57_R)[cs_L57_R$logFC >= 1 & cs_L57_R$FDR < 0.05]
rn_Sen_H <- rownames(cs_L57_H)[r_L57_H$logFC >= 1 & r_L57_H$FDR < 0.05]
rn_Sen_R <- rownames(cs_L57_R)[r_L57_R$logFC >= 1 & r_L57_R$FDR < 0.05]

readr::write_lines(rownames(cs_L57_H), "all_nc_tss.txt")
readr::write_lines(rownames(r_L57_H), "all_nc_rna.txt")

readr::write_lines(cs_Ins_H[cs_Ins_H %in% cs_Ins_R],
  "Exosome_insensitive_tss.txt")
readr::write_lines(rn_Ins_H[rn_Ins_H %in% rn_Ins_R],
  "Exosome_insensitive_rn.txt")

readr::write_lines(cs_Sen_H[cs_Sen_H %in% cs_Sen_R],
  "Exosome_sensitive_tss.txt")
readr::write_lines(rn_Sen_H[rn_Sen_H %in% rn_Sen_R],
  "Exosome_sensitive_rn.txt")

