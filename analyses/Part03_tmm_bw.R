# Benjamin J.M. Tremblay, June 2023
#
# Generate TMM-normalized bigWigs of the csRNA-seq and RNA-seq

library(rtracklayer)

csRNA_tmm <- read.table("csRNAseq.norm_factors.txt")
RNA_tmm <- read.table("RNAseq.norm_factors.txt")

csRNA_tmm$scale <- with(csRNA_tmm, ((norm.factors * lib.size) / 1e6)^-1)
csRNA_tmm$final.lib <- with(csRNA_tmm, lib.size * scale)
csRNA_tmm$final.norm <- with(csRNA_tmm, final.lib / 1e6)
csRNA_tmm <- with(csRNA_tmm, structure(tapply(final.norm, group, mean)))

RNA_tmm$scale <- with(RNA_tmm, ((norm.factors * lib.size) / 1e6)^-1)
RNA_tmm$final.lib <- with(RNA_tmm, lib.size * scale)
RNA_tmm$final.norm <- with(RNA_tmm, final.lib / 1e6)
RNA_tmm <- with(RNA_tmm, structure(tapply(final.norm, group, mean)))

chrom.sizes <- read.table("chrom.sizes")
chrom.sizes <- structure(chrom.sizes$V2, names = chrom.sizes$V1)

for (i in c("DS", "S24", "S72", "L6", "L26", "L57", "hen2", "rrp4")) {
  bg_plus <- import(paste0("csrnseq_bedGraphs/", i, ".merged.plus.bedGraph.gz"))
  bg_minus <- import(paste0("csrnaseq_bedGraphs/", i, ".merged.plus.bedGraph.gz"))
  bg_plus$score <- bg_plus$score * csRNA_tmm[i]
  bg_minus$score <- bg_minus$score * csRNA_tmm[i]
  seqlengths(bg_plus) <- chrom.sizes
  seqlengths(bg_minus) <- chrom.sizes
  export.bw(bg_plus, paste0(i, ".merged.csRNA.plus.bw"))
  export.bw(bg_minus, paste0(i, ".merged.csRNA.minus.bw"))
}

for (i in c("DS", "S24", "S72", "L6", "L26", "L57", "hen2", "rrp4")) {
  bg_plus <- import(paste0("rnseq_bedGraphs/", i, ".merged.plus.bedGraph.gz"))
  bg_minus <- import(paste0("rnseq_bedGraphs/", i, ".merged.plus.bedGraph.gz"))
  bg_plus$score <- bg_plus$score * csRNA_tmm[i]
  bg_minus$score <- bg_minus$score * csRNA_tmm[i]
  seqlengths(bg_plus) <- chrom.sizes
  seqlengths(bg_minus) <- chrom.sizes
  export.bw(bg_plus, paste0(i, ".merged.RNA.plus.bw"))
  export.bw(bg_minus, paste0(i, ".merged.RNA.minus.bw"))
}

