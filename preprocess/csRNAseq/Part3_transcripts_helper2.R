# Benjamin J.M. Tremblay, June 2023
#
# Final adjustments to de novo transcripts

library(rtracklayer)

tss <- import("csrnaseq_TSSs_final/tss.shrunk.bed")
denovo <- import("rnaseq_denovodenovo/stringtie_homer.combined.gtf")
anno <- import("Araport11_GFF3_genes_transposons.gff3")

# Only keep de novo transcripts not overlapping known transcripts, and overlapping a TSS

tssBED <- tss
denovoGTF <- denovo
denovoGTF <- denovoGTF[is.na(findOverlaps(denovoGTF, anno, select = "first"))]
denovoTX <- resize(denovoGTF[as.character(denovoGTF$type) == "transcript"], 50)
denovoTX <- denovoTX[!is.na(findOverlaps(denovoTX, tss, select = "first"))]

denovoTX$TSS <- tss$name[findOverlaps(denovoTX, tss, select = "first")]
denovoTX <- denovoTX[order(width(denovoTX), decreasing = TRUE)]
denovoTX <- denovoTX[!duplicated(denovoTX$TSS)]

denovo <- denovo[denovo$transcript_id %in% denovoTX$transcript_id]

export(denovo, "rnaseq_denovo/novel.transcripts.filtered.gtf")

