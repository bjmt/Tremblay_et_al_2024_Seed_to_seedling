# Benjamin Jean-Marie Tremblay, June 2023
# 
# Generate a final set of de novo transcripts

library(rtracklayer)

stringtie <- import("rnaseq_denovo/stringtie.merged.combined.gtf")
homer <- import("rnaseq_denovo/homer.merged.combined.gtf")
anno <- import("Araport11_GFF3_genes_transposons.gff3")

stringtie_tx <- stringtie[as.character(stringtie$type) == "transcript"]
homer_tx <- homer[as.character(homer$type) == "transcript"]
anno_tx <- anno[as.character(anno$type) == "transcript"]

stringtie_exons <- tapply(stringtie$exon_number, stringtie$transcript_id,
  function(x) sum(!is.na(x)))
str_df <- data.frame(
  row.names = NULL,
  Gene = stringtie_tx$gene_id,
  Transcript = stringtie_tx$transcript_id,
  Width = width(stringtie_tx),
  Exons = stringtie_exons[stringtie_tx$transcript_id]
)
str_df_split <- lapply(unique(str_df$Gene), function(x) str_df[str_df$Gene == x, ])
for (i in seq_along(str_df_split)) {
  if (nrow(str_df_split[[i]]) == 1) next
  str_df_split[[i]] <- str_df_split[[i]][
    str_df_split[[i]]$Exons == max(str_df_split[[i]]$Exons),
  ]
  if (nrow(str_df_split[[i]]) > 1) {
    str_df_split[[i]] <- str_df_split[[i]][
      which.max(str_df_split[[i]]$Width),
    ]
  }
}
str_df_split <- do.call(rbind, str_df_split)
stringtie_tx <- stringtie_tx[stringtie_tx$transcript_id %in% str_df_split$Transcript]

# Only keep transcripts at least 200 bp
homer_tx <- homer_tx[width(homer_tx) >= 200]
stringtie_tx <- stringtie_tx[width(stringtie_tx) >= 200]

homer_tx <- homer_tx[is.na(findOverlaps(homer_tx, anno_tx, select = "first"))]
stringtie_tx <- stringtie_tx[is.na(findOverlaps(stringtie_tx, anno_tx, select = "first"))]
# For overlapping HOMER and StringTie transcripts, prioritize those from StringTie
homer_tx <- homer_tx[is.na(findOverlaps(homer_tx, stringtie_tx, select = "first"))]

homer <- homer[homer$transcript_id %in% homer_tx$transcript_id]
stringtie <- stringtie[stringtie$transcript_id %in% stringtie_tx$transcript_id]

stringtieLocs <- structure(gsub("RNA", "LOC", stringtie_tx$transcript_id),
  names = stringtie_tx$transcript_id)
homerLocs <- structure(gsub("RNA", "LOC", homer_tx$transcript_id),
  names = homer_tx$transcript_id)

homer$gene_id <- homerLocs[homer$gene_id]
stringtie$gene_id <- stringtieLocs[stringtie$gene_id]

export(c(homer, stringtie), "rnaseq_denovodenovo/stringtie_homer.combined.gtf")

