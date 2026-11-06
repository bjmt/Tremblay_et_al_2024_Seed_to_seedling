track_anno <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "anno")
}

build_anno <- function(data, ranges, params) {
  if (!is(data, "TxDb")) {
    stop("[anno] Expected a TxDb object", call. = FALSE)
  }

  # Note: should this be 'any', thus allowing for partial drawings?
  data_g <- subsetByOverlaps(genes(data), ranges)#, type = "within")
  data_e <- exonsBy(data, "gene")[data_g$gene_id]
  
  g_df <- as.data.frame(data_g)
  g_df$strand <- as.character(g_df$strand)
  g_df$start[g_df$start < start(ranges)] <- start(ranges)
  g_df$end[g_df$end > end(ranges)] <- end(ranges)
  e_df <- as.data.frame(data_e)
  e_df$strand <- as.character(e_df$strand)
  e_df <- e_df[e_df$end >= start(ranges) & e_df$start <= end(ranges), ]
  e_df$start[e_df$start < start(ranges)] <- start(ranges)
  e_df$end[e_df$end > end(ranges)] <- end(ranges)
  p <- ggplot()
  # tx_ids <- unique(data_df$tx_id)
  g_p <- g_df[g_df$strand == "+", ]
  g_m <- g_df[g_df$strand == "-", ]
  g_p$ypos <- rep(1, sum(g_df$strand == "+"))
  g_m$ypos <- rep(1, sum(g_df$strand == "-"))
  buff <- (end(ranges) - start(ranges)) * params$anno.buffer
  if (nrow(g_p) > 1) {
    for (i in seq(2, nrow(g_p))) {
      for (j in seq_len(i - 1)) {
        if ((g_p$start[i] - buff) < g_p$end[j]) g_p$ypos[i] <- g_p$ypos[j] + 1
      }
    }
  }
  if (nrow(g_m) > 1) {
    for (i in seq(2, nrow(g_m))) {
      for (j in seq_len(i - 1)) {
        if ((g_m$start[i] - buff) < g_m$end[j]) g_m$ypos[i] <- g_m$ypos[j] + 1
      }
    }
  }
  if (nrow(g_p)) g_p$ypos <- g_p$ypos - 0.5
  if (nrow(g_m)) g_m$ypos <- -(g_m$ypos - 0.5)
  g_pm <- rbind(g_p, g_m)
  g_pm$ypos <- norm_fun(g_pm$ypos, -0.75, 0.75)
  e_df$ypos <- structure(g_pm$ypos, names = g_pm$gene_id)[e_df$group_name]
  p <- p + geom_segment(aes(x = start, xend = end, y = ypos, yend = ypos),
    data = g_pm[g_pm$strand == "+", ], lineend = "butt", colour = params$colour.pos,
    size = params$anno.gene.size)
  p <- p + geom_segment(aes(x = start, xend = end, y = ypos, yend = ypos),
    data = g_pm[g_pm$strand == "-", ], lineend = "butt", colour = params$colour.neg,
    size = params$anno.gene.size)
  p <- p + geom_segment(aes(x = start, xend = end, y = ypos, yend = ypos),
    data = e_df[e_df$strand == "+", ], lineend = "butt", colour = params$colour.pos,
    size = params$anno.exon.size)
  p <- p + geom_segment(aes(x = start, xend = end, y = ypos, yend = ypos),
    data = e_df[e_df$strand == "-", ], lineend = "butt", colour = params$colour.neg,
    size = params$anno.exon.size)
  if (!params$anno.text.hide) {
    # Plot text
  }
  p <- p + 
    scale_y_continuous(expand = c(0, 0), breaks = c(-1, 1),# oob = scales::oob_squish_any,
      limits = c(-1, 1), labels = c(NA, NA)) +
    ylab(NULL) +
    scale_x_continuous(expand = expansion(c(0.01, 0)), limits = c(start(ranges), end(ranges))) +
    theme_bw() +
    build_theme(params)
  if (!is.null(params$aspect.ratio)) {
    p <- p + theme(aspect.ratio = params$aspect.ratio)
  }
  p
}

