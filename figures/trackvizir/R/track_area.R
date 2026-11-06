track_area <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "area")
}

build_area <- function(data, ranges, params) {
  if (is.null(params$nbins) && is.null(params$bin.size)) {
    binScores <- subsetByOverlaps(data, ranges, type = "within")
    bins <- binScores[, 0]
    binScores <- as.data.frame(binScores)
    final <- data.frame(
      x = c(rbind(binScores$start, binScores$end)),
      y = rep(binScores$score, each = 2),
      strand = rep(as.character(binScores$strand), each = 2)
    )
    final <- final[!duplicated(final$x), ]
  } else {
    bins <- bin_ranges(ranges, params$nbins, params$bin.size, "track_area")
    binScores <- subsetByOverlaps(data, ranges, type = "within")
    binScores <- bins_fun(bins, binScores, params$score.col, params$bin.strat, "track_area",
      params$ignore.strand)
    binScores <- as.data.frame(binScores)
    final <- data.frame(
      x = (binScores$start + binScores$end) / 2,
      y = binScores[[6]],
      strand = as.character(binScores$strand)
    )
  }
  if (!params$ignore.strand) {
    final <- fix_data_by_strand(final, final$strand)
  }
  if (!is.null(params$bin.smooth.k)) {
    final$y[final$strand == "+"] <- RcppRoll::roll_mean(final$y[final$strand == "+"],
      n = params$bin.smooth.k, fill = NA)
    final$y[final$strand == "-"] <- RcppRoll::roll_mean(final$y[final$strand == "-"],
      n = params$bin.smooth.k, fill = NA)
  }
  p <- ggplot(final, aes(x, y))
  if (params$zero.line) {
    p <- p +
      geom_segment(aes(y = 0, yend = 0, x = start(ranges)[1], xend = end(ranges)[1]),
        size = params$zero.line.size, colour = params$zero.line.colour)
  }
  barcols <- params$colour.pos
  if (is.null(params$colour.neg)) {
    barcols <- c(barcols, barcols)
  } else {
    barcols <- c(barcols, params$colour.neg)
  }
  p <- p +
    geom_area(data = function(d) d[d$y > 0, ],
      aes(fill = y < 0), na.rm = TRUE) +
    geom_area(data = function(d) d[d$y < 0, ],
      aes(fill = y < 0), na.rm = TRUE) +
    scale_fill_manual(values = barcols, guide = "none") +
    scale_x_continuous(limits = c(start(ranges)[1], end(ranges)[1]),
      # TODO: work out what to do about different expansions and aligning plots
      # expand = expansion(c(if (!params$y.ticks.inside.hide) 0.02 else 0.01, 0)))
      expand = expansion(c(0.01, 0)))
  build_plot_shell(p, final, params, ranges)
}

