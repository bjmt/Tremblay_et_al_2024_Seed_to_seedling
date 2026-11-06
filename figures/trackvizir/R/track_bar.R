track_bar <- function(data, track.name = NULL, ...) {
  # Also useful for DNA methylation data, using a small 'width'
  prep_track(c(as.list(environment()), list(...)), data, "bar")
}

build_bar <- function(data, ranges, params) {
  if (is.null(params$nbins) && is.null(params$bin.size)) {
    binScores <- subsetByOverlaps(data, ranges, type = "within")
    bins <- binScores[, 0]
    binScores <- as.data.frame(binScores)
  } else {
    bins <- bin_ranges(ranges, params$nbins, params$bin.size, "track_bar")
    binScores <- subsetByOverlaps(data, ranges, type = "within")
    binScores <- bins_fun(bins, binScores, params$score.col, params$bin.strat, "track_bar",
      params$ignore.strand)
    binScores <- as.data.frame(binScores)
  }
  final <- data.frame(
    x = (binScores$start + binScores$end) / 2,
    y = binScores[[6]],
    strand = as.character(binScores$strand),
    width = width(bins)
  )
  if (!params$ignore.strand) {
    final <- fix_data_by_strand(final, final$strand)
  }
  if (!is.null(params$bin.smooth.k)) {
    final$y[final$strand == "+"] <- RcppRoll::roll_mean(final$y[final$strand == "+"],
      n = params$bin.smooth.k, fill = NA)
    final$y[final$strand == "-"] <- RcppRoll::roll_mean(final$y[final$strand == "-"],
      n = params$bin.smooth.k, fill = NA)
  }
  final <- final[final$y != 0, ]
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
    geom_bar(aes(fill = y < 0), stat = "identity", na.rm = TRUE,
      width = final$width, alpha = params$alpha) +
    scale_fill_manual(values = barcols, guide = "none") +
    scale_x_continuous(limits = c(start(ranges)[1], end(ranges)[1]),
      # TODO: work out what to do about different expansions and aligning plots
      # expand = expansion(c(if (!params$y.ticks.inside.hide) 0.02 else 0.01, 0)))
      expand = expansion(c(0.01, 0)))
  build_plot_shell(p, final, params, ranges)
}

