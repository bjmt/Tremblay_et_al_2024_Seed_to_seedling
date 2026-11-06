track_coords <- function(data = NULL, track.name = NULL, rel.size = 0.5,
  zero.line.size = 1, zero.line.colour = "black", compact = TRUE, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "coords")
}

build_coords <- function(data, ranges, params) {
  params$y.axis.line.hide <- TRUE
  params$y.ticks.lines.hide <- TRUE
  params$y.ticks.labels.outside.hide <- TRUE
  p <- ggplot() +
    geom_segment(aes(y = 1, yend = 1, x = 1, xend = 2),
      size = params$zero.line.size, colour = params$zero.line.colour) +
    geom_segment(aes(y = 0.5, yend = 1, x = 1, xend = 1),
      size = params$zero.line.size/2, colour = params$zero.line.colour) +
    geom_segment(aes(y = 0.5, yend = 1, x = 2, xend = 2),
      size = params$zero.line.size/2, colour = params$zero.line.colour) +
    geom_segment(aes(y = 0.5, yend = 1, x = 1.5, xend = 1.5),
      size = params$zero.line.size/2, colour = params$zero.line.colour) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-1, 1),# oob = scales::oob_squish_any,
      limits = c(-1, 1), labels = c(NA, NA)) +
    annotate("text", x = 1, y = 0, vjust = 1, hjust = 0,
      size = params$coord.text.size, colour = params$coord.text.colour,
      family = params$coord.text.family, fontface = params$coord.text.fontface,
      label = format(start(ranges), scientific = FALSE, big.mark = ",")) +
    annotate("text", x = 1.5, y = 0, vjust = 1, hjust = 0.5, 
      size = params$coord.text.size, colour = params$coord.text.colour,
      family = params$coord.text.family, fontface = params$coord.text.fontface,
      label = format(round(mean(c(end(ranges), start(ranges)))), scientific = FALSE, big.mark = ",")) +
    annotate("text", x = 2, y = 0, vjust = 1, hjust = 1,
      size = params$coord.text.size, colour = params$coord.text.colour,
      family = params$coord.text.family, fontface = params$coord.text.fontface,
      label = format(end(ranges), scientific = FALSE, big.mark = ",")) +
    ylab(NULL) +
    scale_x_continuous(expand = expansion(c(0.01, 0)), limits = c(1, 2)) +
    theme_bw() +
    build_theme(params)
  if (!is.null(params$aspect.ratio)) {
    p <- p + theme(aspect.ratio = params$aspect.ratio)
  }
  p
}

