prep_track <- function(args, data, trackType) {
  params <- check_params(collect_args(args), paste0("track_", trackType))
  ranges <- check_ranges(args$ranges, paste0("track_", trackType))
  if (!trackType %in% c("coords", "anno")) {
    data <- check_data(data, trackType, trackType == "vspace")
  }
  structure(list(Track.type = trackType, Ranges = ranges, Params = params, Data = data),
    class = "trackvizir.trackData")
}

check_ranges <- function(r, where = "trackvizir") {
  if (!is(r, "GRanges") && !is.null(r))
    stop("[", where, "] Expected NULL or a GRanges object for ranges", call. = FALSE)
  if (!is.null(r) && !length(r))
    warning("[", where, "] Found empty ranges, using NULL", immediate. = TRUE, call. = FALSE)
  if (length(r) > 1)
    stop("[", where, "] trackvizir currently only supports plotting a single range")
  r
}

check_data <- function(d, where = "trackvizir", is.vspace = FALSE) {
  if (!is.vspace) {
    if (!is(d, "GRanges"))
      stop("[", where, "] Expected a GRanges object for data", call. = FALSE)
    if (!length(d))
      stop("[", where, "] Found empty GRanges object for data", call. = FALSE)
  }
  d
}

build_tracks <- function(x) {
  if (!is(x, "trackvizir.trackGroup"))
    stop("Expected a trackvizir.trackGroup object")
  structure(build_plot(x), class = "trackvizir.trackGroupPlot")
}

build_plot <- function(x) {
  p <- lapply(x$Tracks, function(y) build_track(y, x$Params, x$Ranges))
  trackNames <- lapply(p, function(y) y$trackName)
  trackPlots <- lapply(p, function(y) y$trackPlot)
  h <- vapply(p, function(y) y$params$rel.size, numeric(1))
  highlights <- build_highlight(lapply(p, function(y) y$highlight), sum(h))
  list(Plots = trackPlots, Heights = h, Track.names = trackNames, Highlights = highlights)
}

build_trackname_grob <- function(track.name, angle = 0, size = 10,
  hjust = 1, vjust = 0.5, rel.x = 1, colour = "black", fontface = "plain",
  fontfamily = "sans") {
  if (is.null(track.name)) {
    NULL
  } else {
    textGrob(track.name, hjust = hjust, vjust = vjust, x = unit(rel.x, "npc"),
      rot = angle, gp = gpar(col = colour, fontsize = size, fontface = fontface,
        fontfamily = fontfamily))
  }
}

build_track <- function(track, global.params, global.ranges) {
  params <- get_params(track$Params, global.params)
  if (params$y.axis.hide) {
    params$y.axis.line.hide <- TRUE
    params$y.ticks.labels.outside.hide <- TRUE
    params$y.ticks.lines.hide <- TRUE
  }
  if (track$Track.type == "vspace") {
    trackName <- NULL
    ranges <- global.ranges
  } else {
    ranges <- get_ranges(track$Data, track$Ranges, global.ranges, params$prune.seqlevels)
    trackName <- build_trackname_grob(
        params$track.name,
        params$track.name.angle,
        params$track.name.size,
        params$track.name.hjust,
        params$track.name.vjust,
        params$track.name.rel.x,
        params$track.name.colour,
        params$track.name.fontface,
        params$track.name.family
      )
  }
  p <- switch(track$Track.type,
    scalebar = build_scalebar(track$Data, ranges, params),
    coords = build_coords(track$Data, ranges, params),
    histo = build_histo(track$Data, ranges, params),
    anno = build_anno(track$Data, ranges, params),
    point = build_point(track$Data, ranges, params),
    line = build_line(track$Data, ranges, params),
    bar = build_bar(track$Data, ranges, params),
    lollipop = build_bar(track$Data, ranges, params),
    histo = build_histo(track$Data, ranges, params),
    density = build_density(track$Data, ranges, params),
    area = build_area(track$Data, ranges, params),
    hiC = build_hiC(track$Data, ranges, params),
    pileup = build_pileup(track$Data, ranges, params),
    sashimi = build_sashimi(track$Data, ranges, params),
    sequence = build_sequence(track$Data, ranges, params),
    vspace = NULL,
    stop("Unknown Track.type")
  )
  if (!is.null(params$extra.ggplot)) p <- p + params$extra.ggplot
  hl <- build_highlight_df(ranges, params)
  list(trackPlot = p, trackName = trackName, highlight = hl, params = params)
}

build_highlight_df <- function(ranges, params) {
  if (is.null(ranges) || is.null(params$highlight.ranges)) {
    return(c("rel.size" = params$rel.size))
  }
  x1 <- min(start(ranges))
  x2 <- max(end(ranges))
  x1_2 <- x2 - x1
  hlr <- params$highlight.ranges
  data.frame(
    ymin = 0, ymax = params$rel.size,
    xmin = (start(hlr) - x1) / x1_2,
    xmax = (end(hlr) - x1) / x1_2,
    colour = params$highlight.colour,
    alpha = params$highlight.alpha
  )
}

build_highlight <- function(highlightList, total.height) {
  # TODO: deal with multiple chromosomes/ranges
  h_height <- 0
  for (i in rev(seq_along(highlightList))) {
    if (is.data.frame(highlightList[[i]])) {
      highlightList[[i]]$row <- i
      tmp_ymin <- highlightList[[i]]$ymin[1]
      tmp_ymax <- highlightList[[i]]$ymax[1]
      tmp_h <- tmp_ymax - tmp_ymin
      tmp_ymin <- tmp_ymin + h_height
      tmp_ymax <- tmp_ymin + tmp_h
      h_height <- tmp_ymax
      highlightList[[i]]$ymin <- tmp_ymin
      highlightList[[i]]$ymax <- tmp_ymax
    } else {
      h_height <- h_height + highlightList[[i]]["rel.size"]
      highlightList[[i]] <- NULL
    }
  }
  hl <- do.call(rbind, highlightList)
  if (!is.null(hl)) {
    hl <- ggplot(hl, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
        fill = colour, alpha = alpha)) +
      scale_y_continuous(expand = c(0, 0), limits = c(0, total.height)) +
      # TODO: How to deal with different expansion (eg y.tick.inward=TRUE)?
      scale_x_continuous(expand = expansion(c(0.01, 0)), limits = c(0, 1)) +
      geom_rect() +
      scale_alpha_identity() +
      scale_fill_identity() +
      theme(
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin()
      )
  }
  hl
}

build_plot_shell <- function(p, final, params, ranges) {
  y.ticks.labels <- c(
    format(if (is.na(params$ylims[1])) min(final$y) else params$ylims[1]),
    format(if (is.na(params$ylims[2])) max(final$y) else params$ylims[2])
  )
  if (!is.na(params$y.ticks.labels[1])) y.ticks.labels[1] <- params$y.ticks.labels[1]
  if (!is.na(params$y.ticks.labels[2])) y.ticks.labels[2] <- params$y.ticks.labels[2]
  ybreaks <- c(min(final$y), max(final$y))
  if (!is.na(params$ylims[1])) ybreaks[1] <- params$ylims[1]
  if (!is.na(params$ylims[2])) ybreaks[2] <- params$ylims[2]
  p <- p +
    scale_y_continuous(expand = c(0, 0), breaks = ybreaks, oob = scales::oob_squish_any,
      limits = ybreaks, labels = if (params$y.ticks.labels.outside.hide) c(NA, NA)
        else y.ticks.labels) +
    ylab(NULL) +
    theme_bw() +
    build_theme(params)
  if (!params$y.ticks.labels.inside.hide) {
    ytotal2 <- if (is.na(params$ylims[2])) max(final$y) else params$ylims[2]
    ytotal1 <- if (is.na(params$ylims[1])) min(final$y) else params$ylims[1]
    ytotal <- ytotal2 - ytotal1
    xrng <- c(start(ranges)[1], end(ranges)[1])
    p <- p +
      annotate("text", #x = start(ranges)[1],
        x = (xrng[2] - xrng[1]) * params$y.ticks.labels.inside.rel.x + xrng[1],
        y = (ytotal * params$y.ticks.labels.inside.rel.y) + ytotal1,
        vjust = params$y.ticks.labels.inside.vjust,
        hjust = params$y.ticks.labels.inside.hjust,
        size = params$y.ticks.labels.inside.size,
        label = paste0("[",
          if (is.na(params$y.ticks.labels[1]))
            (if (is.na(params$ylims[1])) format(min(final$y)) else params$ylims[1])
          else
            params$y.ticks.labels[1],
          ", ",
          if (is.na(params$y.ticks.labels[2]))
            (if (is.na(params$ylims[2])) format(max(final$y)) else params$ylims[2])
          else
            params$y.ticks.labels[2],
          "]"))
  }
  if (!is.null(params$aspect.ratio))
    p <- p + theme(aspect.ratio = params$aspect.ratio)
  p
}

build_theme <- function(params) {
  with(params, {
    theme <- theme(
      axis.line.y.left = if (y.axis.line.hide) element_blank()
        else element_line(colour = "black", size = y.axis.line.size),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.length.y.left = if (y.ticks.lines.hide) unit(0, "npc")
        else y.ticks.lines.length,
      axis.ticks.x = element_blank(),
      axis.ticks.y.left = if (y.ticks.lines.hide) element_blank()
        else element_line(colour = "black", size = y.ticks.lines.size),
      axis.text.y.left = if (y.ticks.labels.outside.hide) element_blank()
        else element_text(colour = "black", size = y.ticks.labels.outside.size),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.grid = element_blank()
    )
    if (compact) theme <- theme + theme(plot.margin = margin())
    theme
  })
}

