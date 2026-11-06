widthDetailsWithGP <- function (grob) {
    if (!is.null(grob$gp)) {
        pushViewport(viewport(gp=grob$gp))
        on.exit(popViewport())
    }
    return(widthDetails(grob))
}

get_textGrob_width <- function(x) {
  # The point of this function is to deal with cases other than:
  # - track.name.rel.x=1,track.name.hjust=1
  # - track.name.rel.x=0,track.name.hjust=0
  # In these cases, the widthDetails function does not include the extra whitespace.

  # TODO: this is assuming hjust=1!
  # hjust <- x$hjust
  # x <- as.numeric(convertUnit(x$x, "npc"))
  # w <- as.numeric(widthDetails(x))
  # w * (1 / x)

  # TODO: control space between track.names column and plots column
  #   --> just rel.x to control this

  as.numeric(widthDetails(x))
}

plot_tracks <- function(x, track.name.rel.width = "auto",
  track.name.position = c("left", "right"), highlight.ontop = TRUE
) {
  track.name.position <- match.arg(track.name.position)
  if (is(x, "trackvizir.trackGroup"))
    trackPlots <- build_plot(x)
  else if (is(x, "trackvizir.trackGroupPlot"))
    trackPlots <- x
  else
    stop("Expected a trackvizir.trackGroup or trackvizir.trackGroupPlot obejct")
  if (!length(trackPlots$Plots))
    stop("This trackvizir.trackGroupPlot object does not contain any tracks")

  if (track.name.rel.width == "auto") {
    # TODO: take into account rel.x + hjust!
    if (all(vapply(trackPlots$Track.names, is.null, logical(1)))) {
      track.name.rel.width <- 0
    } else {
      if (is.null(dev.list())) {
        message("track.name.rel.width=\"auto\" only works when a plot device is already open!")
        message("Try running the function a second time after the plot opens.")
        track.name.rel.width <- 0.2
      } else {
        pWidth <- dev.size("in")[1]
        nWidths <- vapply(trackPlots$Track.names, get_textGrob_width, numeric(1))
        # For some reason the widthDetails call is wrong when plot_tracks
        # is called when no device is open!
        track.name.rel.width <-  max(nWidths, na.rm = TRUE) / pWidth
      }
    }
    message("Using track.name.rel.width = ",
      track.name.rel.width * 1 / (1 - track.name.rel.width))
    if (track.name.position == "left") {
      track.name.rel.width <- c(track.name.rel.width, 1 - track.name.rel.width)
    } else if (track.name.position == "right") {
      track.name.rel.width <- c(1 - track.name.rel.width, track.name.rel.width)
    }
  } else if (track.name.position == "left") {
    track.name.rel.width <- c(track.name.rel.width, 1)
  } else if (track.name.position == "right") {
    track.name.rel.width <- c(1, track.name.rel.width)
  }

  trackPlots$Track.names <- cowplot::plot_grid(
    plotlist = trackPlots$Track.names, rel_heights = trackPlots$Heights,
    ncol = 1, align = "v", axis = "lr")

  if (is.null(trackPlots$Highlights)) {
    trackPlots$Plots <- cowplot::plot_grid(
      plotlist = trackPlots$Plots, rel_heights = trackPlots$Heights,
      ncol = 1, align = "v", axis = "lr")
  } else {
    if (highlight.ontop)
      tP <- c(list(trackPlots$Highlights), trackPlots$Plots)
    else
      tP <- c(trackPlots$Plots, list(trackPlots$Highlights))
    tP <- cowplot::align_plots(plotlist = tP, align = "v", axis = "lr")
    trackPlots$Plots <- cowplot::plot_grid(
      plotlist = tP[-1], rel_heights = trackPlots$Heights,
      ncol = 1)
    trackPlots$Plots <- cowplot::ggdraw(trackPlots$Plots) +
      cowplot::draw_plot(tP[[1]])
  }

  if (track.name.position == "left") {
    trackPlotsList <- list(trackPlots$Track.names, trackPlots$Plots)
  } else if (track.name.position == "right") {
    trackPlotsList <- list(trackPlots$Plots, trackPlots$Track.names)
  }

  cowplot::plot_grid(plotlist = trackPlotsList, ncol = 2,
    rel_widths = track.name.rel.width)
}

save_tracks <- function(x, filename, track.name.rel.width = "auto",
  track.name.position = c("left", "right"), highlight.ontop = TRUE,
  base_height = 3, base_width = 9, base_asp = NULL, units = "in", ...
) {
  if (is.null(base_height) && is.null(base_width))
    stop("Both base_height and base_width cannot be NULL")
  if (is.null(base_height) && !is.null(base_width)) {
    if (is.null(base_asp))
      stop("base_asp cannot be NULL if one of base_height, base_width are NULL")
    base_height <- base_width / base_asp
  }
  else if (!is.null(base_height) && is.null(base_width)) {
    if (is.null(base_asp))
      stop("base_asp cannot be NULL if one of base_height, base_width are NULL")
    base_width <- base_height * base_asp
  }
  if (is(x, "trackvizir.trackGroup") || is(x, "trackvizir.trackGroupPlot")) {
    if (is(x, "trackvizir.trackGroup")) {
      x <- build_plot(x)
    }
    if (track.name.rel.width == "auto") {
      if (all(vapply(x$Track.names, is.null, logical(1)))) {
        track.name.rel.width <- 0
      } else {
        nWidths <- max(vapply(x$Track.names, get_textGrob_width, numeric(1)), na.rm = TRUE)
        nWidths <- as.numeric(convertUnit(unit(nWidths, "in"), units[1]))
        track.name.rel.width <-  nWidths / base_width
        track.name.rel.width <- track.name.rel.width * 1 / (1 - track.name.rel.width)
      }
      message("Using track.name.rel.width = ", track.name.rel.width)
    }
    x <- plot_tracks(x, track.name.rel.width, track.name.position, highlight.ontop)
  }
  else if (!is(x, "ggplot"))
    stop("Expected a trackvizir.trackGroup, trackvizir.trackGroupPlot, or ggplot object")
  cowplot::save_plot(filename = filename, plot = x, units = units,
    base_height = base_height, base_width = base_width, ...)
}
