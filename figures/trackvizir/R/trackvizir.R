trackvizir <- function(ranges = NULL,
  alpha = 1,
  anno.buffer = 0.001,
  anno.exon.size = 3,
  anno.gene.size = 1,
  anno.text.hide = FALSE,
  anno.text.size = 2,
  anno.text.colour = "black",
  anno.text.fontface = "plain",
  anno.text.family = "sans",
  aspect.ratio = NULL,
  bin.size = NULL,
  bin.smooth.k = NULL,
  ## bin.strat=sum: the tracks own bins are resized to 1
  bin.strat = "avg",
  # TODO: implement this
  colour.lut = NULL,
  colour.neg = "red3",
  colour.pos = "blue3",
  compact = FALSE,
  coord.text.size = 3,
  coord.text.colour = "black",
  coord.text.fontface = "plain",
  coord.text.family = "sans",
  equal.sized.chrs = FALSE,
  extra.ggplot = NULL,
  groups = NULL,
  groups.col = NULL,
  highlight.ranges = NULL,
  highlight.alpha = 0.3,
  highlight.colour = "khaki",
  # Maybe consider adding the ability to add an outline?
  # It can be done manually quite easily though:
  # x$Highlights <- x$Highlights + geom_rect(colour = "black", lty = 2)
  ignore.strand = FALSE,
  inter.chr.rel.size = 0.05,
  nbins = NULL,
  prune.seqlevels = FALSE,
  rel.size = 1,
  score.col = NULL,
  track.name = NULL,
  track.name.angle = 0,
  track.name.colour = "black",
  track.name.fontface = "plain",
  track.name.family = "sans",
  track.name.hjust = 1,
  track.name.rel.x = 1,
  track.name.size = 10,
  track.name.vjust = 0.5,
  # TODO: implement this
  # Apply this transformation to data by group/ranges
  transformation.fun = NULL,
  x.scale.expansion = c(0.01, 0),
  y.axis.hide = FALSE,
  y.axis.line.hide = FALSE,
  y.axis.line.size = 0.45,
  y.ticks.breaks = NULL,
  y.ticks.labels.inside.hide = TRUE,
  y.ticks.labels.inside.hjust = 0,
  y.ticks.labels.inside.size = 2,
  y.ticks.labels.inside.rel.x = 0,
  y.ticks.labels.inside.rel.y = 1,
  y.ticks.labels.inside.vjust = 1,
  y.ticks.labels.outside.hide = FALSE,
  y.ticks.labels.outside.size = 10,
  y.ticks.labels = c(NA, NA),
  y.ticks.lines.length = unit(0.25, "cm"),
  y.ticks.lines.hide = FALSE,
  y.ticks.lines.size = 0.25,
  ylims = c(NA, NA),
  zero.line = TRUE,
  zero.line.colour = "grey50",
  zero.line.size = 0.25
) {
  ranges <- check_ranges(ranges)
  if (!is.null(ranges)) ranges <- prune_seqlevels(ranges[, 0])
  if (prune.seqlevels && !is.null(ranges)) {
    message("Please note that prune.seqlevels=TRUE is ignored when ranges are provided")
    # stop("pruneseqlevels cannot be TRUE when ranges is not NULL")
    prune.seqlevels <- FALSE
  }
  params <- check_params(collect_args(as.list(environment())))
  tracks <- structure(list(), class = "trackvizir.trackDataList")
  structure(list(Ranges = ranges, Params = params, Tracks = tracks), class = "trackvizir.trackGroup")
}

DEFAULT_ARGS <- lapply(formals("trackvizir"), eval, envir = environment())

check_params <- function(x, where = "trackvizir") {
  badParams <- names(x)[!names(x) %in% names(DEFAULT_ARGS)]
  if (length(badParams)) {
    for (i in seq_along(badParams)) {
      warning("[", where, "] Ignoring the following unknown parameter: ",
        paste0(badParams[i], collapse = ", "),
        immediate. = TRUE, call. = FALSE)
    }
  }
  # TODO: add type checks here
  x
}

collect_args <- function(x) {
  x <- x[!vapply(x, is.null, logical(1))]
  x[!names(x) %in% c("data", "ranges")]
}

get_params <- function(local.params, global.params) {
  params <- DEFAULT_ARGS
  for (i in names(DEFAULT_ARGS)) {
    params[[i]] <- get_arg(i, params[[i]], local.params, global.params)
  }
  structure(params, class = "trackvizir.argList")
}

get_arg <- function(arg, default.param, local.params, global.params) {
  if (!is.null(local.params[[arg]]))
    local.params[[arg]]
  else if (!is.null(global.params[[arg]]))
    global.params[[arg]]
  else
    default.param
}

trackvizir_defaults <- function() {
  for (i in seq_along(DEFAULT_ARGS))
    cat(names(DEFAULT_ARGS)[i], ": ",
      if (is.null(DEFAULT_ARGS[[i]])) "NULL"
      else paste0(DEFAULT_ARGS[[i]], collapse = ", "),
      "\n", sep = "")
  invisible(DEFAULT_ARGS)
}
