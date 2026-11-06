"+.trackvizir.trackGroup" <- function(x, y) {
  if (   !is(y, "trackvizir.trackData")
      && !is(y, "trackvizir.vspace")
     )
    stop("Expected a trackvizir.trackData or trackvizir.vspace object")
  x$Tracks <- structure(c(x$Tracks, list(y)), class = "trackvizir.trackDataList")
  structure(x, class = "trackvizir.trackGroup")
}

"-.trackvizir.trackGroup" <- function(x, y) {
  y <- as.integer(y)
  if (anyNA(y)) stop("NA values cannot be used")
  x$Tracks <- structure(x$Tracks[-y], class = "trackvizir.trackDataList")
  x
}

"+.trackvizir.trackGroupPlot" <- function(x, y) {
  stop("Tracks cannot be added to trackvizir.trackGroupPlot objects", call. = FALSE)
}

"-.trackvizir.trackGroupPlot" <- function(x, y) {
  stop("Tracks cannot be removed from trackvizir.trackGroupPlot objects", call. = FALSE)
}

"$.trackvizir.argList" <- function(x, y) {
  if (is.character(y) && !y %in% names(DEFAULT_ARGS)) {
    stop("[internal] Couldn't find ", y, " in DEFAULT_ARGS")
  }
  x[[y]]
}

print.trackvizir.trackGroup <- function(x, ...) {
  cat("[A trackvizir.trackGroup object]\n")
  cat("$Ranges\n")
  if (is.null(x$Ranges)) {
    cat("NULL\n")
  } else {
    rngs <- as.data.frame(x$Ranges)[, c("seqnames", "start", "end")]
    if (!nrow(rngs)) {
      cat("NULL\n")
    } else {
      for (i in 1:nrow(rngs)) {
        cat(as.character(rngs$seqnames[i]), ":",
          rngs$start[i], "..", rngs$end[i], "\n",  sep = "")
        if (i > 10) {
          cat("<...>\n")
          break
        }
      }
    }
  }
  cat("$Params\n")
  if (!length(x$Params)) {
    cat("NULL\n")
  } else {
    cat("<...>\n")
  }
  cat("$Tracks\n")
  print(x$Tracks, printHeader = FALSE)
  cat("\nTo build a trackvizir.trackGroupPlot object, run build_tracks(thisobject).\n")
  cat("To plot the tracks, run plot_tracks(thisobject).\n")
  cat("To save the tracks as an image, run save_tracks(thisobject).\n")
  invisible(x)
}

print.trackvizir.trackData <- function(x, ...) {
  cat("[A trackvizir.trackData object]\n")
  cat("$Track.type\n")
  cat(x$Track.type, "\n", sep = "")
  cat("$Ranges\n")
  if (is.null(x$Ranges)) {
    cat("NULL\n")
  } else {
    rngs <- as.data.frame(x$Ranges)[, c("seqnames", "start", "end")]
    if (!nrow(rngs)) {
      cat("NULL\n")
    } else {
      if (nrow(rngs) > 10) {
        print(rngs[1:10, ])
        cat("<...>\n")
      } else {
        print(rngs)
      }
    }
  }
  cat("$Params\n")
  if (!length(x$Params)) {
    cat("  NULL\n")
  } else {
    cat("<...>\n")
  }
  cat("$Data\n")
  if (is.null(x$Data)) {
    cat("NULL\n")
  } else {
    cat("<", class(x$Data)[1], ">\n", sep = "")
  }
  invisible(x)
}

print.trackvizir.trackDataList <- function(x, printHeader = TRUE, ...) {
  if (printHeader) cat("[A trackvizir.trackDataList object]\n")
  if (!length(x) || is.null(x)) {
    cat("<none>")
  } else {
    trackTypes <- lapply(x, function(x) x$Track.type)
    trackTypesCat <- paste0("[[", seq_along(trackTypes), "]] <trackvizir.trackData: ",
        trackTypes, ">")
    trackTypesCat[trackTypes == "vspace"] <- paste0("[[",
      which(trackTypes == "vspace"), "]] <trackvizir.vspace: ",
      lapply(x[trackTypes == "vspace"], function(x) x$Params$rel.size), ">")
    cat(trackTypesCat, sep = "\n")
  }
  invisible(x)
}

print.trackvizir.trackGroupPlot <- function(x, ...) {
  cat("[A trackvizir.trackGroupPlot object]\n")
  cat("$Plots\n")
  if (is.null(x$Plots)) {
    cat("NULL\n")
  } else {
    for (i in seq_along(x$Plots)) {
      if (is.null(x$Plots[[i]])) {
        cat("[[", i, "]] NULL\n", sep = "")
      } else {
        cat("[[", i, "]] <ggplot>\n", sep = "")
      }
      if (i > 10) {
        cat("<...>\n")
        break
      }
    }
  }
  cat("$Heights\n")
  if (is.null(x$Heights)) {
    cat("NULL\n")
  } else if (length(x$Heights)) {
    h <- x$Heights
    if (length(h) > 10) {
      h <- c(as.character(h[1:10]), "...")
    }
    cat(paste0(h, collapse = ", "), "\n", sep = "")
  }
  cat("$Track.names\n")
  if (is.null(x$Track.names)) {
    cat("NULL\n")
  } else {
    for (i in seq_along(x$Track.names)) {
      if (is.null(x$Track.names[[i]])) {
        cat("[[", i, "]] NULL\n", sep = "")
      } else {
        cat("[[", i, "]] <", as.character(x$Track.names[[i]]), ">\n", sep = "")
      }
      if (i > 10) {
        cat("<...>\n")
        break
      }
    }
  }
  cat("$Highlights\n")
  if (is.null(x$Highlights)) {
    cat("NULL\n")
  } else if (length(x$Highlights)) {
    cat("<ggplot>\n")
  }
  cat("\nTo plot the tracks, run plot_tracks(thisobject).\n")
  cat("To save the tracks as an image, run save_tracks(thisobject).\n")
  invisible(x)
}
