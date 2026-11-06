read_bw <- function(plus, minus = NULL, ranges = NULL, minusIsNegative = FALSE) {
  track <- rtracklayer::import.bw(plus, which = ranges)
  strand(track) <- "+"
  if (!is.null(minus)) {
    minus <- rtracklayer::import.bw(minus, which = ranges)
    strand(minus) <- "-"
    if (minusIsNegative) minus$score <- minus$score * -1
    track <- c(track, minus)
  }
  sort(track)
}

read_bg <- function(plus, minus = NULL, ranges = NULL, minusIsNegative = FALSE) {
  track <- rtracklayer::import.bedGraph(plus, which = ranges)
  strand(track) <- "+"
  if (!is.null(minus)) {
    minus <- rtracklayer::import.bedGraph(minus, which = ranges)
    if (minusIsNegative) minus$score <- minus$score * -1
    strand(minus) <- "-"
    track <- c(track, minus)
  }
  sort(track)
}

