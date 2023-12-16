# Benjamin J.M. Tremblay, June 2023
#
# The heatmaps shown in Figures 6d and 7a were created in part using
# these custom functions and ComplexHeatmap. Metagene plots such as
# Figures S1f-i, S3f-g, S5e, S6a-c, e-h, S7c-g, S9a-d, g-h, and S10b
# were also created using these functions in combination with the
# base plot() function.
#
# Package requirements: matrixStats, zoo, genomation and rtracklayer

nn0 <- function(x, k = 10, norm1 = FALSE) {
  x <- zoo::rollmean(x, k = k, fill = NA)
  if (norm1) x <- (x - min(x, na.rm = TRUE)) /
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  x[is.na(x)] <- 0
  x
}

read_signal <- function(pos, neg = NULL, windows = NULL, use.neg.sign = FALSE) {
  if (!is.null(windows)) {
    pos <- rtracklayer::import(pos, selection = windows)
    if (!is.null(neg)) {
      neg <- rtracklayer::import(neg, selection = windows)
    }
  } else {
    pos <- rtracklayer::import(pos)
    if (!is.null(neg)) {
      neg <- rtracklayer::import(neg)
    }
  }
  if (!is.null(neg) && use.neg.sign) {
    neg$score <- -abs(neg$score)
  } else if (!is.null(neg)) {
    neg$score <- abs(neg$score)
  }
  if (!is.null(neg)) {
    strand(pos) <- "+"
    strand(neg) <- "-"
    pos <- sort(c(pos, neg))
  }
  pos
}

get_unstranded_signal <- function(signal, windows, score.col = "score",
  keep = c(0, 1), return.matrix = FALSE, use.pctile = 0, row.norm = TRUE,
  merge.signal = FALSE) {

  sig <- genomation::ScoreMatrix(signal, windows, strand.aware = TRUE, weight.col = score.col)
  
  sig <- sig@.Data

  sig[sig < 0 | is.na(sig)] <- 0

  if (keep[1] != 0 || keep[2] != 1) {
    qL <- quantile(as.numeric(sig), p = keep, na.rm = TRUE)
    sig[sig < qL[1]] <- qL[1]
    sig[sig > qL[2]] <- qL[2]
  }

  if (use.pctile > 0) {
    thresh <- quantile(matrixStats::rowMaxs(sig, na.rm = TRUE) -
      matrixStats::rowMins(sig, na.rm = TRUE), p = use.pctile, na.rm = TRUE)
    sig <- sig[(rowMaxs(sig, na.rm = TRUE) - rowMins(sig, na.rm = TRUE)) >= thresh, ]
  }

  if (row.norm) {
    for (i in 1:nrow(sig)) {
      s_i <- max(sig[i, ], na.rm = TRUE)
      b_i <- min(sig[i, ], na.rm = TRUE)
      if (s_i == b_i) sig[i, ] <- 0
      else if (s_i > 0) sig[i, ] <- (sig[i, ] - b_i) / (s_i - b_i)
    }
  }

  if (!return.matrix) {

  sig_avg <- colMeans(sig, na.rm = TRUE)
  if (row.norm) {
    sig_avg <- (sig_avg - min(sig_avg, na.rm = TRUE)) /
      (max(sig_avg, na.rm = TRUE) - min(sig_avg, na.rm = TRUE))
  }

  data.frame(row.names = NULL,
    Coord = seq_along(sig_avg),
    Avg = sig_avg
  )

  } else {

    as.matrix(sig)

  }

}

get_stranded_signal <- function(signal, windows, score.col = "score",
  keep = c(0, 1), return.matrix = FALSE, use.pctile = 0, row.norm = TRUE,
  merge.signal = FALSE) {

  windows$order <- 1:length(windows)

  pos_win <- windows[as.character(strand(windows)) == "+"]
  neg_win <- windows[as.character(strand(windows)) == "-"]
  pos_sig <- signal[as.character(strand(signal)) == "+"]
  neg_sig <- signal[as.character(strand(signal)) == "-"]

  if (!length(pos_win) && !length(neg_win))
    stop("Couldn't find any windows on +/- strands")

  if (length(pos_win)) {
    pos_pos <- genomation::ScoreMatrix(pos_sig, pos_win, strand.aware = TRUE, weight.col = score.col)
    neg_pos <- genomation::ScoreMatrix(neg_sig, pos_win, strand.aware = TRUE, weight.col = score.col)
  } else {
    pos_pos <- matrix(nrow = 0, ncol = width(windows)[1])
    neg_pos <- matrix(nrow = 0, ncol = width(windows)[1])
  }
  if (length(neg_win)) {
    pos_neg <- genomation::ScoreMatrix(pos_sig, neg_win, strand.aware = TRUE, weight.col = score.col)
    neg_neg <- genomation::ScoreMatrix(neg_sig, neg_win, strand.aware = TRUE, weight.col = score.col)
  } else {
    pos_neg <- matrix(nrow = 0, ncol = width(windows)[1])
    neg_neg <- matrix(nrow = 0, ncol = width(windows)[1])
  }

  if (nrow(pos_pos) && nrow(neg_pos)) {
    good_pos <- dimnames(pos_pos)[[1]][dimnames(pos_pos)[[1]] %in% dimnames(neg_pos)[[1]]]
  } else {
    good_pos <- character()
  }
  if (nrow(pos_neg) && nrow(neg_neg)) {
    good_neg <- dimnames(pos_neg)[[1]][dimnames(pos_neg)[[1]] %in% dimnames(neg_neg)[[1]]]
  } else {
    good_neg <- character()
  }

  if (!length(good_pos) && !length(good_neg))
    stop("Couldn't recover any ranges")

  if (length(good_pos)) {
    pos_pos <- pos_pos[good_pos, ]@.Data
    neg_pos <- neg_pos[good_pos, ]@.Data
  }
  if (length(good_neg)) {
    pos_neg <- pos_neg[good_neg, ]@.Data
    neg_neg <- neg_neg[good_neg, ]@.Data
  }

  pos_pos[pos_pos < 0 | is.na(pos_pos)] <- 0
  pos_neg[pos_neg < 0 | is.na(pos_neg)] <- 0
  neg_pos[neg_pos < 0 | is.na(neg_pos)] <- 0
  neg_neg[neg_neg < 0 | is.na(neg_neg)] <- 0

  ## Trim:

  if (keep[1] != 0 || keep[2] != 1) {
    qL <- quantile(c(pos_pos, pos_neg, neg_pos, neg_neg), p = keep, na.rm = TRUE)
    pos_pos[pos_pos < qL[1]] <- qL[1]
    pos_neg[pos_neg < qL[1]] <- qL[1]
    neg_pos[neg_pos < qL[1]] <- qL[1]
    neg_neg[neg_neg < qL[1]] <- qL[1]
    pos_pos[pos_pos > qL[2]] <- qL[2]
    pos_neg[pos_neg > qL[2]] <- qL[2]
    neg_pos[neg_pos > qL[2]] <- qL[2]
    neg_neg[neg_neg > qL[2]] <- qL[2]
  }

  ## Only keep top percentile:

  if (use.pctile > 0) {
    thresh_sense <- quantile(
      c(matrixStats::rowMaxs(pos_pos, na.rm = TRUE) -
        matrixStats::rowMins(pos_pos, na.rm = TRUE),
        matrixStats::rowMaxs(neg_neg, na.rm = TRUE) -
          matrixStats::rowMins(neg_eg, na.rm = TRUE)),
      p = use.pctile, na.rm = TRUE)
    thresh_antisense <- quantile(
      c(matrixStats::rowMaxs(neg_pos, na.rm = TRUE) -
        matrixStats::rowMins(neg_pos, na.rm = TRUE),
        matrixStats::rowMaxs(pos_neg, na.rm = TRUE) -
          matrixStats::rowMins(pos_neg, na.rm = TRUE)),
      p = use.pctile, na.rm = TRUE)
    pos_pos <- pos_pos[(matrixStats::rowMaxs(pos_pos, na.rm = TRUE) -
      matrixStats::rowMins(pos_pos, na.rm = TRUE)) >= thresh_sense, ]
    neg_neg <- neg_neg[(matrixStats::rowMaxs(neg_neg, na.rm = TRUE) -
      matrixStats::rowMins(neg_neg, na.rm = TRUE)) >= thresh_sense, ]
    neg_pos <- neg_pos[(matrixStats::rowMaxs(neg_pos, na.rm = TRUE) -
      matrixStats::rowMins(neg_pos, na.rm = TRUE)) >= thresh_antisense, ]
    pos_neg <- pos_neg[(matrixStats::rowMaxs(pos_neg, na.rm = TRUE) -
      matrixStats::rowMins(pos_neg, na.rm = TRUE)) >= thresh_antisense, ]
  }

  ## Normalize signal between 0 and 1 per window:

  if (row.norm) {

    for (i in seq_len(nrow(pos_pos))) {
      s_i <- max(pos_pos[i, ], na.rm = TRUE)
      b_i <- min(pos_pos[i, ], na.rm = TRUE)
      if (s_i == b_i) pos_pos[i, ] <- 0
      else if (s_i > 0) pos_pos[i, ] <- (pos_pos[i, ] - b_i) / (s_i - b_i)
    }

    for (i in seq_len(nrow(neg_pos))) {
      s_i <- max(neg_pos[i, ], na.rm = TRUE)
      b_i <- min(neg_pos[i, ], na.rm = TRUE)
      if (s_i == b_i) neg_pos[i, ] <- 0
      else if (s_i > 0) neg_pos[i, ] <- (neg_pos[i, ] - b_i) / (s_i - b_i)
    }

    for (i in seq_len(nrow(pos_neg))) {
      s_i <- max(pos_neg[i, ], na.rm = TRUE)
      b_i <- min(pos_neg[i, ], na.rm = TRUE)
      if (s_i == b_i) pos_neg[i, ] <- 0
      else if (s_i > 0) pos_neg[i, ] <- (pos_neg[i, ] - b_i) / (s_i - b_i)
    }

    for (i in seq_len(nrow(neg_neg))) {
      s_i <- max(neg_neg[i, ], na.rm = TRUE)
      b_i <- min(neg_neg[i, ], na.rm = TRUE)
      if (s_i == b_i) neg_neg[i, ] <- 0
      else if (s_i > 0) neg_neg[i, ] <- (neg_neg[i, ] - b_i) / (s_i - b_i)
    }

  }

  if (!return.matrix) {

    pos_avg <- colMeans(rbind(pos_pos, neg_neg), na.rm = TRUE)
    neg_avg <- colMeans(rbind(pos_neg, neg_pos), na.rm = TRUE)

    ## Normalize to between 0 and 1 again for the average:

    if (row.norm) {
      if (any(pos_avg > 0))
        pos_avg <- (pos_avg - min(pos_avg, na.rm = TRUE)) /
          (max(pos_avg, na.rm = TRUE) - min(pos_avg, na.rm = TRUE))
      if (any(neg_avg > 0))
        neg_avg <- (neg_avg - min(neg_avg, na.rm = TRUE)) /
          (max(neg_avg, na.rm = TRUE) - min(neg_avg, na.rm = TRUE))
    }

    data.frame(row.names = NULL,
      Coord = seq_along(pos_avg),
      PosAvg = pos_avg,
      NegAvg = neg_avg,
      MergedAvg = pos_avg - neg_avg
    )

  } else {

    good_pos <- pos_win$order[as.integer(good_pos)]
    good_neg <- neg_win$order[as.integer(good_neg)]

    sense_m <- rbind(pos_pos, neg_neg)[order(c(good_pos, good_neg)), ]
    antisense_m <- rbind(neg_pos, pos_neg)[order(c(good_pos, good_neg)), ]

    if (!merge.signal) {
      list(Sense = sense_m, Antisense = antisense_m)
    } else {
      sense_m - antisense_m
    }

  }
  
}

