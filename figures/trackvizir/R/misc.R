norm_fun <- function(x, newMin = 0.05, newMax = 0.95) {
  if (length(unique(x)) == 1) {
    if (unique(x) < 0) newMin else newMax
  } else {
    ((x - min(x)) / (max(x) - min(x))) * (newMax - newMin) + newMin
  }
}

bin_ranges <- function(xranges, nbins = NULL, bin.size = NULL, where = "trackvizir") {
  if (!is.null(nbins) && !is.null(bin.size))
    stop("[", where, "] Please only provide one of nbins, bin.size")
  else if (is.null(nbins) && is.null(bin.size)) {
    message("[", where, "] Both nbins and bin.size are NULL, using bin.size = 1")
    bin.size <- 1
  }
  if (!is.null(nbins))
    xranges <- unlist(tile(xranges, n = nbins))
  else if (!is.null(bin.size))
    xranges <- unlist(tile(xranges, width = bin.size))
  xranges
}

# binnedSum2 <- function(bins, numvar, varname, na.rm = FALSE) {
#   binned <- binnedAverage(bins, numvar, varname, na.rm)
#   mcols(binned)[["varname"]] <- width(binned) * mcols(binned)[["varname"]]
#   binned
# }

binnedSum <- function(bins, numvar, varname, na.rm = FALSE) {
  if (!is(bins, "GRanges"))
    stop("'x' must be a GRanges object")
  if (!is(numvar, "RleList"))
    stop("'numvar' must be an RleList object")
  if (!identical(seqlevels(bins), names(numvar)))
    stop("'seqlevels(bins)' and 'names(numvar)' must be identical")
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar), function(seqname) {
    v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
    viewSums(v, na.rm = na.rm)
  })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[varname]] <- new_mcol
  bins
}

binnedMin <- function(bins, numvar, varname, na.rm = FALSE) {
  if (!is(bins, "GRanges"))
    stop("'x' must be a GRanges object")
  if (!is(numvar, "RleList"))
    stop("'numvar' must be an RleList object")
  if (!identical(seqlevels(bins), names(numvar)))
    stop("'seqlevels(bins)' and 'names(numvar)' must be identical")
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar), function(seqname) {
    v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
    viewMins(v, na.rm = na.rm)
  })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[varname]] <- new_mcol
  bins
}

binnedMax <- function(bins, numvar, varname, na.rm = FALSE) {
  if (!is(bins, "GRanges"))
    stop("'x' must be a GRanges object")
  if (!is(numvar, "RleList"))
    stop("'numvar' must be an RleList object")
  if (!identical(seqlevels(bins), names(numvar)))
    stop("'seqlevels(bins)' and 'names(numvar)' must be identical")
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar), function(seqname) {
    v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
    viewMaxs(v, na.rm = na.rm)
  })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[varname]] <- new_mcol
  bins
}

bins_fun <- function(bins, data, score.col = NULL, bin.strat = "avg",
  where = "trackvizir", ignore.strand = TRUE) {
  fixed <- fix_seqlevels(bins, data)
  bins <- fixed$ranges
  data <- fixed$data
  if (is.null(score.col)) {
    score.col <- colnames(mcols(data))
    if (!length(score.col))
      stop("[", where, "] Need a data column to compute bin data, couldn't find any")
    else {
      if (length(score.col) > 1) {
        message("[", where,
          "] Found more than one metadata column in GRanges object, using first one",
          " (set score.col to use a specific column)")
        score.col <- score.col[1]
      }
    }
  }
  if (ignore.strand) {
    binned <- switch(bin.strat,
      avg = binnedAverage(bins, coverage(data, weight = score.col), score.col),
      sum = binnedSum(bins, coverage(data, weight = score.col), score.col),
      max = binnedMax(bins, coverage(data, weight = score.col), score.col),
      min = binnedMin(bins, coverage(data, weight = score.col), score.col),
      stop("[", where, "] Unknown bin.strat")
    )
    strand(binned) <- "+"
    binned
  } else {
    data_p <- data[as.character(strand(data)) %in% c("+", "*")]
    data_m <- data[as.character(strand(data)) %in% c("-", "*")]
    binned_p <- switch(bin.strat,
      avg = binnedAverage(bins, coverage(data_p, weight = score.col), score.col),
      sum = binnedSum(bins, coverage(data_p, weight = score.col), score.col),
      max = binnedMax(bins, coverage(data_p, weight = score.col), score.col),
      min = binnedMin(bins, coverage(data_p, weight = score.col), score.col),
      stop("[", where, "] Unknown bin.strat")
    )
    strand(binned_p) <- "+"
    binned_m <- switch(bin.strat,
      avg = binnedAverage(bins, coverage(data_m, weight = score.col), score.col),
      sum = binnedSum(bins, coverage(data_m, weight = score.col), score.col),
      max = binnedMax(bins, coverage(data_m, weight = score.col), score.col),
      min = binnedMin(bins, coverage(data_m, weight = score.col), score.col),
      stop("[", where, "] Unknown bin.strat")
    )
    strand(binned_m) <- "-"
    sort(c(binned_p, binned_m))
  }
}

# bins_avg <- function(bins, data, score.col = NULL, where = "trackvizir") {
#   fixed <- fix_seqlevels(bins, data)
#   bins <- fixed$ranges
#   data <- fixed$data
#   if (is.null(score.col)) {
#     score.col <- colnames(mcols(data))
#     if (!length(score.col))
#       stop("[", where, "] Need a data column to compute bin data, couldn't find any")
#     else {
#       if (length(score.col) > 1) {
#         message("[", where,
#           "] Found more than one metadata column in GRanges object, using first one",
#           " (set score.col to use a specific column)")
#         score.col <- score.col[1]
#       }
#     }
#   }
#   binnedAverage(bins, coverage(data, weight = score.col), score.col)
# }

bins_count <- function(bins, data) {
  fixed <- fix_seqlevels(bins, data)
  bins <- fixed$ranges
  data <- fixed$data
  bins$count <- countOverlaps(bins, data)
  bins
}

fix_data_by_strand <- function(data, strand) {
  data$y[strand == "-"] <- data$y[strand == "-"] * -1
  if (any(strand == "*")) {
    data2 <- data[strand == "*", ]
    data2$y <- data2$y * -1
    data <- rbind(data, data2)
  }
  data
}

fix_seqlevels <- function(ranges, data) {
  xlvls <- seqlevels(ranges)
  ylvls <- seqlevels(data)
  if (all(!xlvls %in% ylvls))
    stop(
      paste0("Found no overlapping seqnames between ranges and data\n",
        "- seqlevels(ranges): ", paste0(xlvls, collapse = ", "), "\n",
        "- seqlevels(data): ", paste0(ylvls, collapse = ", "), "\n"),
      call. = FALSE)
  if (any(!xlvls %in% ylvls))
    warning(
      paste0("Some seqnames in ranges are absent from data\n",
        "- seqlevels(ranges): ", paste0(xlvls, collapse = ", "), "\n",
        "- seqlevels(data): ", paste0(ylvls, collapse = ", "), "\n"),
      call. = FALSE, immediate. = TRUE)
  lvls <- sort(unique(c(xlvls, ylvls)))
  seqlevels(ranges) <- lvls
  seqlevels(data) <- lvls
  list(ranges = ranges, data = data)
}

seqlengths2GRanges <- function(x) {
  sl <- seqlengths(x)
  if (anyNA(sl))
    stop("seqlengths(x) cannot contain NA values")
  x <- GRanges(names(sl), IRanges(1, unname(sl)))
  seqlengths(x) <- sl
  x
}

coords2GRanges <- function(chr, start, end) {
  GRanges(chr, IRanges(start, end))
}

TxDb2GRanges <- function(x) {
  txdb <- x
  x <- as.list(x)
  x$transcripts <- makeGRangesFromDataFrame(x$transcripts,
    keep.extra.columns = TRUE, seqnames.field = "tx_chrom",
    start.field = "tx_start", end.field = "tx_end",
    strand.field = "tx_strand")
  x$transcripts$annotation = "Transcript"
  x$transcripts$gene_name <- x$genes$gene_id[x$transcripts$tx_id]
  x$transcripts$exon_rank <- NA_integer_
  x$transcripts$exon_id <- NA_integer_
  x$transcripts$exon_name <- NA_character_
  x$transcripts$cds_id <- NA_integer_
  x$transcripts$cds_start <- NA_integer_
  x$transcripts$cds_end <- NA_integer_
  x$splicings <- makeGRangesFromDataFrame(x$splicings,
    keep.extra.columns = TRUE, seqnames.field = "exon_chrom",
    start.field = "exon_start", end.field = "exon_end",
    strand.field = "exon_strand")
  x$splicings$tx_name <- x$transcripts$tx_name[x$splicings$tx_id]
  x$splicings$tx_type <- x$transcripts$tx_type[x$splicings$tx_id]
  x$splicings$gene_name <- x$transcripts$gene_name[x$splicings$tx_id]
  x$splicings$annotation <- "Region"
  x <- c(x$transcripts, x$splicings)
  seqlevels(x) <- seqlevels(txdb)
  seqinfo(x) <- seqinfo(txdb)
  x
}

get_ranges <- function(data, local.ranges, global.ranges, prune.seqlevels = FALSE) {
  if (!is.null(local.ranges))
    ranges <- local.ranges
  else if (!is.null(global.ranges))
    ranges <- global.ranges
  else
    ranges <- calc_ranges_from_data(data)
  ranges <- ranges[, 0]
  if (prune.seqlevels)
    ranges <- prune_seqlevels(ranges)
  if (length(ranges) > 1)
    stop("Currently only single ranges are supported")
  sort(ranges)
}

prune_seqlevels <- function(x) {
  seqlevels(x) <- sort(unique(as.character(seqnames(x))))
  x
}

calc_ranges_from_data <- function(data, prune.seqlevels = FALSE) {
  if (prune.seqlevels)
    data <- prune_seqlevels(data)
  if (anyNA(seqlengths(data)))
    stop("To compute ranges from data, NA values cannot be present in seqlengths(data)")
  seqlengths2GRanges(data)
}
