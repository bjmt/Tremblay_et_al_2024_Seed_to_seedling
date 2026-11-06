# Add an extras option so users can supply their own ggplot commands to each
# track in the form of a function, eg: function() geom_hline(yintercept = 1)

track_scalebar <- function(data = NULL, track.name = NULL, rel.size = 0.2, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "scalebar")
}

track_point <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "point")
}

track_line <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "line")
}

track_lollipop <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "lollipop")
}

track_histo <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "histo")
}

track_density <- function(data, track.name = NULL, ...) {
  # based on colour, single height bar plot
  prep_track(c(as.list(environment()), list(...)), data, "density")
}

track_hiC <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "hiC")
}

track_pileup <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "pileup")
}

track_sashimi <- function(data, track.name = NULL, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "sashimi")
}

track_sequence <- function(data, track.name = NULL, rel.size = 0.2, ...) {
  prep_track(c(as.list(environment()), list(...)), data, "sequence")
}

build_scalebar <- function(data, ranges, params) {
  ggplot()
}

build_histo <- function(data, ranges, params) {
  ggplot()
}

build_point <- function(data, ranges, params) {
  ggplot()
}

build_line <- function(data, ranges, params) {
  ggplot()
}

build_lollipop <- function(data, ranges, params) {
  ggplot()
}

build_histo <- function(data, ranges, params) {
  ggplot()
}

build_density <- function(data, ranges, params) {
  ggplot()
}

build_hiC <- function(data, ranges, params) {
  ggplot()
}

build_pileup <- function(data, ranges, params) {
  ggplot()
}

build_sashimi <- function(data, ranges, params) {
  ggplot()
}

build_sequence <- function(data, ranges, params) {
  ggplot()
}
