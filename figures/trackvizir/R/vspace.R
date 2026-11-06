vspace <- function(rel.size = 0.2, ...) {
  if (abs(rel.size) > 1)
    stop("[vspace] rel.size cannot be less than -1 or greater than 1")
  prep_track(c(as.list(environment()), list(...)), NULL, "vspace")
}
