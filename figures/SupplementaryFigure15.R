# Note: the original code for this figure was lost, this is a recreation to reproduce the figure near exactly
#
# Supplementary Figure 15
#
# Package requirements: ggplot2, readr

library(ggplot2)

dir.create("SupplementaryFigure15", showWarnings = FALSE)

# Inputs -----------------------------------------------------------------------

# Input-normalized RNAPII ChIP-qPCR enrichment, relative to the promoter primers
ACT7_CHIP <- "RNAPII_ChIP_qPCR_ACT7.tsv"
DOG1_CHIP <- "RNAPII_ChIP_qPCR_DOG1.tsv"

# Functions --------------------------------------------------------------------

# One-sided Student's t-tests against IGN5, in the direction of the difference
test_vs_ign5 <- function(x, primers) {
  do.call(rbind, lapply(c("DS", "S72"), function(s) {
    xs <- x[startsWith(x$Sample, s), ]
    do.call(rbind, lapply(primers, function(p) {
      a <- xs[[p]]
      b <- xs$IGN5_SetI
      data.frame(Sample = s, Primer = p, P = t.test(a, b, var.equal = TRUE,
        alternative = ifelse(mean(a) > mean(b), "greater", "less"))$p.value)
    }))
  }))
}

plot_chip <- function(x, primers, ylab) {
  xl <- reshape(as.data.frame(x), direction = "long", varying = primers, v.names = "FC",
    timevar = "Primer", times = primers, idvar = "Sample")
  xl$Rep <- xl$Sample
  xl$Sample <- sub("[.]\\d+$", "", xl$Sample)
  xl$Primer <- factor(xl$Primer, levels = primers)
  xs <- aggregate(FC ~ Sample + Primer, xl, function(v) c(Mean = mean(v), SEM = sd(v) / sqrt(length(v))))
  xs <- data.frame(Sample = xs$Sample, Primer = xs$Primer, Mean = xs$FC[, "Mean"], SEM = xs$FC[, "SEM"])
  pv <- test_vs_ign5(x, primers[!primers %in% c("IGN5_SetI", primers[2])])
  pv$Label <- ifelse(pv$P < 0.01, "**", ifelse(pv$P < 0.05, "*", ifelse(pv$P < 0.1, ".", "")))
  xs <- merge(xs, pv, all.x = TRUE)
  xs$Label[is.na(xs$Label)] <- ""
  print(pv)
  ggplot(xs, aes(Primer, Mean)) +
    geom_col(fill = "grey35", width = 0.75) +
    geom_errorbar(aes(ymin = Mean - SEM, ymax = Mean + SEM), width = 0.2) +
    geom_point(data = xl, aes(Primer, FC), size = 0.3, colour = "grey60") +
    geom_text(aes(y = Mean + SEM, label = Label), vjust = -0.2, size = 5) +
    facet_wrap(~Sample, nrow = 1) +
    xlab(NULL) +
    ylab(ylab) +
    theme_bw() +
    theme(
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black"),
      strip.background = element_blank(),
      strip.text = element_text(size = 11)
    )
}

# Supplementary Figure 15a -----------------------------------------------------

act7 <- readr::read_tsv(ACT7_CHIP)

plot_chip(act7, c("IGN5_SetI", "Act7_-995", "Act7_+55", "Act7_+886", "Act7_+2561"),
  "FC to Act7_-995")
ggsave("SupplementaryFigure15/SuppFig15a_ACT7.pdf", width = 3.2, height = 2.6)

# Supplementary Figure 15b -----------------------------------------------------

dog1 <- readr::read_tsv(DOG1_CHIP)

plot_chip(dog1, c("IGN5_SetI", "DOG1_P3", "DOG1_C1", "DOG1_C2", "DOG1_C3"),
  "FC to DOG1_P3")
ggsave("SupplementaryFigure15/SuppFig15b_DOG1.pdf", width = 3.2, height = 2.6)
