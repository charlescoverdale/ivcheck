# make_slide_figures.R
# Generate slide-sized figures for the conference deck.
#
# The paper's figures (in ../rj/figures/) are rendered at 5.5 x 3.4 in
# for R Journal columns. That aspect ratio works on 16:9 slides but the
# hero figure on slide 10 benefits from a taller render (9 x 5 in).
#
# Usage:
#   Rscript make_slide_figures.R
#
# Output:
#   figures/*.pdf     (copied or re-rendered from ../rj/figures/)
#   figures/qrcode_paper.png

suppressPackageStartupMessages({
  library(ggplot2)
  library(qrcode)
})

# ------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------

fig_dir <- "figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Chart typography: Helvetica sans-serif (matches the paper)
theme_slides <- function(base_size = 14) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      plot.title = element_blank(),                  # caption-led
      plot.subtitle = element_blank(),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      legend.position = "top",
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Okabe-Ito palette (matches the paper)
okabe_ito <- c(
  "#0072B2", "#E69F00", "#009E73", "#CC79A7",
  "#D55E00", "#56B4E9", "#F0E442", "#000000"
)

# ------------------------------------------------------------------
# Re-render the hero figure at slide size
# ------------------------------------------------------------------
# Edit this block for each package. The hero figure should be the one
# that lands on slide 10.
#
# Example: source the paper's figure script with a flag that re-renders
# only the hero figure at 9 x 5 inches instead of 5.5 x 3.4 inches.

# source("../rj/make_figures.R")          # re-use paper data
# render_hero_figure(                     # package-specific
#   output_path = file.path(fig_dir, "hero_figure.pdf"),
#   width = 9, height = 5
# )

# ------------------------------------------------------------------
# QR code to the paper PDF on the publications page
# ------------------------------------------------------------------

paper_url <- "https://charlescoverdale.github.io/files/coverdale_PACKAGE_2026.pdf"

qr <- qr_code(paper_url, ecl = "M")
png(
  filename = file.path(fig_dir, "qrcode_paper.png"),
  width = 800, height = 800, res = 300, bg = "white"
)
par(mar = rep(0, 4))
plot(qr)
dev.off()

cat("Slide figures written to", fig_dir, "\n")
