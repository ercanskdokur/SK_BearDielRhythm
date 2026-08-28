# ==============================================================================
# 07c_fig01_replot.R  —  FigV3_01 (zero-calibration) replot ONLY
# ==============================================================================
source(file.path(Sys.getenv("PROJECT_ROOT", unset = "/path/to/project"), "BrownBearDielAct", "scripts", "00_setup.R"))
RS <- get_response_spec()
init_log(sprintf("07c_fig01_replot_%s", RS$key))
log_step("=== 07c_fig01_replot | RESP=%s ===", RS$key)

fp <- file.path(tbl_path, paste0(tag_out("Table_V3_02_PPC_Stats"), ".csv"))
if (!file.exists(fp)) stop("Missing: ", fp)
ppc <- utils::read.csv(fp, stringsAsFactors = FALSE)
zr <- ppc[ppc$Statistic == "prop_zero", ]
if (nrow(zr) == 0) stop("No prop_zero row")
zr$Model <- factor(zr$Model, levels = c("H0","H1","H2","H3","H4","H5"))
ok6 <- c(H0="#999999", H1="#E69F00", H2="#56B4E9",
         H3="#009E73", H4="#D55E00", H5="#CC79A7")

p <- ggplot2::ggplot(zr, ggplot2::aes(x = Observed, y = PostMean,
                                      colour = Model, label = Model)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  ggplot2::geom_point(size = 3.5) +
  ggplot2::scale_colour_manual(values = ok6, name = "Hypothesis") +
  ggplot2::labs(x = "Observed zero proportion", y = "Model-predicted zero proportion",
                title = "Zero-inflation calibration (correct prop_zero)",
                subtitle = "Near line = good fit; far = zi component under/over-estimated") +
  theme_bear()
if (requireNamespace("ggrepel", quietly = TRUE)) {
  p <- p + ggrepel::geom_text_repel(size = 3.4, show.legend = FALSE,
            max.overlaps = Inf, box.padding = 0.7, point.padding = 0.5,
            min.segment.length = 0, seed = 42, segment.colour = "grey55")
  log_step("ggrepel present -> labels spread with connector lines")
} else {
  p <- p + ggplot2::geom_text(vjust = -0.9, hjust = -0.2, size = 3.4, show.legend = FALSE)
  log_step("ggrepel absent -> colour+legend distinguishes; label-nudge fallback")
}
save_fig(p, tag_out("FigV3_01_PPC_zero_calibration"), w = 7.5, h = 6)
log_step("=== 07c_fig01_replot DONE: %s ===", tag_out("FigV3_01_PPC_zero_calibration"))
sink(type = "message"); sink()
