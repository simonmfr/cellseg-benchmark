plot_benchmark_heatmap <- function(df_funky, out_file, metric_labels, family = "Arial", font_file, score_label = "Overall", row_height = 1.02, width = NULL, height = NULL, note = NULL) {
  library(funkyheatmap)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)
  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  sysfonts::font_add(family, path.expand(font_file))
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 600)
  update_geom_defaults("text", list(family = family))
  for (s in c("sep0", "sep1", "sep2", "sep3")) df_funky[[s]] <- ""
  # width (inches): text sized for print at this width (6 pt; group headers 5, scores 5.5, densities 4)
  # height (inches): extra space between rows so the table prints at this height (funkyrect shapes distort with row_height)
  s <- 1
  fs <- function(native, pt = 6) if (is.null(width)) native else pt / .pt / s
  build <- function(name_w) {

  box <- list(width = 1.02)
  col_spec <- tribble(
    ~id,                          ~group,       ~geom,       ~palette,      ~name,          ~opts,
    "id",                         "",           "text",      NA,            "",             list(width = name_w, hjust = 0, size = fs(4.2)),
    "sep0",                       "",           "text",      NA,            "",             list(width = 0.5),
    "total_score",                "overall",    "bar",       "overall",     "",             list(width = 2.9, draw_outline = FALSE),
    "sep1",                       "",           "text",      NA,            "",             list(width = 0.3),
    "MECR_adj",                   "expression", "funkyrect", "expression",  "MECR",         box,
    "marker_f1",                  "expression", "funkyrect", "expression",  "Marker F1",    box,
    "silhouette_score",           "expression", "funkyrect", "expression",  "Silhouette",   box,
    "sep2",                       "",           "text",      NA,            "",             list(width = 0.3),
    "ficture_f1_score",           "spatial",    "funkyrect", "spatial",     "Ficture F1",   box,
    "low_z_integrity_vs_tissue",  "spatial",    "funkyrect", "spatial",     "Z integrity",  box,
    "sep3",                       "",           "text",      NA,            "",             list(width = 0.6),
    "pct_assigned_qced",          "general",    "funkyrect", "transcripts", "% Transcripts", box,
    "volume_final",               "general",    "funkyrect", "general",     "Volume",       box,
    "elongation",                 "general",    "funkyrect", "general",     "Elongation",   box,
    "intensities_DAPI",           "general",    "funkyrect", "general",     "DAPI",         box,
    "cells_per_mm2",              "general",    "rect",      "general_lt",  "Cells/mm²",    list(width = 1.1),
    "density_label",              "general",    "text",      NA,            "",             list(overlay = TRUE, size = fs(3.3, 4), hjust = 0.5),
  )

  column_info <- col_spec %>%
    filter(id %in% names(df_funky)) %>%
    left_join(metric_labels, by = "id", suffix = c("", "_label")) %>%
    mutate(name = coalesce(name_label, name), name = ifelse(id == "cells_per_mm2", "Cells/mm²\n(×10³)", name),
           options = opts, legend = FALSE) %>%
    select(-name_label, -opts)

  column_groups <- tribble(
    ~group,       ~level1,      ~palette,
    "overall",    score_label,  "overall",
    "expression", "Expression", "expression",
    "spatial",    "Spatial",    "spatial",
    "general",    "General", "general"
  )

  pal <- function(name, i = 2:7) brewer.pal(9, name)[i]
  palettes <- list(
    overall     = colorRampPalette(c("#FFF4D6", "#FCCB7E", "#EE9E55", "#DA8248"))(6),
    expression  = colorRampPalette(c("#EAF1F9", "#4F86C6"))(6),
    spatial     = colorRampPalette(c("#FBE9E5", "#D9644F"))(6),
    transcripts = colorRampPalette(c("#F1EDE8", "#7D6853"))(5),
    general     = colorRampPalette(c("#F4F2EF", "#A39A90"))(5),
    general_lt  = colorRampPalette(c("#FBFAF9", "#EAE6E1"))(3)
  )

  legends <- list(
    list(title = if (score_label == "Overall") "Overall score" else score_label, palette = "overall", geom = "bar", labels = c("0", "", "0.5", "", "1"),
         size = c(0, 0.25, 0.5, 0.75, 1), color = palettes$overall[c(1, 2, 4, 5, 6)]),
    list(title = "Scaled metric", palette = "expression", geom = "funkyrect",
         labels = c("min", "", "", "", "max"), size = c(0.1, 0.3, 0.5, 0.75, 1),
         color = pal("Greys", c(2, 3, 5, 6, 7)), label_hjust = 0),
    list(palette = "spatial",     enabled = FALSE),
    list(palette = "transcripts", enabled = FALSE),
    list(palette = "general",     enabled = FALSE),
    list(palette = "general_lt",  enabled = FALSE)
  )

  g <- funky_heatmap(
    df_funky,
    column_info   = column_info,
    column_groups = column_groups,
    palettes      = palettes,
    legends       = legends,
    scale_column  = FALSE,
    position_args = position_arguments(
      row_height = row_height, row_space = row_space, col_space = 0.08,
      expand_xmax = 2.6, expand_ymin = if (is.null(note)) 0 else 3.2,
      expand_ymax = 0.3, col_annot_offset = 3.85, col_annot_angle = 45
    )
  )

  p <- g[[1]]
  rl <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomRect") && "column_id" %in% names(l$data)))
  p$layers[[rl]]$data <- mutate(p$layers[[rl]]$data, d = ifelse(column_id %in% "cells_per_mm2", 0.08, 0),
                                xmin = xmin + d, xmax = xmax - d, ymin = ymin + d, ymax = ymax - d)
  pos <- bind_rows(lapply(p$layers, function(l) mutate(l$data[intersect(names(l$data), c("label_value", "row_id", "xmin", "xmax", "ymin", "ymax"))], across(everything(), as.character)))) %>%
    mutate(across(c(xmin, xmax, ymin, ymax), as.numeric))
  col_x <- function(i) with(filter(pos, label_value == column_info$name[column_info$id == i])[1, ], c(xmin, xmax))
  rows <- pos %>% filter(row_id %in% df_funky$id) %>% group_by(row_id) %>% summarise(y = (min(ymin) + max(ymax)) / 2)
  z <- col_x("low_z_integrity_vs_tissue")
  na <- bind_rows(tibble(row_id = df_funky$id[df_funky$na_ficture], x = mean(col_x("ficture_f1_score"))),
                  tibble(row_id = df_funky$id[df_funky$na_z], x = mean(z))) %>% left_join(rows, by = "row_id")
  hdr <- which(sapply(p$layers, function(l) inherits(l$geom, "GeomText")))
  p$layers[[hdr]]$data <- mutate(p$layers[[hdr]]$data, size = case_when(label_value %in% column_groups$level1 ~ fs(4.2, 5), size == fs(3.3, 4) ~ size, TRUE ~ fs(size)))
  # two-line density header: centre the lines on each other, keep the start at the column
  half <- strwidth("Cells/mm²", units = "inches", cex = with(p$layers[[hdr]]$data, size[grepl("\n", label_value)][1]) * .pt / par("ps")) / 0.25 / 2
  p$layers[[hdr]]$data <- mutate(p$layers[[hdr]]$data, two = grepl("\n", label_value), hjust = ifelse(two, 0.5, hjust), vjust = ifelse(two, 0.75, vjust),
                                 x = x + two * half * cos(angle * pi / 180), y = y + two * half * sin(angle * pi / 180))
  bar <- with(filter(pos, label_value == score_label)[1, ], c(xmin, xmax))
  len <- df_funky$total_score * diff(bar)
  dims <- c("3D" = "#1f1f1f", "2D" = "#b3b3b3")
  ry <- rows$y[match(df_funky$id, rows$row_id)]
  pitch <- (row_height + row_space) / 2
  dash <- function(y, h) c(outer(y, (seq_len(4) - 2.5) * h / 2, `+`))
  bl <- dash(ry[df_funky$dim == "Baseline"], pitch)
  g[[1]] <- p +
    annotate("rect", xmin = bar[1] - 0.5, xmax = bar[1] - 0.35, ymin = ry - pitch, ymax = ry + pitch, fill = dims[df_funky$dim],
             colour = NA) +
    annotate("rect", xmin = bar[1] - 0.5, xmax = bar[1] - 0.35, ymin = bl - 0.08, ymax = bl + 0.08, fill = "#b3b3b3", colour = NA) +
    annotate("text", x = bar[1] + ifelse(len > 1.3, 0.08, len + 0.2), y = rows$y[match(df_funky$id, rows$row_id)],
             label = sprintf("%.2f", df_funky$total_score), hjust = 0, size = fs(3.85, 5.5), family = family) +
    annotate("segment", x = mean(c(z[2], col_x("pct_assigned_qced")[1])), y = min(rows$y) - 0.5, yend = max(rows$y) + 0.5,
             colour = "grey60", linewidth = 0.3) +
    annotate("text", x = na$x, y = na$y, label = "–", colour = "grey60", size = fs(3.3, 5.5), family = family) +
    annotate("text", x = min(pos$xmin, na.rm = TRUE), y = min(rows$y) - 0.9, label = if (is.null(note)) "" else note,
             hjust = 0, vjust = 1, lineheight = 1, colour = "grey30", size = fs(3.3, 5.5), family = family)
  legs <- lapply(seq_len(length(g[[3]])), function(i) {
    l <- g[[3]][[i]]
    for (j in which(sapply(l$layers, function(x) inherits(x$geom, "GeomText"))))
      l$layers[[j]]$data <- mutate(l$layers[[j]]$data, size = fs(size), y = y + ifelse(y < -3, 0.15, 0))
    for (j in which(sapply(l$layers, function(x) inherits(x$geom, "GeomRect") && "i" %in% names(x$data))))
      l$layers[[j]]$data <- arrange(l$layers[[j]]$data, is.na(i))
    if (i == 1) {
      for (j in seq_along(l$layers))
        l$layers[[j]]$data <- mutate(l$layers[[j]]$data, across(any_of(c("x", "xmin", "xmax")), ~ .x * 0.7))
      l$width <- l$width * 0.7
    }
    l
  })
  ky <- c(-1, -2, -3)
  kd <- dash(ky[3], 0.35)
  kw <- 0.4 + strwidth("Baseline", units = "inches", cex = fs(3.96) * .pt / par("ps")) / 0.25
  key <- ggplot() +
    annotate("rect", xmin = 0, xmax = 0.15, ymin = ky[1:2] - 0.35, ymax = ky[1:2] + 0.35, fill = dims) +
    annotate("rect", xmin = 0, xmax = 0.15, ymin = kd - 0.05, ymax = kd + 0.05, fill = "#b3b3b3") +
    annotate("text", x = 0.4, y = ky, label = c(names(dims), "Baseline"), hjust = 0, vjust = 0.5, size = fs(3.96), family = family) +
    annotate("blank", x = c(0, kw), y = c(-3.5, -0.5)) + theme_void()
  key$width <- 0.25 * kw
  legs <- c(legs, list(key))
  lw <- sapply(legs, `[[`, "width") / 0.25
  gaps <- c(1.2, 1.4)
  leg <- patchwork::wrap_plots(head(c(rbind(legs, list(patchwork::plot_spacer()))), -1), nrow = 1, widths = head(c(rbind(lw, c(gaps, 0))), -1)) &
    coord_equal(expand = FALSE, clip = "off") & theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_blank())
  top <- max(rows$y) + 1.1
  g <- g[[1]] + annotation_custom(patchwork::patchworkGrob(leg), xmin = min(pos$xmin, na.rm = TRUE), xmax = min(pos$xmin, na.rm = TRUE) + sum(lw) + sum(gaps),
                                  ymin = top, ymax = top + 3)
  g <- g + coord_equal(expand = FALSE, clip = "off") + theme(plot.margin = margin(2, 2, 2, 2))

  g <- g & theme(plot.background = element_rect(fill = NA, colour = NA))
  list(g = g, w = p$width, h = p$height)
  }

  name_w <- 8.3
  row_space <- 0.08
  f <- 6 / .pt / 4.2
  fit <- function() {
    b <- build(name_w)
    if (is.null(width)) return(b)
    w0 <- b$w - 0.25 * (name_w - 8.3)
    name_w <<- 8.3 * f * (w0 - 2.075) / (width - 2.075 * f)
    s <<- width / (w0 + 0.25 * (name_w - 8.3))
    build(name_w)
  }
  b <- fit()
  if (!is.null(height)) {
    row_space <- row_space + (height / s - b$h) / (0.25 * nrow(df_funky))
    b <- build(name_w)
  }
  g <- b$g; w <- b$w; h <- b$h
  showtext::showtext_opts(dpi = 600 * s)
  ggsave(paste0(out_file, ".png"), g, width = w, height = h, dpi = 600 * s, bg = "transparent", type = "cairo")
  svg <- paste0(out_file, ".svg")
  showtext::showtext_auto(FALSE)
  try(systemfonts::register_font(family, plain = path.expand(font_file)), silent = TRUE)
  ggsave(svg, g, width = w, height = h, bg = "transparent", device = svglite::svglite)
  if (s != 1) writeLines(sub("width='[0-9.]+pt' height='[0-9.]+pt'", sprintf("width='%.2fpt' height='%.2fpt'", 72 * w * s, 72 * h * s), readLines(svg)), svg)
}
