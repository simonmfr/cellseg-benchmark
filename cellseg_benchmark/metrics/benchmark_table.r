plot_benchmark_heatmap <- function(df_funky, out_file, metric_labels, family = "Arial", font_file, score_label = "Overall", row_height = 1.02, width = NULL) {
  library(funkyheatmap)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(RColorBrewer)

  sysfonts::font_add(family, path.expand(font_file))
  showtext::showtext_auto()
  showtext::showtext_opts(dpi = 600)
  update_geom_defaults("text", list(family = family))
  for (s in c("sep1", "sep2", "sep3")) df_funky[[s]] <- ""
  # width (inches): text sized for print at this width (6 pt; group headers 5, scores 5.5, densities 4.5)
  s <- 1
  fs <- function(native, pt = 6) if (is.null(width)) native else pt / .pt / s
  build <- function(name_w) {

  box <- list(width = 1.02)
  col_spec <- tribble(
    ~id,                          ~group,       ~geom,       ~palette,      ~name,          ~opts,
    "id",                         "",           "text",      NA,            "",             list(width = name_w, hjust = 0, size = fs(4.2)),
    "total_score",                "overall",    "bar",       "overall",     "",             list(width = 2.9, draw_outline = FALSE),
    "sep1",                       "",           "text",      NA,            "",             list(width = 0.3),
    "MECR_adj",                   "expression", "funkyrect", "expression",  "MECR",         box,
    "marker_f1",                  "expression", "funkyrect", "expression",  "Marker F1",    box,
    "silhouette_score",           "expression", "funkyrect", "expression",  "Silhouette",   box,
    "sep2",                       "",           "text",      NA,            "",             list(width = 0.3),
    "ficture_f1_score",           "spatial",    "funkyrect", "spatial",     "Ficture F1",   box,
    "low_z_integrity_vs_tissue",  "spatial",    "funkyrect", "spatial",     "Z integrity",  box,
    "sep3",                       "",           "text",      NA,            "",             list(width = 0.6),
    "pct_assigned_qced",          "general",    "funkyrect", "transcripts", "% Assigned",   box,
    "volume_final",               "general",    "funkyrect", "general",     "Volume",       box,
    "elongation",                 "general",    "funkyrect", "general",     "Elongation",   box,
    "intensities_DAPI",           "general",    "funkyrect", "general",     "DAPI",         box,
    "cells_per_mm2",              "general",    "rect",      "general_lt",  "Cells/mm²",    list(width = 1.1),
    "density_label",              "general",    "text",      NA,            "",             list(overlay = TRUE, size = fs(3.3, 4.5), hjust = 0.5),
  )

  column_info <- col_spec %>%
    filter(id %in% names(df_funky)) %>%
    left_join(metric_labels, by = "id", suffix = c("", "_label")) %>%
    mutate(name = coalesce(name_label, name), name = ifelse(id == "cells_per_mm2", "Cells/mm² (×10³)", name),
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
    expression  = pal("Purples", 2:6),
    spatial     = colorRampPalette(c("#FCE6DF", "#D67A6B"))(6),
    transcripts = pal("Blues", 2:6),
    general     = colorRampPalette(c("#ECEFF4", "#7B8BA3"))(5),
    general_lt  = colorRampPalette(c("#F8F9FB", "#DDE2EA"))(3)
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
      row_height = row_height, row_space = 0.08, col_space = 0.08,
      expand_xmax = 4.5, expand_ymin = 0,
      expand_ymax = 0.3, col_annot_offset = 4, col_annot_angle = 45
    )
  )

  p <- g[[1]]
  for (j in which(sapply(p$layers, function(l) inherits(l$geom, "GeomRoundedRect"))))
    p$layers[[j]]$data <- mutate(p$layers[[j]]$data, cx = (xmin + xmax) / 2, cy = (ymin + ymax) / 2,
                                 xmin = cx + 0.92 * (xmin - cx), xmax = cx + 0.92 * (xmax - cx), ymin = cy + 0.92 * (ymin - cy), ymax = cy + 0.92 * (ymax - cy))
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
  p$layers[[hdr]]$data <- mutate(p$layers[[hdr]]$data, size = case_when(label_value %in% column_groups$level1 ~ fs(4.2, 5), size == fs(3.3, 4.5) ~ size, TRUE ~ fs(size)))
  bar <- with(filter(pos, label_value == score_label)[1, ], c(xmin, xmax))
  len <- df_funky$total_score * diff(bar)
  g[[1]] <- p +
    annotate("text", x = bar[1] + ifelse(len > 1.3, 0.08, len + 0.2), y = rows$y[match(df_funky$id, rows$row_id)],
             label = sprintf("%.2f", df_funky$total_score), hjust = 0, size = fs(3.85, 5.5), family = family) +
    annotate("segment", x = mean(c(z[2], col_x("pct_assigned_qced")[1])), y = min(rows$y) - 0.5, yend = max(rows$y) + 0.5,
             colour = "grey60", linewidth = 0.3) +
    annotate("text", x = na$x, y = na$y, label = "–", colour = "grey60", size = fs(3.3, 5.5), family = family)
  legs <- lapply(seq_len(length(g[[3]])), function(i) {
    l <- g[[3]][[i]]
    for (j in which(sapply(l$layers, function(x) inherits(x$geom, "GeomText"))))
      l$layers[[j]]$data <- mutate(l$layers[[j]]$data, size = fs(size))
    for (j in which(sapply(l$layers, function(x) inherits(x$geom, "GeomRect") && "i" %in% names(x$data))))
      l$layers[[j]]$data <- arrange(l$layers[[j]]$data, is.na(i))
    l
  })
  lw <- sapply(legs, `[[`, "width") / 0.25
  leg <- patchwork::wrap_plots(c(legs[1], list(patchwork::plot_spacer()), legs[-1]), nrow = 1, widths = c(lw[1], 1.5, lw[-1])) &
    coord_equal(expand = FALSE, clip = "off") & theme(plot.margin = margin(0, 0, 0, 0), plot.background = element_blank())
  top <- max(rows$y) + 1.1
  g <- g[[1]] + annotation_custom(patchwork::patchworkGrob(leg), xmin = min(pos$xmin, na.rm = TRUE), xmax = min(pos$xmin, na.rm = TRUE) + sum(lw) + 1.5,
                                  ymin = top, ymax = top + 3)
  g <- g + coord_equal(expand = FALSE, clip = "off") + theme(plot.margin = margin(2, 2, 2, 2))

  g <- g & theme(plot.background = element_rect(fill = NA, colour = NA))
  list(g = g, w = p$width, h = p$height)
  }

  b <- build(8.3)
  if (!is.null(width)) {
    f <- 6 / .pt / 4.2
    name_w <- 8.3 * f * (b$w - 2.075) / (width - 2.075 * f)
    s <- width / (b$w + 0.25 * (name_w - 8.3))
    b <- build(name_w)
  }
  g <- b$g; w <- b$w; h <- b$h
  showtext::showtext_opts(dpi = 600 * s)
  ggsave(paste0(out_file, ".png"), g, width = w, height = h, dpi = 600 * s, bg = "transparent", type = "cairo")
  svg <- paste0(out_file, ".svg")
  ggsave(svg, g, width = w, height = h, bg = "transparent")
  if (s != 1) writeLines(sub("width='[0-9.]+pt' height='[0-9.]+pt'", sprintf("width='%.2fpt' height='%.2fpt'", 72 * w * s, 72 * h * s), readLines(svg)), svg)
}
