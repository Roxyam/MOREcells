# =============================================================================
# MOREcells plot style: themes, session options and figure export
#
# Look: minimalist, transparent backgrounds (panel, plot and legend have
# fill = NA, so figures sit on any slide or page colour), thin dark axes or an
# optional marked panel border, no grid by default, legend at the bottom.
#
# Every MOREcells plot function ends with `MOREcellsTheme()`, so the user can:
#   * keep the default look           (theme_MOREcells),
#   * add any other theme afterwards  (plot + ggplot2::theme_bw()),
#   * change the theme of every MOREcells plot at once with MOREcellsSetTheme(),
#   * use the package themes in their own ggplots (+ theme_MOREcells()).
# =============================================================================

.MOREcellsInk   <- "#2B2B2B"
.MOREcellsMid   <- "#7F7F7F"
.MOREcellsFaint <- "#E5E5E5"

# Transparent backgrounds, text, legend, title and strip formatting shared by
# every MOREcells theme
.MOREcellsCommon <- function(base_size, legend, legend_title = TRUE) {
  rel   <- ggplot2::rel
  blank <- ggplot2::element_rect(fill = NA, colour = NA)
  ggplot2::theme(
    text               = ggplot2::element_text(colour = .MOREcellsInk),
    plot.title         = ggplot2::element_text(face = "bold", size = rel(1.15), hjust = 0,
                                               margin = ggplot2::margin(b = base_size * 0.5)),
    plot.subtitle      = ggplot2::element_text(colour = .MOREcellsMid, size = rel(0.95), hjust = 0,
                                               margin = ggplot2::margin(b = base_size * 0.5)),
    plot.caption       = ggplot2::element_text(colour = .MOREcellsMid, size = rel(0.8), hjust = 1,
                                               margin = ggplot2::margin(t = base_size * 0.5)),
    plot.title.position   = "plot",
    plot.caption.position = "plot",
    plot.margin        = ggplot2::margin(base_size * 0.7, base_size * 0.9,
                                         base_size * 0.7, base_size * 0.7),
    panel.background   = blank,
    plot.background    = blank,
    legend.background  = blank,
    legend.box.background = blank,
    legend.key         = ggplot2::element_blank(),
    legend.position    = legend,
    legend.title       = if (legend_title) ggplot2::element_text(face = "bold", size = rel(0.9))
                         else ggplot2::element_blank(),
    legend.text        = ggplot2::element_text(size = rel(0.85)),
    strip.background   = blank,
    strip.text         = ggplot2::element_text(face = "bold", size = rel(0.9), colour = .MOREcellsInk,
                                               margin = ggplot2::margin(4, 4, 4, 4)),
    panel.grid.minor   = ggplot2::element_blank(),
    panel.grid.major   = ggplot2::element_blank(),
    panel.spacing      = ggplot2::unit(base_size * 0.9, "pt")
  )
}

# Axis treatment: L-shaped axis lines, a marked panel border, or nothing
.MOREcellsAxes <- function(axes) {
  tick <- ggplot2::element_line(colour = .MOREcellsInk, linewidth = 0.4)
  switch(axes,
    lines = ggplot2::theme(
      axis.line   = ggplot2::element_line(colour = .MOREcellsInk, linewidth = 0.4, lineend = "square"),
      axis.ticks  = tick,
      axis.ticks.length = ggplot2::unit(3, "pt"),
      panel.border = ggplot2::element_blank()
    ),
    box = ggplot2::theme(
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = tick,
      axis.ticks.length = ggplot2::unit(3, "pt"),
      panel.border = ggplot2::element_rect(colour = .MOREcellsInk, fill = NA, linewidth = 0.7),
      strip.background = ggplot2::element_rect(fill = NA, colour = .MOREcellsInk, linewidth = 0.5)
    ),
    none = ggplot2::theme(
      axis.line   = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    )
  )
}

#' MOREcells ggplot2 themes
#'
#' The themes shared by every MOREcells plot. All of them are minimalist, with
#' transparent panel, plot and legend backgrounds, no minor grid and the legend
#' at the bottom.
#'
#' \itemize{
#'   \item \code{theme_MOREcells()}: general purpose (scatter, bars, volcano, forest).
#'     Thin dark axes with square corners, or a marked panel border with
#'     \code{axes = "box"}.
#'   \item \code{theme_MOREcells_min()}: minimal. No axis lines or ticks, no axis
#'     titles and no legend title by default. Use \code{axes = "box"} for a
#'     marked border and \code{axis_text = FALSE} to hide the axis labels.
#'   \item \code{theme_MOREcells_umap()}: \code{theme_MOREcells_min()} without axis
#'     labels, for embeddings (UMAP, PCA).
#'   \item \code{theme_MOREcells_tile()}: heatmaps and tile plots. No axes or grid,
#'     rotated x labels.
#'   \item \code{theme_MOREcells_void()}: networks, diagrams and the palette gallery.
#'     No axes at all.
#' }
#'
#' Text sizes are in points. At the "single column" size of
#' \code{\link{MOREcellsFigureSize}} (85 mm) \code{base_size = 12} looks large and
#' \code{base_size = 8} is more usual for print.
#'
#' @param base_size Base font size (pt).
#' @param base_family Font family.
#' @param grid Major grid lines: \code{"none"}, \code{"y"}, \code{"x"} or \code{"xy"}.
#' @param axes \code{"lines"} (L-shaped axis lines), \code{"box"} (marked panel
#'   border) or \code{"none"}.
#' @param legend Legend position.
#' @param axis_text,axis_title Show the axis labels / titles.
#' @param legend_title Show the legend title.
#' @param x_angle Rotation of the x labels in \code{theme_MOREcells_tile()}.
#' @param ... Further arguments to \code{\link[ggplot2]{theme}}, applied last.
#' @return A ggplot2 theme.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) +
#'   geom_point() +
#'   theme_MOREcells()
#' @family MOREcells plot style
#' @export
theme_MOREcells <- function(base_size = 12, base_family = "",
                         grid = c("none", "y", "x", "xy"),
                         axes = c("lines", "box", "none"),
                         legend = c("bottom", "right", "top", "left", "none"),
                         ...) {
  grid   <- match.arg(grid)
  axes   <- match.arg(axes)
  legend <- match.arg(legend)

  th <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    .MOREcellsCommon(base_size, legend) +
    .MOREcellsAxes(axes) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = ggplot2::rel(0.95)),
      axis.text  = ggplot2::element_text(size = ggplot2::rel(0.85), colour = .MOREcellsInk)
    )

  gl <- ggplot2::element_line(colour = .MOREcellsFaint, linewidth = 0.3)
  if (grid %in% c("y", "xy")) th <- th + ggplot2::theme(panel.grid.major.y = gl)
  if (grid %in% c("x", "xy")) th <- th + ggplot2::theme(panel.grid.major.x = gl)
  th + ggplot2::theme(...)
}

#' @rdname theme_MOREcells
#' @export
theme_MOREcells_min <- function(base_size = 12, base_family = "",
                             axes = c("none", "box", "lines"),
                             axis_text = TRUE, axis_title = FALSE,
                             legend = c("bottom", "right", "top", "left", "none"),
                             legend_title = FALSE, ...) {
  axes   <- match.arg(axes)
  legend <- match.arg(legend)
  th <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    .MOREcellsCommon(base_size, legend, legend_title) +
    .MOREcellsAxes(axes) +
    ggplot2::theme(axis.text = ggplot2::element_text(size = ggplot2::rel(0.85), colour = .MOREcellsInk))
  if (!axis_text)  th <- th + ggplot2::theme(axis.text  = ggplot2::element_blank(),
                                             axis.ticks = ggplot2::element_blank())
  if (!axis_title) th <- th + ggplot2::theme(axis.title = ggplot2::element_blank())
  th + ggplot2::theme(...)
}

#' @rdname theme_MOREcells
#' @export
theme_MOREcells_umap <- function(base_size = 12, base_family = "",
                              axes = c("none", "box", "lines"),
                              legend = c("bottom", "right", "top", "left", "none"),
                              legend_title = FALSE, ...) {
  theme_MOREcells_min(base_size, base_family, axes = match.arg(axes), axis_text = FALSE,
                   axis_title = FALSE, legend = match.arg(legend),
                   legend_title = legend_title, ...)
}

#' @rdname theme_MOREcells
#' @export
theme_MOREcells_tile <- function(base_size = 12, base_family = "",
                              legend = c("right", "bottom", "top", "left", "none"),
                              x_angle = 45, ...) {
  theme_MOREcells_min(base_size, base_family, axes = "none", axis_text = TRUE,
                   axis_title = FALSE, legend = match.arg(legend), legend_title = TRUE) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_angle, hjust = if (x_angle == 0) 0.5 else 1,
                                          vjust = 1)
    ) +
    ggplot2::theme(...)
}

#' @rdname theme_MOREcells
#' @export
theme_MOREcells_void <- function(base_size = 12, base_family = "",
                              legend = c("bottom", "right", "top", "left", "none"),
                              ...) {
  legend <- match.arg(legend)
  ggplot2::theme_void(base_size = base_size, base_family = base_family) +
    .MOREcellsCommon(base_size, legend) +
    ggplot2::theme(...)
}

# ---- session options ---------------------------------------------------------

#' Theme used by all MOREcells plots
#'
#' \code{MOREcellsTheme()} returns the theme that MOREcells plot functions apply. By
#' default it is one of the package themes (\code{variant} picks which).
#' \code{MOREcellsSetTheme()} replaces it for the whole session, with either a
#' theme object or a function returning one.
#'
#' A theme set with \code{MOREcellsSetTheme()} replaces the package theme of
#' \emph{every} MOREcells plot, whatever its \code{variant}. To change a single
#' plot, add a theme to it (\code{plot + ggplot2::theme_bw()}).
#'
#' @param theme A ggplot2 theme, a function returning one (e.g.
#'   \code{function() theme_MOREcells(base_size = 9)}), or \code{NULL} to restore
#'   the default.
#' @param variant Which package theme to use when no user theme is set:
#'   \code{"default"}, \code{"min"}, \code{"umap"}, \code{"tile"} or \code{"void"}.
#' @param ... Passed to the package theme (or to the user function).
#' @return \code{MOREcellsTheme()}: a ggplot2 theme. \code{MOREcellsSetTheme()}: the
#'   previous setting, invisibly.
#' @examples
#' MOREcellsSetTheme(function() theme_MOREcells(base_size = 9, grid = "y"))
#' MOREcellsSetTheme(NULL)
#' @family MOREcells plot style
#' @export
MOREcellsTheme <- function(variant = c("default", "min", "umap", "tile", "void"), ...) {
  variant <- match.arg(variant)
  user <- getOption("MOREcells.theme")
  if (!is.null(user)) {
    if (is.function(user)) return(if (length(list(...))) user(...) else user())
    return(user)
  }
  switch(variant,
         default = theme_MOREcells(...),
         min     = theme_MOREcells_min(...),
         umap    = theme_MOREcells_umap(...),
         tile    = theme_MOREcells_tile(...),
         void    = theme_MOREcells_void(...))
}

#' @rdname MOREcellsTheme
#' @export
MOREcellsSetTheme <- function(theme = NULL) {
  if (!is.null(theme) && !is.function(theme) && !inherits(theme, "theme")) {
    cli::cli_abort("{.arg theme} must be a ggplot2 theme, a function returning one, or {.code NULL}.")
  }
  old <- options(MOREcells.theme = theme)
  invisible(old$MOREcells.theme)
}

#' Apply the MOREcells theme globally
#'
#' Convenience for interactive work: sets \code{\link[ggplot2]{theme_set}} to
#' \code{\link{theme_MOREcells}} so that your own ggplots also use it.
#'
#' @param ... Passed to \code{\link{theme_MOREcells}}.
#' @return The previous theme, invisibly.
#' @family MOREcells plot style
#' @export
MOREcellsUseTheme <- function(...) {
  invisible(ggplot2::theme_set(theme_MOREcells(...)))
}

#' Reset the MOREcells style
#'
#' Undoes every session change: the colour overrides (\code{MOREcellsSetColors()}) and the
#' theme (\code{MOREcellsSetTheme()}).
#'
#' @return \code{NULL}, invisibly.
#' @family MOREcells plot style
#' @export
MOREcellsResetStyle <- function() {
  .MOREcellsEnv$colors <- list()
  options(MOREcells.theme = NULL)
  invisible(NULL)
}

# ---- figure size and export --------------------------------------------------

#' Standard figure sizes
#'
#' Width and height in millimetres of the figure presets used by
#' \code{\link{MOREcellsSave}}: \code{"single"} (one journal column), \code{"onehalf"},
#' \code{"double"} (full width), \code{"square"}, \code{"wide"} and \code{"tall"}.
#'
#' @param size Preset name. \code{NULL} returns all of them as a matrix.
#' @return A named numeric vector \code{c(width, height)} in mm (or a matrix).
#' @family MOREcells plot style
#' @export
MOREcellsFigureSize <- function(size = NULL) {
  presets <- rbind(
    single  = c(width = 85,  height = 70),
    onehalf = c(width = 120, height = 85),
    double  = c(width = 180, height = 100),
    square  = c(width = 100, height = 100),
    wide    = c(width = 180, height = 70),
    tall    = c(width = 85,  height = 120)
  )
  if (is.null(size)) return(presets)
  size <- match.arg(size, rownames(presets))
  presets[size, ]
}

#' Save a plot with a standard size and a transparent background
#'
#' Thin wrapper around \code{\link[ggplot2]{ggsave}} that takes a size preset
#' (see \code{\link{MOREcellsFigureSize}}) instead of a width and a height, and
#' saves with a transparent background by default (formats without
#' transparency, such as jpeg, get a white one).
#'
#' @param plot A ggplot object.
#' @param filename File name; the extension picks the format (png, pdf, svg...).
#' @param size Size preset. Ignored if \code{width} and \code{height} are given.
#' @param width,height Custom size, in \code{units}.
#' @param units Units of \code{width} and \code{height}.
#' @param dpi Resolution for raster formats.
#' @param transparent Save with a transparent background.
#' @param ... Passed to \code{\link[ggplot2]{ggsave}}.
#' @return The file name, invisibly.
#' @family MOREcells plot style
#' @export
MOREcellsSave <- function(plot, filename, size = "double", width = NULL, height = NULL,
                       units = "mm", dpi = 300, transparent = TRUE, ...) {
  if (is.null(width) && is.null(height)) {
    d      <- MOREcellsFigureSize(size)
    width  <- d[["width"]]
    height <- d[["height"]]
    units  <- "mm"
  }
  ext <- tolower(tools::file_ext(filename))
  bg  <- if (transparent && !ext %in% c("jpg", "jpeg", "tif", "tiff", "bmp")) "transparent" else "white"
  ggplot2::ggsave(filename = filename, plot = plot, width = width, height = height,
                  units = units, dpi = dpi, bg = bg, ...)
  invisible(filename)
}
