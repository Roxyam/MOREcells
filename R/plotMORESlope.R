# =============================================================================
# plotSlopeDumbbell(): the slope of each side of a contrast, for many pairs
#
# Colours and theme are the scMore ones (plotPalettes.R, plotTheme.R). It does not compute
# anything: it draws columns that contrastPatterns() already returns.
# =============================================================================

#' Slope of each side of a contrast, for many regulator-target pairs
#'
#' Overview of many pairs of a contrast. Each pair is one position on the x axis, in panels
#' (one per contrast and level of \code{by}); the vertical axis is the slope. It reads the
#' output of \code{contrastPatterns()}, so the slope of each side and whether it is detected
#' are already there.
#'
#' @section How to read it:
#' \itemize{
#'   \item A \strong{diamond} is the slope of side 1 and a \strong{square} the slope of side 2,
#'     joined by a grey segment. The pairs are ordered by the contrast estimate.
#'   \item The \strong{dash} is the contrast estimate (side 1 minus side 2): dark if the
#'     contrast is significant, light grey if not (\code{fadeNonSignificant}).
#'   \item A \strong{filled} point means that slope is detected on its own (its confidence
#'     interval excludes 0); a \strong{hollow} one means it is not detected, which is not the
#'     same as equal to zero.
#' }
#' So a segment that crosses the dashed zero line is a sign switch; a filled and a hollow
#' point is "detected in one level only"; two filled points on the same side of zero are
#' "in both".
#'
#' @section Style:
#' The two sides take the colours of group A and group B of the contrast networks
#' (\code{\link{scMoreColors}("owner")}), so a group has the same colour in the dumbbell and in
#' the network, and \code{scMoreSetColors("owner", ...)} changes both. The theme is
#' \code{\link{scMoreTheme}()}.
#'
#' @param x Output of \code{contrastPatterns()} for a contrast between two levels (not an
#'   interaction).
#' @param top \code{NULL} (all the pairs) or the number of pairs per panel; the most
#'   significant ones, by \code{p.adj}.
#' @param onlySignificant Keep only the significant contrasts.
#' @param labels \code{"auto"} (labels only if a panel has at most \code{maxLabels} pairs),
#'   \code{"yes"} or \code{"no"}.
#' @param maxLabels Largest number of pairs per panel that still gets labels with
#'   \code{labels = "auto"}.
#' @param labelBy What the x labels say: \code{"pair"}, \code{"regulator"} or
#'   \code{"target"}.
#' @param order Order of the pairs inside a panel: \code{"estimate"} (the contrast estimate)
#'   or \code{"p.adj"}.
#' @param decreasing Reverse the order.
#' @param fadeNonSignificant Draw the dash of the non-significant contrasts in light grey.
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' x <- contrastPatterns(contrastRegulation(fit, compare = "Strain", by = "AnnotClass"))
#' plotSlopeDumbbell(x, top = 18)
#' plotSlopeDumbbell(x, onlySignificant = TRUE, labels = "no")
#' plotSlopeDumbbell(x, top = 15, labelBy = "target", order = "p.adj")
#' scMoreSave(plotSlopeDumbbell(x, top = 18), "dumbbell.png", size = "double")
#' }
#' @family scMore plots
#' @importFrom rlang .data
#' @importFrom utils head
#' @export
plotSlopeDumbbell <- function(x, top = NULL, onlySignificant = FALSE,
                              labels = c("auto", "yes", "no"), maxLabels = 30,
                              labelBy = c("pair", "regulator", "target"),
                              order = c("estimate", "p.adj"), decreasing = FALSE,
                              fadeNonSignificant = TRUE) {

  labels <- match.arg(labels); labelBy <- match.arg(labelBy); order <- match.arg(order)
  if (!exists("scMoreColors") || !exists("scMoreTheme")) {
    cli::cli_abort("{.fn plotSlopeDumbbell} needs {.file plotPalettes.R} and {.file plotTheme.R} (scMoreColors and scMoreTheme).")
  }
  a <- attr(x, "args")
  if (is.null(a) || !"pattern" %in% names(x)) cli::cli_abort("{.arg x} must be the output of {.fn contrastPatterns}.")
  if (isTRUE(a$interaction)) cli::cli_abort("For interaction contrasts use {.fn plotSlopeInteraction}.")

  # ---- Data: one row per pair, in panels (one per contrast and level of `by`) ----
  d <- x[!is.na(x$estimate) & !is.na(x$slope_1) & !is.na(x$slope_2), ]
  if (onlySignificant) d <- d[d$significant %in% TRUE, ]
  if (nrow(d) == 0) cli::cli_abort("There are no pairs to plot.")
  d$panel <- do.call(paste, c(lapply(d[c("contrast", a$by)], as.character), sep = " | "))
  if (!is.null(top)) d <- do.call(rbind, lapply(split(d, d$panel), function(z) utils::head(z[base::order(z$p.adj), ], top)))
  d$pair  <- switch(labelBy, pair = paste(d$target, d$regulator, sep = " - "), regulator = d$regulator, target = d$target)
  key     <- if (order == "estimate") d$estimate else d$p.adj
  d       <- d[base::order(d$panel, if (decreasing) -key else key), ]
  d$xid   <- factor(paste(d$pair, d$panel, sep = "\r"), levels = unique(paste(d$pair, d$panel, sep = "\r")))   # order inside each panel
  n_panel <- max(table(d$panel)); multi <- length(unique(d$panel)) > 1

  # ---- Colours: the two sides are group A and group B of the contrast networks ----
  owner <- scMoreColors("owner")
  col   <- c("1" = owner[["A"]], "2" = owner[["B"]])
  ink   <- "#2B2B2B"; ink_light <- "#BDBDBD"          # the dash of the contrast: dark, or light grey if not significant
  one   <- length(unique(d$contrast)) == 1
  side_names <- if (one) c(d$side_1[1], d$side_2[1]) else c("First side of the contrast", "Second side of the contrast")

  pts <- rbind(
    data.frame(xid = d$xid, panel = d$panel, side = "1", slope = d$slope_1, filled = d$state_1 %in% c("+", "-")),
    data.frame(xid = d$xid, panel = d$panel, side = "2", slope = d$slope_2, filled = d$state_2 %in% c("+", "-")))
  pts$fillc <- ifelse(pts$filled, col[pts$side], "white")

  # Sizes shrink as the number of pairs grows
  pt <- if (n_panel > 200) 1.0 else if (n_panel > 80) 1.6 else if (n_panel > 30) 2.3 else 3
  lw <- if (n_panel > 200) 0.35 else if (n_panel > 80) 0.5 else 0.8
  sig_c <- d$significant %in% TRUE

  show_labels <- labels == "yes" || (labels == "auto" && n_panel <= maxLabels)
  caption <- paste(strwrap(paste0("Dash: contrast estimate (side 1 - side 2)", if (fadeNonSignificant) ", light grey if not significant", ". ",
                                  "Filled point: slope detected on its own (interval excludes 0); hollow: not detected."), width = 80), collapse = "\n")
  dash <- function(data, colour) {
    ggplot2::geom_errorbar(data = data, ggplot2::aes(x = .data$xid, ymin = .data$estimate, ymax = .data$estimate),
                           colour = colour, width = 0.75, linewidth = lw + 0.3)
  }

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
    ggplot2::geom_segment(data = d, ggplot2::aes(x = .data$xid, xend = .data$xid, y = .data$slope_1, yend = .data$slope_2),
                          colour = "grey60", linewidth = lw) +
    dash(d[sig_c | !fadeNonSignificant, ], ink) +
    {if (fadeNonSignificant && any(!sig_c)) dash(d[!sig_c, ], ink_light)} +
    ggplot2::geom_point(data = pts, ggplot2::aes(x = .data$xid, y = .data$slope, colour = .data$side, shape = .data$side, fill = .data$fillc),
                        size = pt, stroke = 0.9) +
    ggplot2::scale_fill_identity() +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(fill = unname(col)))) +
    ggplot2::scale_colour_manual(values = col, labels = side_names, name = if (one) NULL else "Side") +
    ggplot2::scale_shape_manual(values = c("1" = 23, "2" = 22), labels = side_names, name = if (one) NULL else "Side") +
    ggplot2::scale_x_discrete(labels = function(v) sub("\r.*$", "", v),
                              guide = ggplot2::guide_axis(n.dodge = if (show_labels && n_panel <= 8) 2 else 1)) +
    ggplot2::labs(title = "Slope of each side of the contrast",
                  subtitle = sprintf("%d pairs per panel at most, ordered by %s", n_panel, if (order == "estimate") "the contrast estimate" else "p.adj"),
                  x = NULL, y = "Slope", caption = caption) +
    scMoreTheme() +
    ggplot2::theme(plot.caption = ggplot2::element_text(size = 8, colour = "grey40", hjust = 0),
                   panel.grid.major.y = ggplot2::element_line(colour = "#E5E5E5", linewidth = 0.3),
                   panel.grid.major.x = ggplot2::element_blank())
  if (!show_labels) {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank()) +
      ggplot2::labs(x = sprintf("%d pairs ordered by %s", n_panel, if (order == "estimate") "contrast estimate" else "p.adj"))
  } else if (n_panel > 8) {                      # many labels: vertical, so that they do not overlap
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = if (n_panel > 20) 7 else 8))
  } else {
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8))
  }
  if (multi) p <- p + ggplot2::facet_wrap(ggplot2::vars(.data$panel), scales = "free_x")
  p
}
