#' Profile plot of a regulator and its target across cells
#'
#' Single-cell version of the \code{plotMORE()} profile plot. Cells are on the
#' x axis, ordered within each level of \code{by} (blocks separated by dashed
#' lines). The regulator is drawn on the left axis and the target on the right
#' axis. Each donor gets a thin smoothed curve and each block a thick smoothed
#' curve, so the slope of the relationship (same or opposite direction) can be
#' read for every condition.
#'
#' The \code{MOREcellsFit} object only stores slopes, not the data, so the original
#' matrices and metadata must be supplied again. The fit is used to resolve
#' names, to recover the cells, the design (\code{fixedEffects},
#' \code{randomEffects}) and the layer of the regulator, and to annotate the
#' slope estimated by the model in each block.
#'
#' @param fit An object of class \code{MOREcellsFit} returned by \code{fitMore()}.
#' @param targetData Target matrix (features x cells) used in \code{fitMore()}.
#' @param regulatoryData Named list of regulatory matrices (regulators x cells).
#' @param metadata Cell metadata (row names = cell IDs).
#' @param target Name of the target feature.
#' @param regulator Name of the regulator.
#' @param layer Regulatory layer of the regulator. If \code{NULL} it is
#'   detected (it must be unambiguous).
#' @param by Metadata column(s) that define the blocks. With one variable the
#'   blocks are drawn side by side in one panel, separated by dashed lines. With
#'   two variables the plot is a grid of panels: the first goes in the rows and
#'   the second in the columns (or name them, \code{c(rows = "cell_type", cols =
#'   "condition")}), one block per panel. Defaults to the first fixed effect of
#'   the fit. For slopes to be annotated in every panel, the model must include
#'   both variables in its design.
#' @param at Named list that filters cells, as in \code{testRegulation()}, for
#'   example \code{list(cell_type = "Tcell")}. Several variables and several
#'   levels per variable are allowed. If a \code{by} variable is one of the
#'   names, its levels also set the order of the blocks. Defaults to no filter.
#' @param donor Metadata column that identifies the donors. Defaults to the
#'   first random effect of the fit.
#' @param orderBy How cells are ordered within each block. \code{"regulator"}
#'   (default): the regulator curve is monotone and the target curve shows how
#'   the target follows it. \code{"target"}: the other way round, the target
#'   curve is monotone and the regulator curve shows how the regulator changes
#'   with it. It can also be the name of a numeric metadata column such as a
#'   pseudotime.
#' @param use \code{"raw"} plots the regulator as supplied. \code{"residuals"}
#'   uses the sample-residualized regulator stored by
#'   \code{fitMore(keepResiduals = TRUE)}, which is what the slopes of a
#'   \code{twoSteps} fit refer to.
#' @param minCells Minimum number of cells of a donor in a block to draw its
#'   curve.
#' @param k Basis dimension of the GAM smoother.
#' @param showPoints Logical. Draw the cells as faint points.
#' @param collapseZeros Logical. Display only: the cells with value 0 in the
#'   ordering variable (typically the regulator) are drawn at a single x
#'   position, with one summary point per variable (mean of those cells), and
#'   the curves start from there. The model and the annotated slopes are not
#'   affected. It needs zeros in the ordering variable, so it has no effect
#'   with \code{use = "residuals"}.
#' @param annotateSlopes Logical. Write the model slope and p-value of each
#'   block (needs \code{emmeans}).
#' @param removeNA Logical. A block has no valid slope when the regulator or
#'   the target is constant in it, or when the model slope or its p-value is
#'   \code{NA}/\code{NaN}. If \code{FALSE} (default) those blocks are kept, a
#'   warning is printed and the cause is written in the block label. If
#'   \code{TRUE} they are removed from the plot.
#' @param cols Named colours for \code{regulator} and \code{target}.
#' @param title Plot title. Defaults to \code{target ~ regulator (layer)}.
#' @param seed Seed used to break ties in the cell order.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plotMORERegulation(fit, targetData, regulatoryData, metadata,
#'                 target = "gene_1", regulator = "TF1",
#'                 by = "condition",
#'                 at = list(cell_type = "Tcell", condition = c("ctrl", "treat")),
#'                 donor = "sample", orderBy = "target", removeNA = TRUE)
#'
#' # Grid: cell types in rows, conditions in columns
#' plotMORERegulation(fit, targetData, regulatoryData, metadata,
#'                 target = "gene_1", regulator = "TF1",
#'                 by = c("cell_type", "condition"), donor = "sample")
#' }
#'
#' @export
plotMORERegulation <- function(fit, targetData, regulatoryData, metadata,
                               target, regulator, layer = NULL,
                               by = NULL, at = NULL, donor = NULL,
                               orderBy = "regulator",
                               use = c("raw", "residuals"),
                               minCells = 30, k = 8,
                               showPoints = FALSE, collapseZeros = FALSE, pointAlpha = 0.3,
                               donorAlpha = 0.45, donorLinewidth = 0.3,
                               globalLinewidth = 1.4,
                               annotateSlopes = TRUE,
                               removeNA = FALSE,
                               cols = c(regulator = "#6A1B9A", target = "#D81B60"),
                               title = NULL, seed = 1808) {

  use <- match.arg(use)
  if (!is.null(at)) {
    at <- as.list(at)
    if (is.null(names(at)) || any(!nzchar(names(at)))) {
      cli::cli_abort(c(
        "{.arg at} must be a named list.",
        "i" = "For example {.code at = list(cell_type = \"Tcell\")}."
      ))
    }
  }
  if (!inherits(fit, "MOREcellsFit")) {
    cli::cli_abort("{.arg fit} must be an object of class {.cls MOREcellsFit}.")
  }
  gs   <- fit$GlobalSummary
  args <- gs$args

  # ~~~~~~~~~~ 1. Names, layer and design ~~~~~~~~~~ #
  # Accept original or sanitized names (the fit stores sanitized ones)
  resolveName <- function(nm, pool, what) {
    if (nm %in% pool) return(nm)
    hit <- which(make.names(pool) == make.names(nm))
    if (length(hit) == 1) return(pool[hit])
    cli::cli_abort("{what} {.val {nm}} was not found in the data.")
  }

  tgtRow <- resolveName(target, rownames(targetData), "Target")
  tgtKey <- make.names(tgtRow)
  if (!tgtKey %in% names(fit$Models)) {
    cli::cli_abort("Target {.val {target}} has no model in {.arg fit}.")
  }

  if (is.null(layer)) {
    layer <- names(regulatoryData)[vapply(regulatoryData, function(m) {
      any(make.names(rownames(m)) == make.names(regulator))
    }, logical(1))]
    if (length(layer) != 1) {
      cli::cli_abort("Could not detect a unique layer for {.val {regulator}}; set {.arg layer}.")
    }
  }
  regRow <- resolveName(regulator, rownames(regulatoryData[[layer]]), "Regulator")
  regKey <- make.names(regRow)
  if (!regKey %in% names(fit$Models[[tgtKey]]$fit)) {
    cli::cli_abort("Regulator {.val {regulator}} is not in the model of {.val {target}}.")
  }

  if (is.null(by))    by    <- args$fixedEffects[1]
  if (length(by) > 2) cli::cli_abort("{.arg by} must have at most two variables.")
  facet <- length(by) == 2
  if (facet) {
    by <- if (all(c("rows", "cols") %in% names(by))) c(by[["rows"]], by[["cols"]]) else unname(by)
    if (by[1] == by[2]) cli::cli_abort("The two variables of {.arg by} must be different.")
  }
  if (is.null(donor)) donor <- args$randomEffects[1]
  if (is.null(donor)) cli::cli_abort("{.arg donor} must be given (the fit has no random effects).")
  for (v in c(by, donor)) {
    if (!v %in% colnames(metadata)) cli::cli_abort("{.val {v}} is not a column of {.arg metadata}.")
  }

  # ~~~~~~~~~~ 2. Cells and values ~~~~~~~~~~ #
  cells <- if (is.character(gs$cells)) gs$cells else rownames(metadata)
  cells <- Reduce(intersect, list(cells, colnames(targetData),
                                  colnames(regulatoryData[[layer]]), rownames(metadata)))
  if (!length(cells)) cli::cli_abort("No cells are shared by the data and {.arg metadata}.")

  getRow <- function(m, r, cells) as.numeric(as.matrix(m[r, cells, drop = FALSE]))
  tgtVal <- getRow(targetData, tgtRow, cells)

  if (use == "raw") {
    regVal <- getRow(regulatoryData[[layer]], regRow, cells)
  } else {
    res <- gs$residuals[[layer]]
    if (is.null(res)) {
      cli::cli_abort("No residuals stored: run {.fn fitMore} with {.code keepResiduals = TRUE}.")
    }
    if (!regKey %in% rownames(res) && regKey %in% colnames(res)) res <- t(res)
    regVal <- getRow(res, regKey, cells)
  }

  if (isTRUE(gs$omicType[[layer]] == 1)) {
    cli::cli_alert_info("{.val {regRow}} is binary: its curve is the smoothed proportion of cells with the regulator active.")
  }

  md <- metadata[cells, , drop = FALSE]
  d <- data.frame(
    donor = as.character(md[[donor]]),
    grp   = if (is.null(by)) "all" else if (facet) paste(md[[by[1]]], md[[by[2]]], sep = " | ") else as.character(md[[by]]),
    rowv  = if (facet) as.character(md[[by[1]]]) else NA_character_,
    colv  = if (facet) as.character(md[[by[2]]]) else NA_character_,
    reg   = regVal,
    tgt   = tgtVal,
    stringsAsFactors = FALSE
  )
  if (identical(orderBy, "regulator")) {
    d$ord <- d$reg
  } else if (identical(orderBy, "target")) {
    d$ord <- d$tgt
  } else {
    if (!orderBy %in% colnames(md)) cli::cli_abort("{.arg orderBy} must be {.val regulator}, {.val target} or a column of {.arg metadata}.")
    d$ord <- as.numeric(md[[orderBy]])
  }

  # ~~~~~~~~~~ 3. Blocks, filtering and ordering ~~~~~~~~~~ #
  keep <- rep(TRUE, nrow(d))
  for (v in names(at)) {
    if (!v %in% colnames(md)) cli::cli_abort("{.val {v}} in {.arg at} is not a column of {.arg metadata}.")
    miss <- setdiff(as.character(at[[v]]), as.character(md[[v]]))
    if (length(miss) > 0) {
      cli::cli_abort("Level(s) {.val {miss}} not found in {.field {v}}.")
    }
    keep <- keep & as.character(md[[v]]) %in% as.character(at[[v]])
  }
  d <- d[keep, , drop = FALSE]
  lvl <- function(v) {
    if (v %in% names(at)) {
      as.character(at[[v]])
    } else if (is.factor(md[[v]])) {
      levels(droplevels(md[[v]][keep]))
    } else {
      sort(unique(as.character(md[[v]][keep])))
    }
  }
  blockLevels <- if (is.null(by)) {
    "all"
  } else if (facet) {
    as.vector(t(outer(lvl(by[1]), lvl(by[2]), paste, sep = " | ")))  # row-major
  } else {
    lvl(by)
  }
  d$grp <- factor(d$grp, levels = blockLevels)
  d <- d[!is.na(d$grp), , drop = FALSE]
  d <- d[stats::ave(seq_len(nrow(d)), d$grp, d$donor, FUN = length) >= minCells, , drop = FALSE]
  if (!nrow(d)) cli::cli_abort("No donor has at least {minCells} cells in the selected blocks.")
  d$grp <- droplevels(d$grp)
  if (facet) {
    d$rowv <- factor(d$rowv, levels = intersect(lvl(by[1]), unique(d$rowv)))
    d$colv <- factor(d$colv, levels = intersect(lvl(by[2]), unique(d$colv)))
  }

  # Report what is actually plotted, to check that the filter was applied
  for (g in levels(d$grp)) {
    dg  <- d[d$grp == g, , drop = FALSE]
    msg <- sprintf("Block %s: %d cells from %d donors", g, nrow(dg), length(unique(dg$donor)))
    cli::cli_alert_info("{msg}")
  }

  # ~~~~~~~~~~ 3b. Model slopes and degenerate blocks ~~~~~~~~~~ #
  slTab <- NULL
  if (annotateSlopes && !is.null(by) && requireNamespace("emmeans", quietly = TRUE)) {
    slTab <- tryCatch({
      emm <- emmeans::as.emmGrid(fit$Models[[tgtKey]]$fit[[regKey]])
      # Restrict the stored slopes to the filtered cells (same logic as `at`)
      grid <- emm@grid
      notInModel <- setdiff(names(at), names(grid))
      if (length(notInModel) > 0) {
        msg <- paste(notInModel, collapse = ", ")
        cli::cli_alert_warning("{msg} is not in the model design: the annotated slope is the same for all its levels.")
      }
      rows <- rep(TRUE, nrow(grid))
      for (v in intersect(names(at), names(grid))) {
        rows <- rows & as.character(grid[[v]]) %in% as.character(at[[v]])
      }
      if (!all(rows)) emm <- emm[which(rows)]
      s   <- as.data.frame(summary(emmeans::emmeans(emm, by),
                                   infer = c(TRUE, TRUE), null = 0))
      skip <- c(by, "SE", "df", "asymp.LCL", "asymp.UCL", "lower.CL", "upper.CL",
                "z.ratio", "t.ratio", "p.value")
      est <- setdiff(names(s), skip)[1]
      data.frame(grp = if (facet) paste(s[[by[1]]], s[[by[2]]], sep = " | ") else as.character(s[[by]]),
                 est = s[[est]],
                 p = if ("p.value" %in% names(s)) s$p.value else NA_real_,
                 hasP = "p.value" %in% names(s))
    }, error = function(e) NULL)
    if (is.null(slTab)) {
      cli::cli_alert_warning("The slopes of this pair could not be read from {.arg fit}; they are not annotated.")
    }
  }

  # A block has no valid slope if the regulator or the target is constant in
  # it, or if the model slope or its p-value is NA/NaN
  diagBlock <- function(g) {
    dg <- d[d$grp == g, , drop = FALSE]
    i  <- if (is.null(slTab)) NA_integer_ else match(g, slTab$grp)
    slopeBad <- !is.null(slTab) &&
      (is.na(i) || !is.finite(slTab$est[i]) || (slTab$hasP[i] && !is.finite(slTab$p[i])))
    regConst <- isTRUE(diff(range(dg$reg)) == 0)
    tgtConst <- isTRUE(diff(range(dg$tgt)) == 0)
    if (!(regConst || tgtConst || slopeBad)) return(list(short = "", long = ""))

    # Presence (value != 0) of the regulator and the target in the block
    n    <- nrow(dg)
    nzR  <- sum(dg$reg != 0)
    nzT  <- sum(dg$tgt != 0)
    nDon <- length(unique(dg$donor))
    nVar <- sum(vapply(split(dg$reg, dg$donor), function(z) diff(range(z)) > 0, logical(1)))
    cause <- if (regConst) "regulator constant" else if (tgtConst) "target constant" else "slope or p-value is NA/NaN"
    list(short = "not estimable",
         long  = sprintf("%s. %d cells; regulator != 0 in %d (%.1f%%), target != 0 in %d (%.1f%%); regulator varies in %d of %d donors",
                         cause, n, nzR, 100 * nzR / n, nzT, 100 * nzT / n, nVar, nDon))
  }
  diagList  <- lapply(stats::setNames(levels(d$grp), levels(d$grp)), diagBlock)
  badShort  <- vapply(diagList, `[[`, character(1), "short")
  badLong   <- vapply(diagList, `[[`, character(1), "long")
  badBlocks <- names(badShort)[nzchar(badShort)]
  if (length(badBlocks) > 0) {
    head <- if (removeNA) "Block(s) without a valid slope, removed from the plot:" else
      "Block(s) without a valid slope, kept in the plot (use removeNA = TRUE to remove them):"
    cli::cli_alert_warning("{head}")
    for (g in badBlocks) {
      msg <- paste0(g, ": ", badLong[[g]])
      cli::cli_alert_info("{msg}")
    }
    if (removeNA) {
      d <- d[!d$grp %in% badBlocks, , drop = FALSE]
      if (!nrow(d)) cli::cli_abort("No block has a valid slope, nothing to plot.")
      d$grp <- droplevels(d$grp)
    }
  }

  set.seed(seed)
  d <- d[order(d$grp, d$ord, stats::runif(nrow(d))), , drop = FALSE]
  d$x <- seq_len(nrow(d))
  d$isZero <- FALSE

  # Display only: zeros of the ordering variable share one x position per block
  if (collapseZeros) {
    if (!any(d$ord == 0)) {
      cli::cli_alert_info("No zeros in the ordering variable: {.arg collapseZeros} has no effect.")
    }
    start <- 0
    for (g in levels(d$grp)) {
      idx <- which(d$grp == g)
      z   <- d$ord[idx] == 0
      if (any(z) && min(d$ord[idx]) >= 0) {
        gap <- max(1, round(0.02 * length(idx)))
        d$isZero[idx] <- z
        d$x[idx[z]]  <- start + 1
        d$x[idx[!z]] <- start + 1 + gap + seq_len(sum(!z))
      } else {
        d$x[idx] <- start + seq_along(idx)
      }
      start <- max(d$x[idx])
    }
  }

  # In a grid every panel is one block: x is rescaled to [0, 1] inside each one
  if (facet) {
    for (g in levels(d$grp)) {
      idx <- which(d$grp == g)
      r   <- range(d$x[idx])
      d$x[idx] <- (d$x[idx] - r[1]) / max(1, diff(r))
    }
  }

  # The target is linearly rescaled to the range of the regulator; the right
  # axis undoes the transformation
  rr <- range(d$reg); tr <- range(d$tgt)
  sc <- if (diff(rr) > 0 && diff(tr) > 0) diff(rr) / diff(tr) else 1
  toLeft  <- function(t) (t - tr[1]) * sc + rr[1]
  toRight <- function(y) (y - rr[1]) / sc + tr[1]

  long <- rbind(
    data.frame(x = d$x, y = d$reg,         var = "regulator", donor = d$donor, grp = d$grp, isZero = d$isZero, rowv = d$rowv, colv = d$colv),
    data.frame(x = d$x, y = toLeft(d$tgt), var = "target",    donor = d$donor, grp = d$grp, isZero = d$isZero, rowv = d$rowv, colv = d$colv)
  )
  long$var <- factor(long$var, levels = c("regulator", "target"))
  long$donorGroup  <- interaction(long$donor, long$grp, long$var, drop = TRUE)
  long$globalGroup <- interaction(long$grp, long$var, drop = TRUE)

  zeroLong <- NULL
  zeroInfo <- NULL
  if (collapseZeros && any(d$isZero)) {
    zs <- d[d$isZero, , drop = FALSE]
    zeroInfo <- do.call(rbind, lapply(split(zs, zs$grp, drop = TRUE), function(z) {
      data.frame(grp = as.character(z$grp[1]), rowv = z$rowv[1], colv = z$colv[1],
                 x = z$x[1], reg = mean(z$reg),
                 tgt = mean(z$tgt),
                 frac = nrow(z) / sum(d$grp == z$grp[1]))
    }))
    zeroLong <- rbind(
      data.frame(x = zeroInfo$x, y = zeroInfo$reg,          var = "regulator", rowv = zeroInfo$rowv, colv = zeroInfo$colv),
      data.frame(x = zeroInfo$x, y = toLeft(zeroInfo$tgt),  var = "target",    rowv = zeroInfo$rowv, colv = zeroInfo$colv)
    )
    zeroLong$var <- factor(zeroLong$var, levels = c("regulator", "target"))
  }

  # ~~~~~~~~~~ 4. Block labels with the slope estimated by the model ~~~~~~~~~~ #
  lab <- unname(vapply(levels(d$grp), function(g) {
    head <- if (facet) character(0) else g   # in a grid the strip already names the block
    body <- if (nzchar(badShort[[g]])) {
      badShort[[g]]
    } else {
      i <- if (is.null(slTab)) NA_integer_ else match(g, slTab$grp)
      if (is.na(i)) character(0) else sprintf("slope = %.3f (p = %.2g)", slTab$est[i], slTab$p[i])
    }
    paste(c(head, body), collapse = "\n")
  }, character(1)))
  if (!is.null(zeroInfo)) {
    for (i in seq_len(nrow(zeroInfo))) {
      j   <- match(zeroInfo$grp[i], levels(d$grp))
      add <- sprintf("%.0f%% zeros collapsed", 100 * zeroInfo$frac[i])
      lab[j] <- if (nzchar(lab[j])) paste0(lab[j], "\n", add) else add
    }
  }
  bnd <- do.call(rbind, lapply(seq_along(levels(d$grp)), function(i) {
    idx <- which(d$grp == levels(d$grp)[i])
    xs  <- d$x[idx]
    data.frame(xmax = max(xs), mid = mean(range(xs)), label = lab[i],
               rowv = d$rowv[idx][1], colv = d$colv[idx][1])
  }))

  # ~~~~~~~~~~ 5. Plot ~~~~~~~~~~ #
  smoothArgs <- list(geom = "line", method = "gam", se = FALSE,
                     formula = y ~ s(x, bs = "cs", k = k))

  p <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y, colour = var))
  if (showPoints) {
    p <- p + ggplot2::geom_point(data = long[!long$isZero, , drop = FALSE],
                                 alpha = pointAlpha, size = 0.8, shape = 16)
  }
  p <- p +
    do.call(ggplot2::stat_smooth,
            c(list(mapping = ggplot2::aes(group = donorGroup),
                   linewidth = donorLinewidth, alpha = donorAlpha), smoothArgs)) +
    do.call(ggplot2::stat_smooth,
            c(list(mapping = ggplot2::aes(group = globalGroup),
                   linewidth = globalLinewidth), smoothArgs))

  if (!is.null(zeroLong)) {
    p <- p + ggplot2::geom_point(data = zeroLong,
                                 ggplot2::aes(x = x, y = y, colour = var),
                                 shape = 18, size = 4, inherit.aes = FALSE)
  }

  if (!facet && nrow(bnd) > 1) {
    p <- p + ggplot2::geom_vline(xintercept = utils::head(bnd$xmax, -1) + 0.5,
                                 linetype = "dashed")
  }

  if (is.null(title)) {
    title <- sprintf("%s ~ %s (%s)", tgtRow, regRow, layer)
    if (length(at) > 0) {
      sub <- paste(vapply(names(at), function(v) {
        paste0(v, " = ", paste(at[[v]], collapse = "/"))
      }, character(1)), collapse = ", ")
      title <- paste0(title, "\n", sub)
    }
  }
  xlab <- sprintf("Cells ordered by %s%s", orderBy,
                  if (is.null(by)) "" else if (facet) sprintf(" within each panel (rows: %s, columns: %s)", by[1], by[2]) else paste0(" within ", by))

  txt <- if (facet) {
    ggplot2::geom_text(data = bnd, ggplot2::aes(x = mid, y = Inf, label = label),
                       inherit.aes = FALSE, vjust = 1.3, size = 3.5)
  } else {
    ggplot2::annotate("text", x = bnd$mid, y = Inf, label = bnd$label,
                      vjust = 1.3, size = 3.5)
  }

  out <- p +
    txt +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.01)) +
    ggplot2::scale_y_continuous(
      name     = regRow,
      expand   = ggplot2::expansion(mult = c(0.05, 0.2)),
      sec.axis = ggplot2::sec_axis(~ toRight(.), name = tgtRow)
    ) +
    ggplot2::labs(x = xlab, title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position    = "none",
      plot.title         = ggplot2::element_text(hjust = 0.5),
      axis.text.x        = ggplot2::element_blank(),
      axis.ticks.x       = ggplot2::element_blank(),
      axis.title.y       = ggplot2::element_text(colour = cols[["regulator"]]),
      axis.text.y        = ggplot2::element_text(colour = cols[["regulator"]]),
      axis.title.y.right = ggplot2::element_text(colour = cols[["target"]]),
      axis.text.y.right  = ggplot2::element_text(colour = cols[["target"]])
    )
  if (facet) {
    out <- out + ggplot2::facet_grid(rows = ggplot2::vars(rowv), cols = ggplot2::vars(colv))
  }
  out
}
