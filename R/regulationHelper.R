# Internal helpers for the functions that analyse the output of fitMore()
# (contrastRegulation() and future exploration functions) --------------------------

#' Map every regulator to its regulatory layer
#'
#' @param fit An object of class \code{scMoreFit}.
#'
#' @return A named character vector: names are regulators and values are the
#'   names of their layers.
#'
#' @keywords internal
layerLookup <- function(fit) {
  priors <- fit$GlobalSummary$priorAssociations
  unlist(lapply(names(priors), function(l) {
    regs <- unique(as.character(priors[[l]]$regulator))
    stats::setNames(rep(l, length(regs)), regs)
  }))
}


#' Rebuild an \code{emmGrid} from the list stored in a model
#'
#' Wrapper of \code{emmeans::as.emmGrid} that also restores the estimability
#' information. \code{as.list} on an \code{emmGrid} drops the slot \code{nbasis},
#' which is what marks the combinations that cannot be estimated, so
#' \code{\link{fitTargetModel}} stores it as the attribute \code{nbasis} of the
#' list. Without it, those combinations would be reported as zero instead of
#' \code{NA}.
#'
#' @param trend List returned by \code{as.list(emm, model.info.slot = TRUE)}, with
#'   the attribute \code{nbasis} when some combinations cannot be estimated.
#'
#' @return An \code{emmGrid} object.
#'
#' @importFrom emmeans as.emmGrid
#' @keywords internal
toEmmGrid <- function(trend) {
  emm <- emmeans::as.emmGrid(trend)
  nb  <- attr(trend, "nbasis")
  if (!is.null(nb)) emm@nbasis <- nb
  emm
}

#' Grid of fixed-effect combinations stored in a scMoreFit
#'
#' Reads the \code{emmGrid} of the first target-regulator pair whose stored
#' slopes can be reconstructed. All pairs share the same grid, so it is used to
#' validate user input (variables and levels) without looping over every pair.
#' Models that cannot be reconstructed (for example rank-deficient models, whose
#' slopes are stored with inconsistent dimensions) are skipped.
#'
#' @param models Named list of models (\code{fit$Models}, possibly subset).
#'
#' @return A data.frame with one row per combination of fixed-effect levels.
#'
#' @importFrom emmeans as.emmGrid
#' @keywords internal
firstTrendGrid <- function(models) {
  last_error <- NULL
  for (m in models) {
    ok <- Filter(Negate(is.null), m$fit)
    if (length(ok) == 0) next
    grid <- tryCatch(toEmmGrid(ok[[1]])@grid,
                     error = function(e) { last_error <<- conditionMessage(e); NULL })
    if (!is.null(grid)) return(grid)
  }
  cli::cli_abort(c(
    "None of the selected models has slopes that can be reconstructed.",
    "i" = if (!is.null(last_error)) "Last error: {last_error}" else "No target has a fitted model."
  ))
}


#' Validate a list of levels against the stored grid
#'
#' @param at Named list of levels, or \code{NULL}.
#' @param fixedEffects Character vector with the fixed effects of the model.
#' @param grid Data.frame returned by \code{\link{firstTrendGrid}}.
#'
#' @return \code{NULL}, invisibly. Aborts with an informative message if a
#'   variable is not a fixed effect or a level does not exist.
#'
#' @keywords internal
checkAt <- function(at, fixedEffects, grid) {
  if (is.null(at)) return(invisible(NULL))
  if (!is.list(at) || is.null(names(at))) {
    cli::cli_abort("{.arg at} must be a named list, for example {.code list(cell_type = \"Tcell\")}.")
  }
  bad_at <- setdiff(names(at), fixedEffects)
  if (length(bad_at) > 0) {
    cli::cli_abort("{.arg at} contains variables that are not fixed effects of the model: {.val {bad_at}}.")
  }
  for (v in names(at)) {
    avail <- unique(as.character(grid[[v]]))
    miss  <- setdiff(as.character(at[[v]]), avail)
    if (length(miss) > 0) {
      cli::cli_abort(c(
        "Level(s) {.val {miss}} not found in {.field {v}}.",
        "i" = "Available levels: {.val {avail}}."
      ))
    }
  }
  invisible(NULL)
}


#' Restrict a trend grid to the requested levels
#'
#' @param emm An \code{emmGrid} reconstructed from the stored slopes.
#' @param at Named list of levels to keep.
#'
#' @return The \code{emmGrid} restricted to the rows that match \code{at}, or
#'   \code{NULL} if no row matches.
#'
#' @keywords internal
subsetTrendGrid <- function(emm, at) {
  keep <- rep(TRUE, nrow(emm@grid))
  for (v in names(at)) {
    keep <- keep & (as.character(emm@grid[[v]]) %in% as.character(at[[v]]))
  }
  if (!any(keep)) return(NULL)
  emm[which(keep)]
}


#' Number of levels of a fixed effect after applying \code{at}
#'
#' @param v Character. Name of the fixed effect.
#' @param at Named list of levels, or \code{NULL}.
#' @param grid Data.frame returned by \code{\link{firstTrendGrid}}.
#'
#' @return An integer.
#'
#' @keywords internal
nLevelsAt <- function(v, at, grid) {
  if (!is.null(at[[v]])) length(unique(at[[v]])) else length(unique(as.character(grid[[v]])))
}


#' Translate sanitized names back to the original identifiers
#'
#' @param res Data.frame with the columns \code{target}, \code{regulator} and
#'   \code{layer}.
#' @param fit An object of class \code{scMoreFit} (uses
#'   \code{GlobalSummary$nameMap}).
#'
#' @return \code{res} with \code{target} and \code{regulator} translated.
#'
#' @keywords internal
restoreOriginalNames <- function(res, fit) {
  nm  <- fit$GlobalSummary$nameMap
  inv <- function(map) stats::setNames(names(map), unname(map))

  if (length(nm$target) > 0) {
    i   <- inv(nm$target)
    hit <- res$target %in% names(i)
    res$target[hit] <- i[res$target[hit]]
  }
  for (l in names(fit$GlobalSummary$priorAssociations)) {
    if (length(nm[[l]]) == 0) next
    i   <- inv(nm[[l]])
    hit <- res$layer == l & res$regulator %in% names(i)
    res$regulator[hit] <- i[res$regulator[hit]]
  }
  res
}


#' Validate the requested regulatory layers
#'
#' @param layers Character vector of layer names, or \code{NULL}.
#' @param fit An object of class \code{scMoreFit}.
#'
#' @return \code{NULL}, invisibly. Aborts if a layer does not exist.
#'
#' @keywords internal
checkLayers <- function(layers, fit) {
  if (is.null(layers)) return(invisible(NULL))
  available  <- names(fit$GlobalSummary$priorAssociations)
  bad_layers <- setdiff(layers, available)
  if (length(bad_layers) > 0) {
    cli::cli_abort("Unknown layer(s): {.val {bad_layers}}. Available: {.val {available}}.")
  }
  invisible(NULL)
}


#' Select the models of the requested targets
#'
#' @param fit An object of class \code{scMoreFit}.
#' @param targets Character vector of targets (sanitized names), or \code{NULL}
#'   for all of them.
#'
#' @return The sublist of \code{fit$Models} for the requested targets. Aborts if
#'   there is none.
#'
#' @keywords internal
selectModels <- function(fit, targets = NULL) {
  models <- fit$Models
  if (!is.null(targets)) models <- models[intersect(targets, names(models))]
  if (length(models) == 0) cli::cli_abort("No targets to analyse.")
  models
}


#' Give common names to the columns returned by emmeans
#'
#' \pkg{emmeans} names columns differently for mixed models (asymptotic z) and
#' for \code{lm} (t with degrees of freedom), so both sets of names are mapped to
#' common ones.
#'
#' @param res Data.frame built from \code{summary()} of \code{emmeans} objects.
#'
#' @return \code{res} with the columns \code{std.error}, \code{conf.low},
#'   \code{conf.high} and \code{statistic}, and without the \code{null} column.
#'
#' @keywords internal
tidyEmmeansColumns <- function(res) {
  rename <- c(SE = "std.error", asymp.LCL = "conf.low", lower.CL = "conf.low",
              asymp.UCL = "conf.high", upper.CL = "conf.high",
              z.ratio = "statistic", t.ratio = "statistic")
  for (old in names(rename)) names(res)[names(res) == old] <- rename[[old]]
  res$null <- NULL
  res
}


#' Apply a test function to every target-regulator pair
#'
#' Loops over the targets (with a progress bar) and, inside, over their
#' regulators. A failure in one pair never stops the loop: it is recorded
#' together with its reason and returned, so the caller can report it.
#'
#' @param models Named list of models (see \code{\link{selectModels}}).
#' @param regulators Optional character vector to restrict the regulators.
#' @param layers Optional character vector to restrict the layers.
#' @param layerOf Named vector regulator -> layer (see \code{\link{layerLookup}}).
#' @param testOne Function with arguments \code{(trend, reg)} that returns a
#'   data.frame for one pair, or fails.
#' @param verb Character. Used in the error message when no pair can be
#'   processed (for example \code{"tested"} or \code{"compared"}).
#'
#' @return A list with \code{res} (data.frame with one block of rows per pair,
#'   preceded by the columns \code{target}, \code{regulator} and \code{layer})
#'   and \code{failed} (data.frame with the pairs that failed and the reason, or
#'   \code{NULL}).
#'
#' @importFrom pbapply pblapply
#' @keywords internal
testPairs <- function(models, regulators, layers, layerOf, testOne, verb = "tested") {

  res_list <- pbapply::pblapply(names(models), function(g) {
    trends <- models[[g]]$fit
    if (is.null(trends)) return(NULL)

    regs <- names(trends)
    if (!is.null(regulators)) regs <- intersect(regs, regulators)
    if (!is.null(layers))     regs <- regs[layerOf[regs] %in% layers]

    out <- lapply(regs, function(reg) {
      if (is.null(trends[[reg]])) return(list(ok = NULL, fail = c(g, reg, "slopes could not be estimated")))
      tryCatch({
        s    <- testOne(trends[[reg]], reg)
        comp <- attr(s, "components")                     # slopes that make up a contrast (contrastRegulation only)
        list(ok   = cbind(target = g, regulator = reg, layer = unname(layerOf[reg]), s,
                          stringsAsFactors = FALSE),
             comp = if (!is.null(comp)) cbind(target = g, regulator = reg, layer = unname(layerOf[reg]), comp,
                                              stringsAsFactors = FALSE),
             fail = NULL)
      }, error = function(e) list(ok = NULL, fail = c(g, reg, conditionMessage(e))))
    })
    list(ok   = do.call(rbind, lapply(out, `[[`, "ok")),
         comp = do.call(rbind, lapply(out, `[[`, "comp")),
         fail = do.call(rbind, lapply(out, `[[`, "fail")))
  })

  res        <- do.call(rbind, lapply(res_list, `[[`, "ok"))
  components <- do.call(rbind, lapply(res_list, `[[`, "comp"))
  failed     <- do.call(rbind, lapply(res_list, `[[`, "fail"))
  if (!is.null(failed)) {
    failed <- as.data.frame(failed, stringsAsFactors = FALSE)
    names(failed) <- c("target", "regulator", "reason")
  }

  if (is.null(res) || nrow(res) == 0) {
    cli::cli_abort(c(
      "No regulator-target pairs could be {verb}.",
      "i" = if (!is.null(failed) && nrow(failed) > 0) "First error: {failed$reason[1]}" else
        "Check {.arg at}, {.arg regulators} and {.arg layers}."
    ))
  }

  list(res = res, failed = failed, components = components)
}


#' Estimable rows of a grid of slopes
#'
#' A row of the stored grid is not estimable when its slope could not be
#' estimated (the group is not usable for the target, the regulator is not
#' estimable in it, or the coefficient was dependent): its estimate is \code{NA}.
#'
#' @param emm An \code{emmGrid} with the slopes of one regulator (see
#'   \code{\link{toEmmGrid}}).
#' @param reg Character. Name of the regulator.
#'
#' @return A logical vector with one element per row of the grid.
#'
#' @keywords internal
estimableRows <- function(emm, reg) {
  s <- suppressMessages(as.data.frame(summary(emm)))
  !is.na(s[[paste0(reg, ".trend")]])
}


#' Rows that take part in a balanced contrast
#'
#' When a contrast averages over other fixed effects, the levels of \code{compare}
#' are compared like with like: an average only keeps the combinations of the
#' averaged fixed effects in which every compared group is estimable. When nothing
#' is averaged over, each estimable row stands on its own, so the contrasts that do
#' not involve a group that could not be estimated are still computed.
#'
#' @param grid Data.frame with the grid of the \code{emmGrid}.
#' @param ok Logical vector with the estimable rows (see
#'   \code{\link{estimableRows}}).
#' @param compare Character vector of the compared fixed effects.
#' @param by Character vector of the fixed effects that define strata.
#' @param fixedEffects Character vector of the fixed effects of the model.
#'
#' @return A logical vector with the rows that are kept.
#'
#' @keywords internal
balancedRows <- function(grid, ok, compare, by, fixedEffects) {
  averaged <- setdiff(fixedEffects, c(compare, by))
  if (length(averaged) == 0) return(ok)
  keyOf <- function(cols) {
    if (length(cols) > 0) do.call(paste, c(lapply(grid[cols], as.character), sep = "\r")) else rep("all", nrow(grid))
  }
  complete <- as.logical(stats::ave(ok, keyOf(by), keyOf(averaged), FUN = all))
  ok & complete
}


#' Average slopes over the rows that are kept
#'
#' Builds the marginal means of a grid of slopes using only the rows in
#' \code{keep}. \code{emmeans::emmeans} cannot do it because it refuses irregular
#' grids (a missing combination of the fixed effects), so the averages are written
#' as linear combinations of the coefficients, with equal weights or with weights
#' proportional to the number of cells (as \code{emmeans} does: the frequency of
#' each combination of the fixed effects that are averaged over). A combination with no row
#' left is reported as not estimable (\code{NA}).
#'
#' @param emm An \code{emmGrid} of slopes built by \code{\link{fitTargetModel}}.
#' @param factors Character vector of the fixed effects of the resulting means
#'   (\code{character(0)} for a single overall mean).
#' @param by Character vector of the fixed effects that define strata; it is
#'   used to flag the combinations with nothing left to average. It must be
#'   contained in \code{factors}.
#' @param keep Logical vector with the rows of the grid that are averaged.
#' @param ok Logical vector with the estimable rows of the grid.
#' @param weights Either \code{"equal"} or \code{"proportional"}.
#'
#' @return An \code{emmGrid} over \code{factors}.
#'
#' @importFrom emmeans emmobj
#' @keywords internal
averageEstimable <- function(emm, factors, by, keep, ok, weights = c("equal", "proportional")) {
  weights <- match.arg(weights)
  grid    <- emm@grid
  rowKey  <- function(d, cols) {
    if (length(cols) > 0) do.call(paste, c(lapply(d[cols], as.character), sep = "\r")) else rep("all", nrow(d))
  }

  if (length(factors) > 0) {
    lv <- lapply(factors, function(f) intersect(emm@levels[[f]], as.character(grid[[f]])))
    names(lv) <- factors
    cells <- expand.grid(lv, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  } else {
    lv    <- list(`1` = "overall")                                  # as emmeans names a single overall mean
    cells <- data.frame(`1` = "overall", check.names = FALSE)
  }
  keyGrid <- rowKey(grid, factors); keyCell <- rowKey(cells, factors)
  byGrid  <- rowKey(grid, by);      byCell  <- rowKey(cells, by)

  # Proportional weights follow emmeans: the frequency of each combination of the fixed effects that are
  # averaged over, counted across all the other combinations (not within every stratum)
  averaged <- setdiff(setdiff(names(grid), ".wgt."), factors)
  w <- if (weights == "equal" || length(averaged) == 0) {
    rep(1, nrow(grid))
  } else {
    stats::ave(grid$.wgt., rowKey(grid, averaged), FUN = sum)
  }
  W     <- matrix(0, nrow(cells), nrow(grid))
  cellN <- numeric(nrow(cells))
  for (i in seq_len(nrow(cells))) {
    rows <- which(keyGrid == keyCell[i] & keep)
    if (length(rows) > 0 && sum(w[rows]) > 0) {
      W[i, rows] <- w[rows] / sum(w[rows])
      cellN[i]   <- sum(grid$.wgt.[rows])
    } else {
      # Nothing to average: the combination points to a row that cannot be estimated, so it is NA.
      # Every empty combination gets a different scale: the contrast between two of them is then
      # still a non-estimable combination (with the same scale it would cancel out to zero)
      bad <- which(!ok & byGrid == byCell[i])
      if (length(bad) == 0) stop("a combination without estimable rows has no row that can be flagged")
      W[i, bad[1]] <- 2^((i - 1) %% 40)
    }
  }

  out <- emmeans::emmobj(bhat = emm@bhat, V = emm@V, levels = lv, linfct = W %*% emm@linfct,
                         dffun = emm@dffun, dfargs = emm@dfargs, estName = emm@misc$estName)
  out@nbasis     <- emm@nbasis
  out@grid$.wgt. <- cellN
  out
}


#' Groups left out of each row of a result
#'
#' For every row of a result (a combination of the levels of \code{by}), lists the
#' groups of the grid that are not in \code{keep}.
#'
#' @param emm An \code{emmGrid} of slopes.
#' @param keep Logical vector with the rows of the grid that were averaged.
#' @param by Character vector of the fixed effects that define the rows of the
#'   result.
#' @param fixedEffects Character vector of the fixed effects of the model.
#' @param s Data.frame with the result of one regulator (it must have the columns
#'   in \code{by}).
#'
#' @return A list with \code{n} (number of groups left out of each row) and
#'   \code{groups} (their names, separated by \code{"; "}, or \code{NA} if none).
#'
#' @keywords internal
leftOutGroups <- function(emm, keep, by, fixedEffects, s) {
  grid  <- emm@grid
  label <- do.call(paste, c(lapply(grid[fixedEffects], as.character), sep = " / "))
  keyOf <- function(d) {
    if (length(by) > 0) do.call(paste, c(lapply(d[by], as.character), sep = "\r")) else rep("all", nrow(d))
  }
  keyGrid <- keyOf(grid)
  out     <- lapply(keyOf(s), function(k) label[keyGrid == k & !keep])
  list(n      = vapply(out, length, integer(1)),
       groups = vapply(out, function(x) if (length(x) > 0) paste(x, collapse = "; ") else NA_character_, character(1)))
}


#' Slopes that make up every contrast
#'
#' A contrast is a linear combination of the slopes of the groups (cells of the
#' grid) with coefficients that add up to zero. This function returns the slopes
#' (or differences of slopes) that make it up, so that the sign of a contrast can
#' be interpreted:
#' \itemize{
#'   \item \code{side}: \code{side_1} and \code{side_2} are the weighted means of
#'     the cells with positive and with negative coefficient. The contrast equals
#'     \code{k * (side_1 - side_2)}, where \code{k} (column \code{coef}) is the sum
#'     of the positive coefficients.
#'   \item \code{cell}: the slope of every cell with a non-zero coefficient. Only
#'     given when a contrast involves more than two cells (for example a level
#'     against the mean of the others) or is an interaction.
#'   \item \code{effect}: for interaction contrasts, the simple effects (the
#'     contrast of the first compared fixed effect inside a level of the second).
#'     The interaction is the sum of the effects weighted by their \code{coef}.
#' }
#'
#' @param sub \code{emmGrid} with the slopes averaged over the fixed effects that are
#'   not in \code{compare} or \code{by}.
#' @param ct \code{emmGrid} with the contrasts, built from \code{sub}.
#' @param emm \code{emmGrid} with the stored slopes of the regulator (one row per
#'   group of the design).
#' @param keep Logical vector with the rows of \code{emm} that were averaged.
#' @param factors Character vector with the fixed effects of \code{sub}
#'   (\code{compare} and \code{by}).
#' @param compare Character vector of the compared fixed effects.
#' @param by Character vector of the fixed effects that define strata.
#' @param interaction Logical. Whether the contrasts are interaction contrasts.
#' @param fixedEffects Character vector of the fixed effects of the model.
#' @param labels Character vector with the tidy label of every contrast, in the order of the
#'   rows of \code{ct}.
#' @param byvals Data.frame with the levels of \code{by} of every contrast (or
#'   \code{NULL}).
#' @param alpha Numeric. Confidence level is \code{1 - alpha}.
#'
#' @return A data.frame with one row per component (see Details), or \code{NULL} if
#'   the contrasts cannot be decomposed.
#'
#' @importFrom emmeans emmobj
#' @keywords internal
contrastComponents <- function(sub, ct, emm, keep, factors, compare, by, interaction,
                               fixedEffects, labels, byvals, alpha) {

  if (isTRUE(interaction) && length(compare) != 2) return(NULL)
  keyOf <- function(d, cols) do.call(paste, c(lapply(d[cols], as.character), sep = "\r"))
  fmt   <- function(x) if (length(x) == 1) x else sprintf("mean(%s)", paste(x, collapse = ", "))
  state <- function(lo, hi) ifelse(is.na(lo) | is.na(hi), NA_character_, ifelse(lo > 0, "+", ifelse(hi < 0, "-", "0")))

  cf <- stats::coef(ct)
  if (is.null(cf)) stop("emmeans returned no contrast coefficients (coef() is NULL): update the emmeans package")
  g  <- sub@grid
  cn <- setdiff(names(cf), c(names(g), ".wgt."))                                   # the columns with the coefficients, whatever their names
  if (length(cn) == 0) stop("coef() of the contrasts has no coefficient columns (columns: ", paste(names(cf), collapse = ", "), ")")
  if (all(factors %in% names(cf))) {
    cf <- cf[match(keyOf(g, factors), keyOf(cf, factors)), , drop = FALSE]        # same order as the cells of `sub`
  } else if (nrow(cf) != nrow(g)) {
    stop("coef() of the contrasts does not match the grid of the means")
  }
  cellLab <- apply(g[compare], 1, paste, collapse = " ")             # the strata (`by`) have their own columns

  # Groups of the design (rows of `emm`) that make up every cell of `sub`
  eg      <- emm@grid
  grpLab  <- do.call(paste, c(lapply(eg[fixedEffects], as.character), sep = " / "))
  cellGroups <- lapply(seq_len(nrow(g)), function(i) {
    m <- keep
    for (f in factors) m <- m & as.character(eg[[f]]) == as.character(g[[f]][i])
    grpLab[m]
  })

  out <- vector("list", length(cn))
  for (k in seq_along(cn)) {
    cvec <- cf[[cn[k]]]
    nz   <- which(abs(cvec) > 1e-12)
    if (length(nz) == 0) next
    W <- list(); info <- list()
    add <- function(kind, comp, label, w, coefv) {
      W[[length(W) + 1]]       <<- w
      info[[length(info) + 1]] <<- data.frame(kind = kind, component = comp, label = label, coef = coefv,
                                              groups = paste(unique(unlist(cellGroups[which(w != 0)])), collapse = "; "),
                                              stringsAsFactors = FALSE)
    }
    cellRows <- integer(0)

    if (!isTRUE(interaction)) {
      pos <- nz[cvec[nz] > 0]; neg <- nz[cvec[nz] < 0]
      w1 <- numeric(length(cvec)); w1[pos] <- cvec[pos] / sum(cvec[pos])
      w2 <- numeric(length(cvec)); w2[neg] <- -cvec[neg] / sum(-cvec[neg])
      kk <- sum(cvec[pos])
      add("side", "side_1", fmt(cellLab[pos]), w1, kk)
      add("side", "side_2", fmt(cellLab[neg]), w2, kk)
      if (length(nz) > 2) cellRows <- nz
    } else {
      # Coefficients of an interaction contrast: a product u (first factor) x v (second factor)
      f1 <- compare[1]; f2 <- compare[2]
      M  <- unclass(stats::xtabs(cvec[nz] ~ as.character(g[[f1]])[nz] + as.character(g[[f2]])[nz]))
      ij <- arrayInd(which.max(abs(M)), dim(M))
      u  <- M[, ij[2]] / max(abs(M[, ij[2]])); v <- M[ij[1], ] / u[ij[1]]
      for (j in which(abs(v) > 1e-12)) {
        w <- numeric(length(cvec))
        for (i in seq_along(u)) {
          r <- which(as.character(g[[f1]]) == rownames(M)[i] & as.character(g[[f2]]) == colnames(M)[j] & seq_along(cvec) %in% nz)
          if (length(r) > 0 && abs(u[i]) > 1e-12) w[r] <- u[i]
        }
        add("effect", paste0("effect_", sum(vapply(info, function(z) z$kind == "effect", logical(1))) + 1),
            sprintf("%s - %s in %s", fmt(rownames(M)[u > 0]), fmt(rownames(M)[u < 0]), colnames(M)[j]), w, v[j])
      }
      cellRows <- nz
    }
    for (i in cellRows) { w <- numeric(length(cvec)); w[i] <- 1; add("cell", cellLab[i], cellLab[i], w, cvec[i]) }

    d <- do.call(rbind, info)
    for (f in compare) d[[f]] <- NA_character_
    for (r in which(d$kind == "cell")) for (f in compare) d[[f]][r] <- as.character(g[[f]][cellRows[r - sum(d$kind != "cell")]])
    bv <- NULL
    if (!is.null(byvals)) { bv <- byvals[k, , drop = FALSE]; rownames(bv) <- NULL }
    out[[k]] <- list(W = do.call(rbind, W), info = d, label = labels[k], by = bv)
  }
  out <- Filter(Negate(is.null), out)
  if (length(out) == 0) return(NULL)

  # Estimate, confidence interval and test of every component of every contrast in one go: they are all
  # linear combinations of the cells of `sub`
  L   <- do.call(rbind, lapply(out, `[[`, "W")) %*% sub@linfct
  obj <- emmeans::emmobj(bhat = sub@bhat, V = sub@V, levels = list(comp = seq_len(nrow(L))), linfct = L,
                         dffun = sub@dffun, dfargs = sub@dfargs, estName = "est")
  obj@nbasis <- sub@nbasis
  s   <- as.data.frame(summary(obj, infer = c(TRUE, TRUE), null = 0, level = 1 - alpha, adjust = "none"))
  lo  <- s[[grep("LCL|lower", names(s))[1]]]; hi <- s[[grep("UCL|upper", names(s))[1]]]

  start <- 0
  res <- lapply(out, function(p) {
    i <- start + seq_len(nrow(p$W)); start <<- start + nrow(p$W)
    pieces <- list(data.frame(contrast = p$label, stringsAsFactors = FALSE), p$by, p$info,
                   data.frame(estimate = s$est[i], std.error = s$SE[i], conf.low = lo[i], conf.high = hi[i],
                              p.value = s$p.value[i], state = state(lo[i], hi[i]), stringsAsFactors = FALSE))
    do.call(cbind, Filter(Negate(is.null), pieces))
  })
  res <- do.call(rbind, res)
  rownames(res) <- NULL
  res
}


#' % of cells with a value of a target in some groups
#'
#' @param fit An object of class \code{scMoreFit}.
#' @param target Character. Target feature.
#' @param groups Character. Groups separated by \code{"; "} (as in the column
#'   \code{groups} of the components of a contrast).
#'
#' @return Numeric. Percentage of cells with a value in the groups in which the
#'   target is usable, or \code{NA}.
#'
#' @keywords internal
detectionOfGroups <- function(fit, target, groups) {
  est <- fit$GlobalSummary$estimability
  if (is.null(est) || !(target %in% rownames(est$targetCounts)) || is.na(groups) || !nzchar(groups)) return(NA_real_)
  gs <- strsplit(groups, "; ", fixed = TRUE)[[1]]
  gs <- gs[gs %in% est$groups]
  gs <- gs[est$targetGroups[target, gs]]
  if (length(gs) == 0) return(NA_real_)
  100 * sum(est$targetCounts[target, gs]) / sum(est$groupSizes[gs])
}
