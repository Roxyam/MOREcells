# Internal helpers for the functions that analyse the output of fitMore()
# (contrastRegulation() and future exploration functions) --------------------------

#' Map every regulator to its regulatory layer
#'
#' @param fit An object of class \code{MOREcellsFit}.
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

#' Grid of fixed-effect combinations stored in a MOREcellsFit
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
#' @param fit An object of class \code{MOREcellsFit} (uses
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
#' @param fit An object of class \code{MOREcellsFit}.
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
#' @param fit An object of class \code{MOREcellsFit}.
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
        s <- testOne(trends[[reg]], reg)
        list(ok = cbind(target = g, regulator = reg, layer = unname(layerOf[reg]), s,
                        stringsAsFactors = FALSE), fail = NULL)
      }, error = function(e) list(ok = NULL, fail = c(g, reg, conditionMessage(e))))
    })
    list(ok   = do.call(rbind, lapply(out, `[[`, "ok")),
         fail = do.call(rbind, lapply(out, `[[`, "fail")))
  })

  res    <- do.call(rbind, lapply(res_list, `[[`, "ok"))
  failed <- do.call(rbind, lapply(res_list, `[[`, "fail"))
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

  list(res = res, failed = failed)
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
