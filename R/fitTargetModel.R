#' Construct the formula and fit the mixed model of one target feature
#'
#' Prepares the cell-level data of a target feature using only the cells of the
#' groups in which the target is usable, builds the model formula with its prior
#' regulators and the groups of the design, adds random intercepts and fits the
#' model with \code{\link{runGLMM}}. Optionally, the slope of every regulator in
#' every group is read from the coefficients of the model and the fitted model is
#' replaced by a lightweight list of \code{emmeans} objects (one per regulator).
#'
#' @details
#' A group is a combination of the levels of the fixed effects (see
#' \code{\link{cellGroups}}). With \code{interaction = TRUE} the model has one
#' intercept per group and one explicit slope column for every regulator and group
#' in which the regulator is estimable:
#' \code{group + reg_1:group_1 + ... + reg_k:group_g}. It is the same model as
#' \code{(reg_1 + ... + reg_k) * (fe_1 * ... * fe_m)} when every group is used, but
#' it stays well defined when some groups are left out (a crossed formula would
#' have empty cells and \code{lme4} would silently drop coefficients). With
#' \code{interaction = FALSE} it is \code{reg_1 + ... + reg_k + fe_1 + ... + fe_m}
#' (one slope per regulator) and without fixed effects it is
#' \code{reg_1 + ... + reg_k}. The random part is
#' \code{(1 | re_1) + ... + (1 | re_j)}.
#'
#' When \code{asEmmGrid = TRUE} the slopes are read from the coefficients and
#' their covariance matrix and stored as \code{emmGrid} lists over the complete
#' grid of the fixed effects (asymptotic degrees of freedom for mixed models).
#' The combinations that cannot be estimated (groups in which the target is not
#' usable, regulators that are not estimable in a group, groups without cells or
#' dependent coefficients) are \code{NA}: the information is kept in the attribute
#' \code{nbasis} of every list, which \code{\link{toEmmGrid}} restores.
#'
#' @param target Character. Name of the target feature (row of \code{targetData}).
#' @param targetData Matrix, sparse \code{Matrix} or data.frame of target
#'   features (features x cells).
#' @param regData Numeric matrix of regulatory features (regulators x cells),
#'   typically the stacked regulatory layers.
#' @param metadata Data.frame with cell metadata. Row names must be cell IDs.
#' @param regulators Character vector of regulators of \code{target}
#'   (rownames of \code{regData}).
#' @param fixedEffects Character vector of \code{metadata} columns used as
#'  fixed terms. They must be factors.
#' @param randomEffects Character vector of \code{metadata} columns used
#'  as random intercepts.
#' @param interaction Logical. If \code{TRUE}, regulators are crossed
#'   with all fixed effects and their interactions.
#' @param targetOmicType Numeric flag: \code{0} = Gaussian (\code{lmer}),
#'   \code{1} = binary (\code{glmer}).
#' @param REML Logical. Restricted maximum likelihood (Gaussian models only).
#' @param verbose Logical. If \code{TRUE}, prints optimization progress.
#' @param asEmmGrid Logical. If \code{TRUE}, the fitted model is converted into a
#'   named list (one element per regulator) of \code{emmGrid} lists.
#' @param cellGroup Factor named by cell with the group of each cell (see
#'   \code{\link{cellGroups}}). If \code{NULL}, it is computed from
#'   \code{fixedEffects}.
#' @param usableGroups Named logical vector with the groups in which the target is
#'   usable. Only the cells of these groups are used. If \code{NULL}, all groups
#'   are used.
#' @param estimableGroups Logical matrix (regulators x groups) with the groups in
#'   which each regulator is estimable. If \code{NULL}, all groups are estimable.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{fit}: The fitted model, or a named list of \code{emmGrid} lists
#'     per regulator if \code{asEmmGrid = TRUE} (an element is \code{NULL}
#'     if its estimation failed; the whole \code{fit} is \code{NULL} if the model
#'     failed).
#'   \item \code{converged}, \code{is_singular}, \code{warnings}, \code{error}:
#'     Diagnostics returned by \code{\link{runGLMM}}.
#'   \item \code{target}: Identifier of the target feature.
#'   \item \code{formula}: The model formula.
#'   \item \code{regulators}: Regulators included in the model.
#'   \item \code{usedGroups}: Groups used to fit the model.
#'   \item \code{nonEstimable}: Data.frame with the regulator-group slopes that
#'     could not be estimated and the reason (\code{NULL} if all could).
#'   \item \code{n_slopes_possible}, \code{n_slopes_estimable}: Number of
#'     regulator-group slopes (in the groups with cells) and how many were
#'     estimated.
#'   \item \code{slopeMap}: Data.frame that links every slope column of the model
#'     with its regulator and its row in the grid (\code{interaction = TRUE} only).
#' }
#'
#' @importFrom stats as.formula
#' @importFrom emmeans emmobj
#' @keywords internal
fitTargetModel <- function(target,
                           targetData,
                           regData,
                           metadata,
                           regulators,
                           fixedEffects,
                           randomEffects,
                           interaction,
                           targetOmicType,
                           REML,
                           verbose,
                           asEmmGrid,
                           cellGroup       = NULL,
                           usableGroups    = NULL,
                           estimableGroups = NULL) {

  # Result returned when the model cannot be built
  noModel <- function(msg, regs = character(0)) {
    list(fit = NULL, converged = FALSE, is_singular = NA, warnings = NULL, error = msg,
         target = target, formula = NULL, regulators = regs)
  }

  regulators <- intersect(unique(regulators), rownames(regData))
  if (length(regulators) == 0) {
    return(noModel("No regulators of this target are available in `regData`."))
  }

  has_fixed <- length(fixedEffects) > 0
  if (has_fixed && !all(vapply(fixedEffects, function(v) is.factor(metadata[[v]]), logical(1)))) {
    cli::cli_abort("Fixed effects must be factors: convert the numeric ones in {.arg metadata}.")
  }

  # ~~~~~~~~~~ 1. Groups of the design ~~~~~~~~~~ #
  # The grid has every combination of the fixed effects, also those without cells
  cells  <- colnames(targetData)
  if (is.null(cellGroup)) cellGroup <- cellGroups(metadata, fixedEffects, cells)
  groups <- levels(cellGroup)

  xlev       <- if (has_fixed) lapply(metadata[fixedEffects], levels) else list(`1` = "overall")
  gridLabels <- if (has_fixed) {
    do.call(paste, c(expand.grid(xlev, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE), sep = " / "))
  } else {
    "all"
  }
  nGrid   <- length(gridLabels)
  gridRow <- match(groups, gridLabels)

  if (is.null(usableGroups)) usableGroups <- stats::setNames(rep(TRUE, length(groups)), groups)
  if (is.null(names(usableGroups)) && length(usableGroups) == length(groups)) names(usableGroups) <- groups
  usableGroups <- usableGroups[groups]
  if (is.null(estimableGroups)) {
    estimableGroups <- matrix(TRUE, length(regulators), length(groups), dimnames = list(regulators, groups))
  }
  estimableGroups <- estimableGroups[regulators, groups, drop = FALSE]
  usedGroups      <- groups[usableGroups]

  # ~~~~~~~~~~ 2. Cells of the groups used ~~~~~~~~~~ #
  cellsUsed <- cells[as.character(cellGroup[cells]) %in% usedGroups]
  if (length(cellsUsed) == 0) {
    return(noModel("The target is not usable in any group.", regulators))
  }

  # ~~~~~~~~~~ 3. Design matrix construction ~~~~~~~~~~ #
  design_vars <- unique(c(fixedEffects, randomEffects))
  target_expr <- targetData[target, cellsUsed]
  reg_sub     <- t(as.matrix(regData[regulators, cellsUsed, drop = FALSE]))
  md_sub      <- droplevels(metadata[cellsUsed, design_vars, drop = FALSE])
  cellRow     <- gridRow[match(as.character(cellGroup[cellsUsed]), groups)]   # grid row of each cell

  re_terms <- if (length(randomEffects) > 0) {
    paste0("(1 | ", randomEffects, ")", collapse = " + ")
  } else {
    NULL
  }

  # ~~~~~~~~~~ 4. Dynamic formula generation ~~~~~~~~~~ #
  slopeMap <- NULL
  if (interaction && has_fixed) {
    # One slope column for every regulator and group in which it is estimable
    gfac <- factor(cellRow)
    cols <- list()
    map  <- list()
    for (reg in regulators) {
      for (g in sort(unique(cellRow))) {
        if (isTRUE(estimableGroups[reg, gridLabels[g]])) {
          nm <- paste0("slope__", length(cols) + 1)
          cols[[nm]] <- reg_sub[, reg] * (cellRow == g)
          map[[nm]]  <- data.frame(column = nm, regulator = reg, grid = g, stringsAsFactors = FALSE)
        }
      }
    }
    if (length(cols) == 0) {
      return(noModel("No regulator is estimable in the groups where the target is usable.", regulators))
    }
    slopeMap <- do.call(rbind, unname(map))

    model_data <- data.frame(
      TargetExpression = target_expr,
      gfac             = gfac,
      md_sub[, randomEffects, drop = FALSE],
      cols,
      check.names = FALSE
    )
    fe_terms <- c(if (nlevels(gfac) > 1) "gfac", names(cols))
  } else {
    # One slope per regulator (shared by all the groups)
    fixed_used <- fixedEffects[vapply(fixedEffects, function(v) nlevels(md_sub[[v]]) > 1, logical(1))]
    model_data <- data.frame(
      TargetExpression = target_expr,
      reg_sub,
      md_sub[, design_vars, drop = FALSE],
      check.names = FALSE
    )
    fe_terms <- c(regulators, fixed_used)
  }
  if (length(fe_terms) == 0) fe_terms <- "1"

  formula_str   <- paste("TargetExpression ~", paste(c(paste(fe_terms, collapse = " + "), re_terms), collapse = " + "))
  model_formula <- stats::as.formula(formula_str, env = globalenv())

  # ~~~~~~~~~~ 5. Model execution ~~~~~~~~~~ #
  fit_res <- runGLMM(
    formula        = model_formula,
    data           = model_data,
    targetOmicType = targetOmicType,
    REML           = REML,
    verbose        = verbose
  )

  fit_res$target     <- target
  fit_res$formula    <- model_formula
  fit_res$regulators <- regulators
  fit_res$usedGroups <- usedGroups
  fit_res$slopeMap   <- slopeMap

  # ~~~~~~~~~~ 6. Slopes as EmmGrid lists ~~~~~~~~~~ #
  if (asEmmGrid && !is.null(fit_res$fit)) {
    fit   <- fit_res$fit
    is_me <- inherits(fit, "merMod")
    cf    <- if (is_me) lme4::fixef(fit, add.dropped = TRUE) else stats::coef(fit)
    cf    <- cf[!is.na(cf)]                                             # coefficients dropped by the model are NA
    V     <- as.matrix(stats::vcov(fit))[names(cf), names(cf), drop = FALSE]

    # Asymptotic degrees of freedom, except for linear models
    isLm <- inherits(fit, "lm") && !inherits(fit, "glm")
    dfv  <- if (isLm) stats::df.residual(fit) else Inf
    dfa  <- if (isLm) list(df = dfv) else list()
    dfun <- if (isLm) function(k, dfargs) dfargs$df else function(k, dfargs) Inf

    wgt      <- tabulate(cellRow, nbins = nGrid)                        # cells of every group (weights of the grid)
    usedRows <- sort(unique(cellRow))

    # Reason why a slope cannot be estimated
    reasonOf <- function(reg, rows) {
      lab <- gridLabels[rows]
      ifelse(!(lab %in% groups), "No cells in the design",
             ifelse(!(lab %in% usedGroups), "Target not detected in the group",
                    ifelse(!estimableGroups[reg, match(lab, groups)], "Regulator not detected in the group",
                           "Linearly dependent on other terms")))
    }

    trends <- list()
    nonEst <- list()
    nOk    <- 0L
    for (reg in regulators) {
      trends[[reg]] <- tryCatch({
        bh <- numeric(nGrid); Vr <- matrix(0, nGrid, nGrid); ok <- rep(FALSE, nGrid)
        if (!is.null(slopeMap)) {
          m <- slopeMap[slopeMap$regulator == reg & slopeMap$column %in% names(cf), , drop = FALSE]
          if (nrow(m) > 0) {
            bh[m$grid] <- cf[m$column]
            Vr[m$grid, m$grid] <- V[m$column, m$column]
            ok[m$grid] <- TRUE
          }
        } else if (reg %in% names(cf)) {
          bh[usedRows] <- cf[reg]
          Vr[usedRows, usedRows] <- V[reg, reg]
          ok[usedRows] <- TRUE
        }

        emm <- emmeans::emmobj(bh, Vr, levels = xlev, linfct = diag(nGrid), df = dfv,
                               dffun = dfun, dfargs = dfa, estName = paste0(reg, ".trend"))
        emm@grid$.wgt. <- wgt
        if (any(!ok)) emm@nbasis <- diag(nGrid)[, !ok, drop = FALSE]

        l <- as.list(emm, model.info.slot = TRUE)
        # as.list() drops the estimability information: it is stored as an attribute
        if (any(!ok)) attr(l, "nbasis") <- emm@nbasis

        inGroups <- gridLabels %in% groups
        nOk      <- nOk + sum(ok & inGroups)
        if (any(!ok & inGroups)) {
          rows <- which(!ok & inGroups)
          nonEst[[reg]] <- data.frame(target = target, regulator = reg, group = gridLabels[rows],
                                      reason = reasonOf(reg, rows), stringsAsFactors = FALSE)
        }
        l
      }, error = function(e) NULL)
    }

    fit_res$fit                <- trends
    fit_res$nonEstimable       <- if (length(nonEst) > 0) do.call(rbind, unname(nonEst)) else NULL
    fit_res$n_slopes_possible  <- length(regulators) * length(groups)
    fit_res$n_slopes_estimable <- nOk
  }

  return(fit_res)
}
