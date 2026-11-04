# =============================================================================
# contrastPatterns(): what does each contrast mean?
#
# A contrast of slopes ("ctrl - treat") only says whether the difference is positive or negative.
# The same sign can come from very different situations: +0.6 against -0.6 (sign switch), +0.8
# against +0.3 (stronger activation), -0.3 against -0.8 (weaker repression), +0.6 against 0 (gain of
# activation)... contrastPatterns() puts next to each contrast the slope of every compared level,
# whether it is detected, a pattern and a sentence. The plots and the networks of a contrast read
# that table.
# =============================================================================

.pattern_levels <- c("sign switch", "activation in both", "repression in both", "detected in one level only", "none detected")
.int_levels     <- c("opposite changes", "same direction", "change detected in one level only", "no change detected in either")

# ---- Classification ----------------------------------------------------------------------------

# Pattern of a pair of slopes, from their states (+, -, 0) and their estimates. `sig` says whether the
# contrast is significant: the text only says "stronger in ..." if it is, and otherwise warns that the
# difference is not significant (the difference between a significant and a non-significant slope is
# not itself significant).
.classify <- function(s1, s2, e1, e2, side1, side2, sig) {
  n <- length(s1)
  pattern <- rep(NA_character_, n); text <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    a <- s1[i]; b <- s2[i]
    if (is.na(a) || is.na(b)) next
    ns <- if (isTRUE(sig[i])) "" else " (diff. not significant)"
    if (a == "+" && b == "+") {
      pattern[i] <- "activation in both"
      text[i] <- if (isTRUE(sig[i])) sprintf("Activation in both; stronger in %s", if (e1[i] > e2[i]) side1[i] else side2[i]) else "Activation in both; no significant difference"
    } else if (a == "-" && b == "-") {
      pattern[i] <- "repression in both"
      text[i] <- if (isTRUE(sig[i])) sprintf("Repression in both; stronger in %s", if (e1[i] < e2[i]) side1[i] else side2[i]) else "Repression in both; no significant difference"
    } else if (a == "+" && b == "-") {
      pattern[i] <- "sign switch"; text[i] <- paste0(sprintf("Activation in %s, repression in %s", side1[i], side2[i]), ns)
    } else if (a == "-" && b == "+") {
      pattern[i] <- "sign switch"; text[i] <- paste0(sprintf("Repression in %s, activation in %s", side1[i], side2[i]), ns)
    } else if (a == "0" && b == "0") {
      pattern[i] <- "none detected"; text[i] <- paste0("No regulation detected in either level", ns)
    } else {                                                           # detected in one side, not in the other
      yes <- if (a != "0") side1[i] else side2[i]; no <- if (a != "0") side2[i] else side1[i]; st <- if (a != "0") a else b
      pattern[i] <- "detected in one level only"
      text[i] <- paste0(sprintf("%s detected in %s, not detected in %s", if (st == "+") "Activation" else "Repression", yes, no), ns)
    }
  }
  data.frame(pattern = factor(pattern, levels = .pattern_levels), interpretation = text, stringsAsFactors = FALSE)
}

# Pattern of an interaction contrast: how the two simple effects relate. A simple effect is the
# difference of slopes between level `a` and level `b` of the first factor, inside one level of the
# second factor. "+" = higher slope in `a`, "-" = higher in `b`, "0" = no difference detected.
.classifyInteraction <- function(st1, st2, d1, d2, a, b, lvl1, lvl2, sig) {
  n <- length(st1)
  pattern <- rep(NA_character_, n); text <- rep(NA_character_, n)
  for (i in seq_len(n)) {
    u <- st1[i]; v <- st2[i]
    if (is.na(u) || is.na(v) || is.na(a[i])) next
    hi <- function(st) if (st == "+") a[i] else b[i]
    lo <- function(st) if (st == "+") b[i] else a[i]
    ns <- if (isTRUE(sig[i])) "" else " (interaction not significant)"
    if (u != "0" && v != "0" && u == v) {
      pattern[i] <- "same direction"
      big <- if (abs(d1[i]) > abs(d2[i])) lvl1[i] else lvl2[i]
      text[i] <- if (isTRUE(sig[i])) sprintf("Slope higher in %s than in %s in both %s and %s; larger difference in %s", hi(u), lo(u), lvl1[i], lvl2[i], big)
      else sprintf("Slope higher in %s than in %s in both %s and %s; similar difference (interaction not significant)", hi(u), lo(u), lvl1[i], lvl2[i])
    } else if (u != "0" && v != "0") {
      pattern[i] <- "opposite changes"
      text[i] <- paste0(sprintf("Slope higher in %s than in %s in %s, but higher in %s than in %s in %s", hi(u), lo(u), lvl1[i], hi(v), lo(v), lvl2[i]), ns)
    } else if (u == "0" && v == "0") {
      pattern[i] <- "no change detected in either"
      text[i] <- paste0(sprintf("No difference between %s and %s detected in either %s or %s", a[i], b[i], lvl1[i], lvl2[i]), ns)
    } else {
      st <- if (u != "0") u else v; yes <- if (u != "0") lvl1[i] else lvl2[i]; no <- if (u != "0") lvl2[i] else lvl1[i]
      pattern[i] <- "change detected in one level only"
      text[i] <- paste0(sprintf("Difference between %s and %s detected in %s (slope higher in %s), not detected in %s", a[i], b[i], yes, hi(st), no), ns)
    }
  }
  data.frame(pattern = factor(pattern, levels = .int_levels), interpretation = text, stringsAsFactors = FALSE)
}

# One key per row, to match the rows of `res` with the rows of its components
.contrastKey <- function(d, cols) do.call(paste, c(lapply(d[cols], as.character), sep = "\r"))

# ---- Main function -------------------------------------------------------------------------------

#' Interpret the contrasts of contrastRegulation()
#'
#' A contrast of slopes (\code{"ctrl - treat"}) only says whether the difference is positive
#' or negative, and the same sign can come from very different situations: \eqn{+0.6} against
#' \eqn{-0.6} (a sign switch), \eqn{+0.8} against \eqn{+0.3} (a stronger activation), \eqn{+0.6}
#' against \eqn{0} (a gain of activation)... \code{contrastPatterns()} reads the slope of each
#' compared level, which \code{contrastRegulation()} keeps in its attribute \code{components},
#' and puts it next to the contrast with its state, a pattern and a sentence. It does not need
#' the \code{fit}, and it does not remove or reorder any row.
#'
#' @section States:
#' A slope is \code{"+"} or \code{"-"} when its confidence interval (1 - alpha, \strong{without}
#' multiplicity correction) excludes 0, and \code{"0"} when it includes 0. \code{"0"} is
#' \strong{not} "equal to zero": it is "not distinguishable from zero", and it depends on how
#' many cells and how much signal the group has. That is why the sentences say "detected" and
#' "not detected", and why a pattern such as "detected in one level only" is only a difference
#' between the two levels if the contrast itself is significant (otherwise the sentence adds
#' "diff. not significant").
#'
#' @section Patterns of a contrast between two levels:
#' From the states of the slope in side 1 and in side 2:
#' \tabular{ll}{
#'   \code{sign switch} \tab \code{+/-} or \code{-/+} \cr
#'   \code{activation in both} \tab \code{+/+} \cr
#'   \code{repression in both} \tab \code{-/-} \cr
#'   \code{detected in one level only} \tab \code{+/0}, \code{0/+}, \code{-/0} or \code{0/-} \cr
#'   \code{none detected} \tab \code{0/0}
#' }
#'
#' @section Patterns of an interaction contrast:
#' The two simple effects (each the difference of slopes between two levels of the first factor
#' inside a level of the second) and how they relate. Only interactions with exactly two simple
#' effects are classified; with more (for example \code{del.eff}) \code{res} is returned
#' unchanged, with a warning.
#' \tabular{ll}{
#'   \code{opposite changes} \tab \code{+/-} or \code{-/+} \cr
#'   \code{same direction} \tab \code{+/+} or \code{-/-} \cr
#'   \code{change detected in one level only} \tab \code{+/0}, \code{0/+}, \code{-/0} or \code{0/-} \cr
#'   \code{no change detected in either} \tab \code{0/0}
#' }
#' For contrasts of a level against a mean of others (\code{eff}, \code{del.eff}) the second
#' side is a mean of several levels, which can hide levels with different signs.
#'
#' @section Columns that are added:
#' For a contrast between two levels: \code{side_1}, \code{side_2} (the levels), \code{slope_1},
#' \code{slope_2}, \code{ci_low_1}, \code{ci_high_1}, \code{ci_low_2}, \code{ci_high_2},
#' \code{p_slope_1}, \code{p_slope_2}, \code{state_1}, \code{state_2}, \code{detection_1},
#' \code{detection_2}, \code{contrast_scale}, \code{pattern}, \code{interpretation},
#' \code{slopes_match} (whether the contrast really is the difference of the two slopes) and
#' \code{detection_warning} (the two levels have very different amounts of data, see
#' \code{detectionRatio}).
#'
#' For an interaction: \code{effect_1}, \code{effect_2}, \code{diff_1}, \code{diff_2} with their
#' intervals (\code{diff_1_low}, \code{diff_1_high}...), \code{p_diff_1}, \code{p_diff_2},
#' \code{state_diff_1}, \code{state_diff_2}, the levels (\code{level_a} to \code{level_d}),
#' \code{pattern}, \code{interpretation}, \code{slopes_text}, \code{detection_min},
#' \code{detection_max}, \code{slopes_match} and \code{detection_warning}.
#'
#' @param res Output of \code{contrastRegulation()} with \code{components = TRUE}: it must carry
#'   the attributes \code{args} and \code{components}. Filtering rows (\code{res[rows, ]},
#'   \code{head()}, \code{dplyr::filter()}) keeps them, but selecting columns
#'   (\code{res[, cols]}, \code{subset()}, \code{dplyr::select()}) or \code{merge()} drops them:
#'   call \code{contrastPatterns()} before that.
#' @param detectionRatio Flag (\code{detection_warning}) the pairs in which the larger of the
#'   detections of the two levels is more than this many times the smaller one.
#' @return \code{res} with the columns above, and the attributes \code{args} and
#'   \code{components} kept.
#' @examples
#' \dontrun{
#' res <- contrastRegulation(fit, compare = "Strain", by = "AnnotClass")
#' x   <- contrastPatterns(res)
#' table(x$pattern, x$significant)
#' x[x$pattern == "sign switch" & x$significant, c("target", "regulator", "interpretation")]
#' summarizePatterns(x)
#' }
#' @seealso \code{\link{summarizePatterns}} for the counts, \code{\link{contrastNetwork}} and
#'   \code{\link{plotSlopeDumbbell}} to draw the result.
#' @family contrast interpretation
#' @export
contrastPatterns <- function(res, detectionRatio = 3) {
  a  <- attr(res, "args"); cm <- attr(res, "components")
  if (is.null(a) || is.null(cm)) cli::cli_abort("{.arg res} must be the output of {.fn contrastRegulation} with {.code components = TRUE}.")
  by   <- a$by
  keys <- c("target", "regulator", by, "contrast")
  kr   <- .contrastKey(res, keys)
  det_ratio <- function(lo, hi) !is.na(lo) & !is.na(hi) & lo > 0 & hi / lo > detectionRatio

  if (!isTRUE(a$interaction)) {
    # Two sides: the weighted means of the slopes with positive and with negative coefficient
    s <- lapply(c("side_1", "side_2"), function(id) { d <- cm[cm$component == id, ]; d[match(kr, .contrastKey(d, keys)), ] })
    for (i in 1:2) {
      d <- s[[i]]
      res[[paste0("side_", i)]]    <- d$label;     res[[paste0("slope_", i)]]   <- d$estimate
      res[[paste0("ci_low_", i)]]  <- d$conf.low;  res[[paste0("ci_high_", i)]] <- d$conf.high
      res[[paste0("p_slope_", i)]] <- d$p.value;   res[[paste0("state_", i)]]   <- d$state
      res[[paste0("detection_", i)]] <- d$detection
    }
    res$contrast_scale <- s[[1]]$coef
    cls <- .classify(res$state_1, res$state_2, res$slope_1, res$slope_2, res$side_1, res$side_2, res$significant)
    res$pattern <- cls$pattern; res$interpretation <- cls$interpretation
    chk <- res$estimate - res$contrast_scale * (res$slope_1 - res$slope_2)
    res$slopes_match <- is.na(chk) | abs(chk) < 1e-6
    res$detection_warning <- det_ratio(pmin(res$detection_1, res$detection_2), pmax(res$detection_1, res$detection_2))
  } else {
    # Interaction: the simple effects (differences of slopes inside a level of the second factor) and the cells.
    # Only contrasts with exactly two simple effects are classified (pairs of levels). With more (for example
    # del.eff or eff on the second factor) the classification would use the wrong effects.
    n_eff <- tapply(cm$kind == "effect", .contrastKey(cm, keys), sum)
    if (any(n_eff > 2, na.rm = TRUE)) {
      cli::cli_alert_warning("Interaction contrasts with more than two simple effects (for example {.code del.eff}) are not classified yet; {.arg res} is returned unchanged. Use {.fn plotFocalInteraction} to see them.")
      return(res)
    }
    e <- lapply(c("effect_1", "effect_2"), function(id) { d <- cm[cm$component == id, ]; d[match(kr, .contrastKey(d, keys)), ] })
    for (i in 1:2) {
      d <- e[[i]]
      res[[paste0("effect_", i)]]  <- d$label;     res[[paste0("diff_", i)]]      <- d$estimate
      res[[paste0("diff_", i, "_low")]] <- d$conf.low; res[[paste0("diff_", i, "_high")]] <- d$conf.high
      res[[paste0("p_diff_", i)]]  <- d$p.value;   res[[paste0("state_diff_", i)]] <- d$state
    }
    pr <- function(lab) { m <- regmatches(lab, regexec("^(.+) - (.+) in (.+)$", lab)); t(vapply(m, function(z) if (length(z) == 4) z[2:4] else rep(NA_character_, 3), character(3))) }
    p1 <- pr(res$effect_1); p2 <- pr(res$effect_2)
    res$level_a <- p1[, 1]; res$level_b <- p1[, 2]; res$level_c <- p1[, 3]; res$level_d <- p2[, 3]
    cls <- .classifyInteraction(res$state_diff_1, res$state_diff_2, res$diff_1, res$diff_2, res$level_a, res$level_b, res$level_c, res$level_d, res$significant)
    res$pattern <- cls$pattern; res$interpretation <- cls$interpretation
    cells <- cm[cm$kind == "cell", ]; ck <- .contrastKey(cells, keys)
    res$slopes_text <- vapply(kr, function(k) { d <- cells[ck == k, ]; paste(sprintf("%s %+.2f", d$label, d$estimate), collapse = ", ") }, character(1))
    res$detection_min <- vapply(kr, function(k) { v <- cells$detection[ck == k]; if (all(is.na(v))) NA_real_ else min(v, na.rm = TRUE) }, numeric(1))
    res$detection_max <- vapply(kr, function(k) { v <- cells$detection[ck == k]; if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE) }, numeric(1))
    chk <- res$estimate - (e[[1]]$coef * res$diff_1 + e[[2]]$coef * res$diff_2)
    res$slopes_match <- is.na(chk) | abs(chk) < 1e-6
    res$detection_warning <- det_ratio(res$detection_min, res$detection_max)
  }
  if (any(!res$slopes_match)) {
    cli::cli_alert_warning("{sum(!res$slopes_match)} row(s): the contrast does not equal the combination of its components (different groups averaged, see {.arg onlyEstimable}).")
  }
  attr(res, "args") <- a; attr(res, "components") <- cm
  res
}

# ---- Summary ---------------------------------------------------------------------------------------

#' How many pairs follow each pattern
#'
#' Counts the pairs of the output of \code{\link{contrastPatterns}} by pattern, for each
#' contrast and level of \code{by} (for an interaction, for each level of \code{by}, because
#' it has one label per pair of levels).
#'
#' @param x Output of \code{contrastPatterns()}.
#' @param onlySignificant Count only the significant contrasts (default). Without it the
#'   counts include every pair of the table.
#' @return A data.frame with the grouping columns, \code{pattern}, \code{n} (number of pairs)
#'   and \code{pct} (percentage inside its group), ordered by group and then by \code{n}.
#' @examples
#' \dontrun{
#' x <- contrastPatterns(contrastRegulation(fit, compare = "Strain", by = "AnnotClass"))
#' summarizePatterns(x)
#' summarizePatterns(x, onlySignificant = FALSE)
#' }
#' @family contrast interpretation
#' @export
summarizePatterns <- function(x, onlySignificant = TRUE) {
  if (!"pattern" %in% names(x)) cli::cli_abort("{.arg x} must be the output of {.fn contrastPatterns}.")
  d <- if (onlySignificant) x[x$significant %in% TRUE, ] else x
  d <- d[!is.na(d$pattern), , drop = FALSE]
  by  <- attr(x, "args")$by
  grp <- if (isTRUE(attr(x, "args")$interaction)) by else c("contrast", by)             # interaction contrasts have one label per pair of levels
  d$pattern <- droplevels(d$pattern)
  tab <- as.data.frame(table(d[, c(grp, "pattern"), drop = FALSE]), stringsAsFactors = FALSE)
  names(tab)[names(tab) == "Freq"] <- "n"
  tab <- tab[tab$n > 0, , drop = FALSE]
  tab$pct <- if (length(grp) > 0) round(100 * tab$n / stats::ave(tab$n, tab[, grp, drop = FALSE], FUN = sum), 1) else round(100 * tab$n / sum(tab$n), 1)
  tab[if (length(grp) > 0) order(tab[[grp[1]]], -tab$n) else order(-tab$n), , drop = FALSE]
}
