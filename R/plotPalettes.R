# =============================================================================
# scMore colours
#
# Only what the plots use:
#   modality   colour of each regulator layer (TF, ATAC, miRNA...)
#   node       colour of the targets
#   sign       colour of the edges of a group network (positive / negative slope)
#   owner      colour of the edges of a contrast network (group A, group B, conserved)
#
# Everything can be changed for the session with scMoreSetColors().
# =============================================================================

# The general palette (the colours of your colour list). The first ones are used for the
# layers and the targets; the rest are used for any extra layer.
.scMoreColorList <- c("#98A6D4", "#a9d498", "#d498b9", "#f4e147", "#98ccd4",
                      "#F3B775", "#b598d4", "#d49898", "#F0789A", "#FFFAAA",
                      "#97DAE1", "#D98181", "#FFB3BA", "#C5B0D5", "#C7E9C0",
                      "#FDD0A2", "#F7B6D2", "#DBDB8D", "#AEC7E8", "#D62728",
                      "#C49C94", "#78BAF0")

.scMoreDefaultColors <- list(
  modality = c(TF = "#98A6D4", ATAC = "#a9d498"),                 # 1st and 2nd colour of the list
  node     = c(target = "#d498b9"),                                # 3rd colour of the list
  sign     = c(positive = "#C8525B", negative = "#6479BE"),
  owner    = c(A = "#CE4257", B = "#9DB484", shared = "#BC836D")
)

# changes made with scMoreSetColors() live here (the environment can be modified even
# though the package namespace is locked)
.scMoreEnv <- new.env(parent = emptyenv())
.scMoreEnv$colors <- list()

# A layer is recognised whatever its capitalisation or plural: "TF", "TFs" and "tf" are
# the same layer, and so are "miRNA", "miRNAs" and "mirna".
.layerKey <- function(x) toupper(sub("s$", "", as.character(x), ignore.case = TRUE))

.toHex <- function(colours) {
  m <- tryCatch(grDevices::col2rgb(colours), error = function(e) NULL)
  if (is.null(m)) cli::cli_abort("{.arg values} must be valid colours (names such as {.val red} or hex codes).")
  hex <- grDevices::rgb(m[1, ], m[2, ], m[3, ], maxColorValue = 255)
  names(hex) <- names(colours)
  hex
}

#' Colours of the networks
#'
#' \code{scMoreColors()} returns the colours the network plots use:
#' \describe{
#'   \item{\code{"modality"}}{one colour per regulator layer. TF and ATAC take the first
#'     and the second colour of the colour list. A layer is recognised whatever its
#'     capitalisation or plural (\code{TF}, \code{TFs}, \code{tf}).}
#'   \item{\code{"node"}}{the colour of the targets (\code{target}).}
#'   \item{\code{"sign"}}{the edges of a group network: \code{positive} and \code{negative}
#'     slope.}
#'   \item{\code{"owner"}}{the edges of a contrast network: group \code{A} (the first
#'     level), group \code{B} (the second) and \code{shared} (conserved in both).}
#' }
#' \code{scMoreSetColors()} changes some of them for the session, and
#' \code{scMoreSetColors(kind, NULL)} restores them.
#' \code{scMoreModalityColors()} gives a colour to each layer present in the data: the
#' known layers keep theirs and any other gets the next free colour of the colour list.
#' \code{scMorePalette()} returns the colour list.
#'
#' @param kind \code{"modality"}, \code{"node"}, \code{"sign"} or \code{"owner"}.
#' @param values A named vector of colours replacing or adding entries of \code{kind}.
#'   \code{NULL} restores the defaults for that kind.
#' @param levels Names of the layers in the data.
#' @param n Number of colours of the colour list (interpolated if there are more than 22);
#'   \code{NULL} for all of them.
#' @return A named character vector (\code{scMorePalette()}: an unnamed one).
#' @examples
#' scMoreColors("owner")
#' scMoreSetColors("owner", c(A = "#CE4257", B = "#9DB484"))
#' scMoreSetColors("modality", c(CNV = "#F0789A"))     # a new layer with its own colour
#' scMoreSetColors("owner", NULL)                      # back to the defaults
#' @export
scMoreColors <- function(kind = c("modality", "node", "sign", "owner")) {
  kind <- match.arg(kind)
  out <- .scMoreDefaultColors[[kind]]
  ov  <- .scMoreEnv$colors[[kind]]
  if (is.null(ov)) return(out)
  if (kind == "modality") {                       # an override of "TFs" replaces "TF"
    out <- out[!.layerKey(names(out)) %in% .layerKey(names(ov))]
    return(c(out, ov))
  }
  out[names(ov)] <- ov                            # replaces the entries it names and adds new ones
  out
}

#' @rdname scMoreColors
#' @export
scMoreSetColors <- function(kind = c("modality", "node", "sign", "owner"), values = NULL) {
  kind <- match.arg(kind)
  if (is.null(values)) {
    .scMoreEnv$colors[[kind]] <- NULL
  } else {
    if (is.null(names(values))) cli::cli_abort("{.arg values} must be a named vector of colours.")
    .scMoreEnv$colors[[kind]] <- .toHex(values)
  }
  invisible(scMoreColors(kind))
}

#' @rdname scMoreColors
#' @export
scMoreModalityColors <- function(levels) {
  known <- scMoreColors("modality")
  key   <- .layerKey(levels)
  out   <- stats::setNames(unname(known[match(key, .layerKey(names(known)))]), levels)
  if (anyNA(out)) {                                  # other layers: the next free colours, one per layer
    used <- unlist(list(known, scMoreColors("node")))
    pool <- .scMoreColorList[!toupper(.scMoreColorList) %in% toupper(used)]
    miss <- is.na(out)
    new  <- unique(key[miss])
    cols <- stats::setNames(rep_len(pool, length(new)), new)
    out[miss] <- unname(cols[key[miss]])
  }
  out
}

#' @rdname scMoreColors
#' @export
scMorePalette <- function(n = NULL) {
  if (is.null(n)) return(.scMoreColorList)
  if (n <= length(.scMoreColorList)) return(.scMoreColorList[seq_len(n)])
  grDevices::colorRampPalette(.scMoreColorList)(n)
}
