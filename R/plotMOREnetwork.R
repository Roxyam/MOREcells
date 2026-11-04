# =============================================================================
# scMore networks
#
# Everything about the regulation networks, in one file. A network is a plain igraph
# object: its style (colours, shapes, sizes, labels, arrows, line types) is stored as vertex
# and edge attributes and the layout in the graph attribute `layout`, so
# plot(g, layout = g$layout) works as for any igraph. The colours come from plotPalettes.R.
#
#   PART 1  networks of a group        selectEdges, buildNodes, styleNetwork, regulationNetwork
#   PART 2  networks of a contrast     classifyDifferential, pairSlopes, contrastEdges, contrastNetwork
#   PART 3  layout, plots, files       networkLayout, plotNetwork, as_visNetwork, networkSave
#
# Both kinds of network share steps 2 to 4 (nodes, style, igraph + layout); only the first
# step, which says which regulations are drawn, is different.
# =============================================================================

# Statuses of a regulation in a contrast (see classifyDifferential())

.statusLevels  <- c("only_A", "only_B", "stronger_A", "stronger_B", "reversed",
                    "conserved", "suggestive", "divergent", "none")
.changedStatus <- c("only_A", "only_B", "stronger_A", "stronger_B", "reversed")


# ===========================================================================
# PART 1: networks of a group
#
# table of results  ->  edges  ->  nodes  ->  style  ->  one igraph per group
# ===========================================================================

# ---- step 1: which regulations are drawn ------------------------------------------

#' Select the regulations to draw
#'
#' Takes the table returned by \code{\link{testRegulation}} and keeps the rows that
#' will become edges. Nothing is drawn or built here, so you can look at the result
#' and decide.
#'
#' The steps, in this order:
#' \enumerate{
#'   \item Keep the rows with a slope different from 0 (and significant, if
#'     \code{onlySignificant}).
#'   \item Keep only the groups in \code{groups}, and only the \code{regulators} /
#'     \code{targets} you name.
#'   \item \strong{Percentile of the effect}: in each group, keep the edges whose
#'     absolute slope is at or above the \code{effect_percentile} of that group
#'     (0 keeps all, 0.9 keeps the strongest 10 percent).
#'   \item Optional, to get hubs with their targets: \code{top_regulators} keeps, in
#'     each group and layer, the regulators with most targets and \code{top_targets}
#'     the strongest targets of each; \code{top} keeps the strongest edges of each
#'     group.
#' }
#'
#' @param res Output of \code{testRegulation()}: columns \code{target},
#'   \code{regulator}, \code{estimate} and, if present, \code{layer},
#'   \code{significant} and the column in \code{group_by}.
#' @param group_by Column that defines the groups. Default \code{"group"}; if it is
#'   not in \code{res}, all rows form one group called \code{"all"}. Several columns
#'   are pasted together.
#' @param groups Names of the groups to keep (default all), in the order you give.
#' @param onlySignificant Keep only significant rows (needs a \code{significant}
#'   column; a warning says so if it is missing).
#' @param effect_percentile Number in [0, 1). Percentile of the absolute slope used
#'   as cut-off within each group.
#' @param top_regulators,top_targets,top Optional limits (see above); \code{NULL}
#'   means no limit.
#' @param regulators,targets Keep only these regulators / targets (exact names).
#' @return A data.frame with one row per edge: \code{group}, \code{layer},
#'   \code{regulator}, \code{target}, \code{estimate}, \code{weight} (absolute
#'   slope) and \code{sign} (\code{"positive"} or \code{"negative"}).
#' @examples
#' \dontrun{
#' res <- testRegulation(fit)
#' e <- selectEdges(res, effect_percentile = 0.9)
#' table(e$group)
#' }
#' @family scMore networks
#' @export
selectEdges <- function(res, group_by = "group", groups = NULL, onlySignificant = TRUE,
                        effect_percentile = 0, top_regulators = NULL, top_targets = NULL, top = NULL,
                        regulators = NULL, targets = NULL) {
  if (!is.data.frame(res)) cli::cli_abort("{.arg res} must be a data.frame.")
  miss <- setdiff(c("target", "regulator", "estimate"), names(res))
  if (length(miss)) cli::cli_abort("{.arg res} needs the columns {.val {miss}} (output of {.fn testRegulation}).")
  if (!is.numeric(effect_percentile) || length(effect_percentile) != 1 || effect_percentile < 0 || effect_percentile >= 1) {
    cli::cli_abort("{.arg effect_percentile} must be a number in [0, 1).")
  }
  if (!is.null(group_by) && !all(group_by %in% names(res))) {
    if (identical(group_by, "group")) group_by <- NULL
    else cli::cli_abort("Column(s) {.val {setdiff(group_by, names(res))}} of {.arg group_by} not found in {.arg res}.")
  }
  group <- if (is.null(group_by)) rep("all", nrow(res)) else do.call(paste, c(lapply(res[group_by], as.character), sep = " / "))

  # 1. which rows are regulations
  keep <- !is.na(res$estimate) & res$estimate != 0
  if (onlySignificant) {
    if ("significant" %in% names(res)) keep <- keep & res$significant %in% TRUE
    else cli::cli_warn("{.arg res} has no {.field significant} column, so all the rows are kept.")
  }
  # 2. groups, regulators, targets
  if (!is.null(groups)) {
    bad <- setdiff(groups, unique(group))
    if (length(bad)) cli::cli_abort("Unknown {.arg groups}: {.val {bad}}. Available: {.val {unique(group)}}.")
    keep <- keep & group %in% groups
  }
  layer <- if ("layer" %in% names(res)) as.character(res$layer) else rep("regulator", nrow(res))
  layer[is.na(layer)] <- "regulator"
  e <- data.frame(group = group, layer = layer, regulator = as.character(res$regulator),
                  target = as.character(res$target), estimate = res$estimate, weight = abs(res$estimate),
                  sign = ifelse(res$estimate > 0, "positive", "negative"), stringsAsFactors = FALSE)[keep, , drop = FALSE]
  if (!is.null(regulators)) e <- e[e$regulator %in% regulators, , drop = FALSE]
  if (!is.null(targets))    e <- e[e$target %in% targets, , drop = FALSE]
  if (nrow(e) == 0) cli::cli_abort("No edges left after filtering.")

  # 3. percentile of the effect, inside each group
  if (effect_percentile > 0) {
    cut <- tapply(e$weight, e$group, stats::quantile, probs = effect_percentile, names = FALSE)
    e <- e[e$weight >= cut[e$group], , drop = FALSE]
  }
  # 4. optional hubs / global limit
  e <- .limitEdges(e, top_regulators, top_targets, top)

  e <- e[order(match(e$group, if (is.null(groups)) unique(e$group) else groups)), , drop = FALSE]
  rownames(e) <- NULL
  n_by <- table(e$group)
  if (is.null(top) && effect_percentile == 0 && is.null(top_regulators) && is.null(top_targets) && max(n_by) > 300) {
    cli::cli_inform(c("i" = "The largest group has {max(n_by)} edges, which is hard to read as a network.",
                      " " = "Try {.code effect_percentile = 0.9}, or {.arg regulators} to focus on some regulators."))
  }
  e
}

# the optional limits of selectEdges() and contrastEdges(); `e` needs the columns
# group, layer, regulator and weight
.limitEdges <- function(e, top_regulators = NULL, top_targets = NULL, top = NULL) {
  if (!is.null(top_regulators) || !is.null(top_targets)) {
    by_layer <- split(seq_len(nrow(e)), list(e$group, e$layer), drop = TRUE)
    idx <- unlist(lapply(by_layer, function(i) {
      n_targets <- sort(table(e$regulator[i]), decreasing = TRUE)
      if (!is.null(top_regulators)) n_targets <- utils::head(n_targets, top_regulators)
      unlist(lapply(names(n_targets), function(r) {
        j <- i[e$regulator[i] == r]
        if (!is.null(top_targets)) j <- j[order(-e$weight[j])][seq_len(min(top_targets, length(j)))]
        j
      }))
    }))
    e <- e[sort(idx), , drop = FALSE]
  }
  if (!is.null(top)) {
    idx <- unlist(lapply(split(seq_len(nrow(e)), e$group), function(i) i[order(-e$weight[i])][seq_len(min(top, length(i)))]))
    e <- e[sort(idx), , drop = FALSE]
  }
  e
}

# ---- step 2: which nodes exist, and what each one is -------------------------------

#' Build the nodes of the networks
#'
#' Gives every edge the id of its two end nodes and builds the node table, one row
#' per node and group.
#'
#' \strong{Ids.} A regulator is identified by its layer and name (\code{TF_STAT1},
#' \code{ATAC_peak1}), so regulators of different layers never get mixed up. A target
#' is identified by its name (\code{STAT1}).
#'
#' \strong{Merging a regulator and a target.} With
#' \code{merge_regulator_target = TRUE} a gene that is a target and also a
#' regulator \emph{in the same group} becomes one node (the regulator's), which
#' sends and receives regulations. With \code{FALSE} they are two nodes. If several
#' layers have a regulator with that name, the target joins the first one (in the
#' order of \code{merge_layers}, or of appearance). \code{merge_layers} limits which
#' layers can be merged (e.g. \code{"TF"}); by default all.
#'
#' \strong{Role.} In each group, a node is a \code{regulator} (it only regulates),
#' a \code{target} (it is only regulated) or \code{both}.
#'
#' @param edges Output of \code{\link{selectEdges}}.
#' @param merge_regulator_target Merge a regulator and a target with the same name.
#' @param merge_layers Layers whose regulators can be merged with a target;
#'   \code{NULL} for all.
#' @return A list with \code{edges} (the input plus \code{from} and \code{to}, the
#'   ids of the end nodes) and \code{nodes} (columns \code{group}, \code{id},
#'   \code{symbol} (the plain name), \code{role}, \code{layer}, \code{n_out} (number
#'   of targets) and \code{n_in} (number of regulators)).
#' @family scMore networks
#' @export
buildNodes <- function(edges, merge_regulator_target = TRUE, merge_layers = NULL) {
  edges$from <- paste(edges$layer, edges$regulator, sep = "_")
  edges$to   <- edges$target
  if (merge_regulator_target) {
    for (g in unique(edges$group)) {
      in_g <- edges$group == g
      ok   <- in_g & (if (is.null(merge_layers)) TRUE else edges$layer %in% merge_layers)
      regs <- edges[ok, c("regulator", "from", "layer"), drop = FALSE]
      if (!is.null(merge_layers)) regs <- regs[order(match(regs$layer, merge_layers)), , drop = FALSE]
      regs <- regs[!duplicated(regs$regulator), , drop = FALSE]
      joined <- regs$from[match(edges$target[in_g], regs$regulator)]   # NA: not a regulator here
      edges$to[in_g] <- ifelse(is.na(joined), edges$target[in_g], joined)
    }
  }
  nodes <- do.call(rbind, lapply(unique(edges$group), function(g) {
    d     <- edges[edges$group == g, , drop = FALSE]
    ids   <- unique(c(d$from, d$to))
    pairs <- unique(d[, c("from", "to")])
    n_out <- as.integer(table(factor(pairs$from, levels = ids)))
    n_in  <- as.integer(table(factor(pairs$to,   levels = ids)))
    symbol <- c(stats::setNames(d$regulator, d$from), stats::setNames(d$target, d$to))[ids]
    role   <- ifelse(n_out > 0 & n_in > 0, "both", ifelse(n_out > 0, "regulator", "target"))
    layer_of <- stats::setNames(d$layer, d$from)
    data.frame(group = g, id = ids, symbol = unname(symbol), role = role,
               layer = ifelse(role == "target", "gene", unname(layer_of[ids])),
               n_out = n_out, n_in = n_in, stringsAsFactors = FALSE)
  }))
  rownames(nodes) <- NULL
  list(edges = edges, nodes = nodes, group_labels = attr(edges, "group_labels"))
}

# ---- step 3: colours, shapes, sizes -------------------------------------------------

# the network files need the colours defined in plotPalettes.R; an old copy of that file gives
# a confusing error, so say it clearly
.checkPalettes <- function() {
  if (!exists(".layerKey") || !"node" %in% eval(formals(scMoreColors)$kind)) {
    cli::cli_abort(c("{.file plotPalettes.R} is an old version.",
                     "i" = "Replace it with the latest copy and reload ALL the files, in this order: plotPalettes.R, plotTheme.R, network.R."))
  }
}

#' Give the networks their style
#'
#' Adds the visual attributes to the node and edge tables, using the names that
#' igraph understands (\code{color}, \code{shape}, \code{size}, \code{label},
#' \code{width}...), so the graphs can be plotted with plain \code{plot()}.
#'
#' \itemize{
#'   \item \strong{Node colour}: by layer, with the fixed colours of
#'     \code{\link{scMoreColors}("modality")} (the first colours of the general
#'     palette: TF, ATAC...); targets have their own colour,
#'     \code{\link{scMoreColors}("node")}.
#'   \item \strong{Node shape}: in igraph, square for ATAC and circle for the rest. In
#'     visNetwork every layer and the targets get a different shape (attribute
#'     \code{vis.shape}): by default TF diamond, ATAC square, miRNA triangle, targets
#'     circle ("dot"), and other layers hexagon, star... \code{shapes} changes them
#'     (and so does the \code{shapes} argument of \code{\link{as_visNetwork}}, which
#'     needs no rebuilding).
#'   \item \strong{Node size}: regulators are larger, and grow with their number of
#'     targets; targets are small. The scale is common to all groups, and it shrinks
#'     automatically when the network has many nodes (see \code{size_scale}); with
#'     big nodes in a crowded network the edges are shorter than the nodes and the
#'     arrowheads disappear underneath them.
#'     (visNetwork does not need the reduction and uses the attribute \code{vis.size}.)
#'   \item \strong{Frame}: a node that is both regulator and target keeps its layer
#'     colour and shape, with a frame in the colour of the targets.
#'   \item \strong{Label}: the plain name, only for the \code{label_top} regulators
#'     with most targets and the \code{label_targets} targets with most regulators of
#'     each group (all the other nodes have no label, to avoid overlapping text).
#'     Both take a number or \code{"all"}.
#'   \item \strong{Edge colour}: by sign (\code{\link{scMoreColors}("sign")}); in a
#'     contrast, by the group the edge belongs to (\code{scMoreColors("owner")}), with a
#'     dashed line for "stronger in one group" and two curved edges for a reversal.
#'   \item \strong{Edge width}: proportional to the absolute slope, with a scale
#'     common to all groups.
#'   \item \strong{Arrowheads}: from regulator to target. A positive slope ends in an
#'     arrowhead and a negative one in a bar (\code{\link{plotNetwork}} draws the bars
#'     in igraph; visNetwork has them built in). The heads are smaller when there are
#'     many edges, which are also made fainter (\code{edge_alpha}).
#' }
#'
#' @param built Output of \code{\link{buildNodes}}.
#' @param directed Draw arrowheads.
#' @param label_top Regulators labelled per group (those with most targets): a number
#'   or \code{"all"}.
#' @param label_targets Targets labelled per group (those with most regulators): a
#'   number, \code{"all"}, or 0 (default) for none.
#' @param shapes Named vector that changes the visNetwork shapes: the names are the
#'   layers, and \code{gene} is the targets; e.g. \code{c(TF = "star", gene = "square")}.
#'   Valid shapes: dot (circle), square, triangle, triangleDown, diamond, hexagon, star.
#' @param size_scale Multiplier of the node sizes. \code{NULL} (default): 1 for up to
#'   40 nodes and smaller as the largest group grows (down to 0.4), so that the
#'   arrowheads stay visible. A number overrides it (1 keeps the sizes as they are).
#' @param edge_alpha Edge transparency, 0 to 1. Default: opaque, 0.55 if a group has
#'   more than 150 edges.
#' @return The same list with the attribute columns added.
#' @family scMore networks
#' @export
styleNetwork <- function(built, directed = TRUE, label_top = 20, label_targets = 0,
                         edge_alpha = NULL, size_scale = NULL, shapes = NULL) {
  .checkPalettes()
  nodes <- built$nodes; edges <- built$edges
  isreg <- nodes$role != "target"
  layer_cols <- scMoreModalityColors(unique(nodes$layer[isreg]))
  max_out <- max(nodes$n_out[isreg], 1); max_in <- max(nodes$n_in, 1)

  target_color <- scMoreColors("node")[["target"]]
  nodes$color <- ifelse(isreg, unname(layer_cols[nodes$layer]), target_color)
  nodes$shape <- ifelse(isreg & .layerKey(nodes$layer) == "ATAC", "square", "circle")
  vis_shapes  <- .visShapes(unique(nodes$layer), shapes)
  nodes$vis.shape <- unname(vis_shapes[nodes$layer])
  if (is.null(size_scale)) size_scale <- min(1, max(0.4, sqrt(40 / max(table(nodes$group)))))
  base_size   <- ifelse(isreg, 7 + 8 * sqrt(nodes$n_out / max_out), 3 + 2 * sqrt(nodes$n_in / max_in))
  nodes$size  <- size_scale * base_size
  # visNetwork zooms to fit, so it does not need the reduction. Regulators and targets get
  # similar sizes (regulators a little larger, and more so the more targets they have);
  # as_visNetwork() adjusts it to the area of each shape.
  nodes$vis.size <- ifelse(isreg, 10 + 6 * sqrt(nodes$n_out / max_out), 8 + 2 * sqrt(nodes$n_in / max_in))
  nodes$frame.color <- ifelse(nodes$role == "both", target_color, "#6B6B6B")
  nodes$frame.width <- ifelse(nodes$role == "both", 2.5, 1)
  # labels: the regulators with most targets and, if asked, the targets with most regulators
  nodes$label <- ""
  n_top <- .howMany(label_top); n_tgt <- .howMany(label_targets)
  for (g in unique(nodes$group)) {
    i <- which(nodes$group == g & isreg)
    i <- i[order(-nodes$n_out[i])][seq_len(min(n_top, length(i)))]
    nodes$label[i] <- nodes$symbol[i]
    j <- which(nodes$group == g & !isreg)
    j <- j[order(-nodes$n_in[j])][seq_len(min(n_tgt, length(j)))]
    nodes$label[j] <- nodes$symbol[j]
  }
  nodes$label.color <- "#2B2B2B"; nodes$label.cex <- 0.7; nodes$label.dist <- 1
  nodes$label.degree <- -pi / 2; nodes$label.family <- "sans"

  max_e <- max(table(edges$group))
  if (is.null(edge_alpha)) edge_alpha <- if (max_e > 150) 0.55 else 1
  sign_cols <- scMoreColors("sign")
  contrast  <- "owner" %in% names(edges) && any(!is.na(edges$owner))
  if (contrast) {
    # a contrast: the colour says which group the edge belongs to
    edges$color <- unname(scMoreColors("owner")[edges$owner])
    edges$lty   <- unname(c(solid = 1, dashed = 2, dotted = 3)[edges$linetype])
    edges$curved <- ifelse(edges$curved, ifelse(edges$side == "A", 0.25, -0.25), 0)   # the two arrows of a reversal
    edge_key <- list(kind = "contrast", owners = unique(edges$owner),
                     linetypes = unique(edges$linetype), reversed = any(!is.na(edges$side)),
                     divergent = any(edges$status == "divergent"), suggestive = any(edges$status == "suggestive"))
  } else {
    edges$color <- unname(sign_cols[edges$sign])
    edges$lty   <- 1
    edges$curved <- 0
    edge_key <- list(kind = "group")
  }
  if (edge_alpha < 1) edges$color <- grDevices::adjustcolor(edges$color, edge_alpha)
  edges$width <- 0.7 + 2.6 * edges$weight / max(edges$weight, 1e-9)
  # a positive slope ends in an arrowhead; a negative one in a bar (drawn by plotNetwork())
  edges$arrow.mode  <- if (directed) ifelse(edges$sign == "positive", 2, 0) else 0
  edges$arrow.size  <- if (max_e <= 50) 0.45 else if (max_e <= 150) 0.28 else 0.2
  edges$arrow.width <- 1.1
  list(edges = edges, nodes = nodes, layer_colors = layer_cols, sign_colors = sign_cols,
       target_color = target_color, vis_shapes = vis_shapes, edge_key = edge_key,
       group_labels = built$group_labels)
}

# how many labels: a number, "all" or NULL (= all)
.howMany <- function(x) {
  if (is.null(x) || identical(x, "all")) return(Inf)
  if (!is.numeric(x) || length(x) != 1 || x < 0) cli::cli_abort("{.arg label_top} and {.arg label_targets} must be a number, {.val all} or {.code NULL}.")
  x
}

# visNetwork shape of each layer (and of the targets, called "gene")
.visShapes <- function(layers, user = NULL) {
  valid <- c("dot", "square", "triangle", "triangleDown", "diamond", "hexagon", "star")
  known <- c(GENE = "dot", TF = "diamond", ATAC = "square", MIRNA = "triangle")   # names as .layerKey() gives them
  if (!is.null(user)) {
    if (is.null(names(user)) || !all(user %in% valid)) {
      cli::cli_abort("{.arg shapes} must be a named vector with shapes among {.val {valid}}.")
    }
    known[.layerKey(names(user))] <- user
  }
  key <- .layerKey(layers)
  out <- stats::setNames(rep(NA_character_, length(layers)), layers)
  hit <- key %in% names(known)
  out[hit] <- known[key[hit]]
  if (any(!hit)) {
    pool <- setdiff(c("hexagon", "star", "triangleDown", "triangle", "square", "diamond"), c(out[hit], "dot"))
    out[!hit] <- rep_len(pool, sum(!hit))
  }
  out
}

# ---- step 4: the igraph objects ------------------------------------------------------

#' Build regulation networks
#'
#' Runs \code{\link{selectEdges}}, \code{\link{buildNodes}} and
#' \code{\link{styleNetwork}} and returns one igraph per group, with the style as
#' vertex and edge attributes and the layout in \code{g$layout}. Read the help of
#' those three functions for what each step decides; the arguments are passed on.
#'
#' @section What you get:
#' With \code{groups = NULL} (default), a named list with one igraph per group
#' (class \code{scMoreNetworks}). With a single name, \code{groups = "Tcell"}, that
#' igraph alone; with several, the list of those. Scales (sizes, widths) and the
#' shared layout are computed on the groups you ask for.
#'
#' @inheritParams selectEdges
#' @inheritParams buildNodes
#' @inheritParams styleNetwork
#' @param layout Layout method: \code{"fr"} (default, Fruchterman-Reingold, as in
#'   MORE), \code{"stars"}, \code{"stress"}, \code{"kk"} or \code{"circle"}. See
#'   \code{\link{networkLayout}}.
#' @param shared_layout \code{FALSE} (default): one layout per group.
#'   \code{TRUE}: one layout for all, so a node has the same position in each.
#' @param seed Random seed for the layout.
#' @return A list of igraph objects or one igraph.
#' @examples
#' \dontrun{
#' res  <- testRegulation(fit)
#' nets <- regulationNetwork(res, effect_percentile = 0.9)   # one igraph per group
#' plotNetwork(nets)
#' g <- regulationNetwork(res, groups = "Tcell")
#' plot(g, layout = g$layout)                                 # a plain igraph
#' }
#' @family scMore networks
#' @export
regulationNetwork <- function(res, group_by = "group", groups = NULL, onlySignificant = TRUE,
                              effect_percentile = 0, top_regulators = NULL, top_targets = NULL, top = NULL,
                              regulators = NULL, targets = NULL,
                              merge_regulator_target = TRUE, merge_layers = NULL,
                              directed = TRUE, label_top = 20, label_targets = 0,
                              edge_alpha = NULL, size_scale = NULL, shapes = NULL,
                              layout = c("fr", "stars", "stress", "kk", "circle"),
                              shared_layout = FALSE, seed = 1) {
  layout <- match.arg(layout)
  edges <- selectEdges(res, group_by, groups, onlySignificant, effect_percentile,
                       top_regulators, top_targets, top, regulators, targets)
  .assemble(edges, groups, merge_regulator_target, merge_layers, directed, label_top, label_targets,
            edge_alpha, size_scale, shapes, layout, shared_layout, seed)
}

# nodes + style + one igraph per group + layout; shared by regulationNetwork() and contrastNetwork()
.assemble <- function(edges, groups, merge_regulator_target, merge_layers, directed, label_top, label_targets,
                      edge_alpha, size_scale, shapes, layout, shared_layout, seed) {
  st <- styleNetwork(buildNodes(edges, merge_regulator_target, merge_layers), directed, label_top, label_targets,
                     edge_alpha, size_scale, shapes)
  group_names <- unique(st$edges$group)
  graphs <- lapply(group_names, function(gp) {
    ed <- st$edges[st$edges$group == gp, , drop = FALSE]
    nd <- st$nodes[st$nodes$group == gp, , drop = FALSE]
    names(nd)[names(nd) == "id"] <- "name"                  # igraph reads the first column as the vertex name
    nd <- nd[, c("name", setdiff(names(nd), c("name", "group")))]
    ed <- ed[, c("from", "to", setdiff(names(ed), c("from", "to", "regulator", "target")))]
    g <- igraph::graph_from_data_frame(ed, directed = directed, vertices = nd)
    g$group <- gp; g$layer_colors <- st$layer_colors; g$sign_colors <- st$sign_colors
    g$target_color <- st$target_color; g$vis_shapes <- st$vis_shapes
    g$edge_key <- st$edge_key; g$group_labels <- st$group_labels
    g
  })
  names(graphs) <- group_names
  nets <- structure(graphs, class = c("scMoreNetworks", "list"))
  nets <- networkLayout(nets, method = layout, shared = shared_layout, seed = seed)
  if (!is.null(groups) && length(groups) == 1) return(nets[[groups]])
  nets
}

# ===========================================================================
# PART 2: networks of a contrast between two groups
#
# Same steps as PART 1; only the first one changes: instead of the significant slopes of each
# group, each regulation is classified by what happened to it between A and B.
#
# How an edge is drawn:
#   colour    the group it belongs to (A, B) or the "shared" colour if it is conserved
#   line      solid, or dashed if it is the same regulation but stronger in one group
#   head      an arrowhead for a positive slope, a bar for a negative one
#   reversal  TWO arrows between the same nodes, one per group, each in its own colour
# ===========================================================================

#' Classify what happened to a regulation between two groups
#'
#' Combines the slope of the regulator in two groups (A and B), whether each is
#' significant and whether the difference between them is, into one status.
#'
#' \describe{
#'   \item{\code{only_A}, \code{only_B}}{significant in one group only, and the
#'     difference is supported.}
#'   \item{\code{stronger_A}, \code{stronger_B}}{significant in both with the same
#'     sign, larger effect in one, difference supported.}
#'   \item{\code{reversed}}{significant in both with opposite signs.}
#'   \item{\code{conserved}}{significant in both with the same sign, no supported
#'     difference.}
#'   \item{\code{suggestive}}{significant in one group, difference not supported
#'     (inconclusive).}
#'   \item{\code{divergent}}{significant in neither, but the difference is.}
#'   \item{\code{none}}{no evidence.}
#' }
#'
#' @param beta_A,beta_B Slopes in each group.
#' @param sig_A,sig_B Whether each slope is significant.
#' @param sig_diff Whether the difference between the groups is significant.
#' @return A factor with the statuses above.
#' @examples
#' classifyDifferential(c(0.5, 0.4, -0.3), c(0, 0.1, 0.4),
#'                      c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE), c(TRUE, TRUE, TRUE))
#' @family scMore networks
#' @export
classifyDifferential <- function(beta_A, beta_B, sig_A, sig_B, sig_diff) {
  sig_A <- sig_A %in% TRUE; sig_B <- sig_B %in% TRUE; sig_diff <- sig_diff %in% TRUE
  bA <- ifelse(is.na(beta_A), 0, beta_A)
  bB <- ifelse(is.na(beta_B), 0, beta_B)
  both <- sig_A & sig_B
  same <- both & sign(bA) == sign(bB)
  st <- rep("none", length(bA))
  st[!sig_A & !sig_B & sig_diff]    <- "divergent"
  st[xor(sig_A, sig_B) & !sig_diff] <- "suggestive"
  st[sig_A & !sig_B & sig_diff]     <- "only_A"
  st[!sig_A & sig_B & sig_diff]     <- "only_B"
  st[same & !sig_diff]              <- "conserved"
  st[same & sig_diff & abs(bA) >= abs(bB)] <- "stronger_A"
  st[same & sig_diff & abs(bA) <  abs(bB)] <- "stronger_B"
  st[both & !same]                  <- "reversed"
  factor(st, levels = .statusLevels)
}

# ---- helper: the slope of each regulation in the two groups -------------------------------

#' Put the slopes of two groups side by side
#'
#' Builds the table that \code{\link{contrastEdges}} needs from the output of
#' \code{\link{testRegulation}}: for every regulator-target pair, its slope and
#' significance in group A and in group B, and whether the difference is significant.
#'
#' \code{res} must contain \strong{all} the rows of \code{testRegulation()} (not only the
#' significant ones): a pair missing in a group counts as slope 0, not significant.
#'
#' The significance of the difference comes from \code{contrast_res} (the output of
#' \code{\link{contrastRegulation}}). Without it, an approximate z test is used with the
#' standard errors of the two slopes, which ignores their covariance, so it is only a
#' rough guide; a warning says so.
#'
#' @param res Output of \code{testRegulation()} with the columns \code{target},
#'   \code{regulator}, \code{estimate}, \code{significant} (or \code{p.adj}) and, for
#'   the approximate test, \code{std.error}; plus the columns named in \code{compare}
#'   and \code{within}.
#' @param a,b The two values of the column \code{compare} to contrast (A and B).
#' @param compare Column that holds the two groups to contrast (e.g. \code{"condition"}).
#' @param within Column(s) that define the panels, i.e. the groups inside which A and B
#'   are compared (e.g. \code{"cell_type"}); \code{NULL} for a single panel.
#' @param contrast_res Optional output of \code{contrastRegulation()} for this contrast.
#' @param contrast If \code{contrast_res} has several contrasts, the one to use (value
#'   of its \code{contrast} column).
#' @param alpha Significance level when only \code{p.adj} is available.
#' @return A data.frame with \code{target}, \code{regulator}, \code{layer},
#'   \code{group}, \code{beta_A}, \code{beta_B}, \code{sig_A}, \code{sig_B} and
#'   \code{sig_diff}.
#' @family scMore networks
#' @export
pairSlopes <- function(res, a, b, compare = "condition", within = NULL, contrast_res = NULL,
                       contrast = NULL, alpha = 0.05) {
  need <- c("target", "regulator", "estimate", compare, within)
  miss <- setdiff(need, names(res))
  if (length(miss)) cli::cli_abort("{.arg res} needs the columns {.val {miss}}.")
  if (!all(c(a, b) %in% res[[compare]])) cli::cli_abort("{.arg a} and {.arg b} must be values of {.field {compare}}: {.val {unique(res[[compare]])}}.")
  sig <- if ("significant" %in% names(res)) res$significant %in% TRUE
  else if ("p.adj" %in% names(res)) res$p.adj < alpha & !is.na(res$p.adj)
  else cli::cli_abort("{.arg res} needs a {.field significant} or {.field p.adj} column.")
  se  <- if ("std.error" %in% names(res)) res$std.error else rep(NA_real_, nrow(res))
  layer <- if ("layer" %in% names(res)) as.character(res$layer) else rep("regulator", nrow(res))
  key <- do.call(paste, c(list(res$target, res$regulator, layer), lapply(within, function(w) res[[w]]), sep = "\r"))
  take <- function(g) { i <- which(res[[compare]] == g); data.frame(key = key[i], beta = res$estimate[i], sig = sig[i], se = se[i], stringsAsFactors = FALSE) }
  A <- take(a); B <- take(b)
  keys <- union(A$key, B$key)
  ia <- match(keys, A$key); ib <- match(keys, B$key)
  parts <- strsplit(keys, "\r", fixed = TRUE)
  out <- data.frame(target = vapply(parts, `[`, "", 1), regulator = vapply(parts, `[`, "", 2), layer = vapply(parts, `[`, "", 3),
                    stringsAsFactors = FALSE)
  out$group <- if (is.null(within)) "all" else vapply(parts, function(p) paste(p[-(1:3)], collapse = " / "), "")
  out$beta_A <- ifelse(is.na(ia), 0, A$beta[ia]); out$beta_B <- ifelse(is.na(ib), 0, B$beta[ib])
  out$sig_A  <- ifelse(is.na(ia), FALSE, A$sig[ia]); out$sig_B <- ifelse(is.na(ib), FALSE, B$sig[ib])
  se_A <- ifelse(is.na(ia), NA, A$se[ia]); se_B <- ifelse(is.na(ib), NA, B$se[ib])
  # significance of the difference
  if (!is.null(contrast_res)) {
    cr <- contrast_res
    if (!is.null(contrast)) cr <- cr[cr$contrast == contrast, , drop = FALSE]
    else if ("contrast" %in% names(cr) && length(unique(cr$contrast)) > 1) {
      cli::cli_abort("{.arg contrast_res} has several contrasts ({.val {unique(cr$contrast)}}): choose one with {.arg contrast}.")
    }
    csig <- if ("significant" %in% names(cr)) cr$significant %in% TRUE else cr$p.adj < alpha & !is.na(cr$p.adj)
    clayer <- if ("layer" %in% names(cr)) as.character(cr$layer) else rep("regulator", nrow(cr))
    ck <- do.call(paste, c(list(cr$target, cr$regulator, clayer), lapply(within, function(w) cr[[w]]), sep = "\r"))
    m <- match(keys, ck)
    out$sig_diff <- ifelse(is.na(m), FALSE, csig[m])
  } else {
    cli::cli_warn("No {.arg contrast_res}: the significance of the difference is approximated with the standard errors (covariance ignored).")
    z <- (out$beta_A - out$beta_B) / sqrt(se_A^2 + se_B^2)
    pv <- 2 * stats::pnorm(-abs(z))
    only_one <- is.na(z)                                  # the pair is missing in one group
    pv[only_one] <- ifelse(out$sig_A[only_one] | out$sig_B[only_one], 0, 1)
    out$sig_diff <- stats::p.adjust(pv, "BH") < alpha
  }
  out[order(match(out$group, unique(out$group))), , drop = FALSE]
}

# ---- the table that contrastPatterns() returns can be used directly ---------------------------

# contrastPatterns(res) puts next to each contrast the slope of the two compared levels
# (slope_1, slope_2), whether each is detected (state_1, state_2: "+", "-" or "0", from its
# confidence interval) and whether the contrast itself is significant. That is everything a
# contrast network needs, so it is translated here instead of being rebuilt with pairSlopes().
.isPatterns <- function(x) all(c("slope_1", "slope_2", "state_1", "state_2", "significant") %in% names(x))

.fromPatterns <- function(x, group_by, contrast) {
  if ("diff_1" %in% names(x)) {
    cli::cli_abort(c("Interaction contrasts are not supported in a network.",
                     "i" = "Use a simple contrast (e.g. {.code compare = \"condition\", by = \"cell_type\"})."))
  }
  if ("contrast" %in% names(x)) {
    cs <- unique(as.character(x$contrast))
    if (is.null(contrast) && length(cs) > 1) {
      cli::cli_abort("The table has several contrasts: choose one with {.arg contrast}. Available: {.val {cs}}.")
    }
    if (!is.null(contrast)) {
      if (!contrast %in% cs) cli::cli_abort("Unknown {.arg contrast}: {.val {contrast}}. Available: {.val {cs}}.")
      x <- x[as.character(x$contrast) == contrast, , drop = FALSE]
    }
  }
  a <- attr(x, "args")
  if (identical(group_by, "group") && !"group" %in% names(x)) group_by <- if (length(a$by)) a$by else NULL
  if (!is.null(group_by)) group_by <- intersect(group_by, names(x))
  if (length(group_by) == 0) group_by <- NULL
  if ("slopes_match" %in% names(x) && any(!x$slopes_match %in% TRUE)) {
    cli::cli_warn("{sum(!x$slopes_match %in% TRUE)} contrast(s) are not the plain difference of the two slopes (groups averaged); their edges use the slopes shown.")
  }
  ok <- !is.na(x$slope_1) & !is.na(x$slope_2)
  if (any(!ok)) {
    cli::cli_inform(c("i" = "{sum(!ok)} of {nrow(x)} rows have no slope in one of the two levels and cannot be drawn.",
                      " " = "Look at them with {.code x[is.na(x$slope_1) | is.na(x$slope_2), ]}."))
  }
  x <- x[ok, , drop = FALSE]
  layer_col <- intersect(c("layer", "omic"), names(x))[1]
  out <- data.frame(target = as.character(x$target), regulator = as.character(x$regulator),
                    layer = if (is.na(layer_col)) "regulator" else as.character(x[[layer_col]]),
                    beta_A = x$slope_1, beta_B = x$slope_2,
                    sig_A = x$state_1 %in% c("+", "-"), sig_B = x$state_2 %in% c("+", "-"),
                    sig_diff = x$significant %in% TRUE, stringsAsFactors = FALSE)
  if (!is.null(group_by)) for (g in group_by) out[[g]] <- x[[g]]
  list(pairs = out, group_by = group_by,
       labels = if (all(c("side_1", "side_2") %in% names(x)) && nrow(x)) c(as.character(x$side_1[1]), as.character(x$side_2[1])) else NULL)
}

# ---- what has each regulation of the table become? --------------------------------------------

#' Look at the status of every regulation of a contrast
#'
#' Classifies each row of the table with \code{\link{classifyDifferential}} and says how many
#' edges it would give, so you can see what \code{\link{contrastNetwork}} will draw and what
#' it will leave out. It accepts the same tables as \code{\link{contrastEdges}}.
#'
#' @inheritParams contrastEdges
#' @return A data.frame with one row per regulation: \code{group}, \code{layer},
#'   \code{regulator}, \code{target}, \code{beta_A}, \code{beta_B}, \code{sig_A},
#'   \code{sig_B}, \code{sig_diff}, \code{status} and \code{edges} (how many edges it gives
#'   if drawn: 2 for a reversal and for a divergent regulation with opposite signs, 0 for
#'   \code{none}, 1 otherwise).
#' @examples
#' \dontrun{
#' s <- contrastStatus(x)
#' table(s$status)
#' s[s$status == "divergent", ]
#' }
#' @family scMore networks
#' @export
contrastStatus <- function(pairs, group_by = "group", contrast = NULL) {
  if (.isPatterns(pairs)) {
    conv <- .fromPatterns(pairs, group_by, contrast)
    pairs <- conv$pairs; group_by <- conv$group_by
  }
  if (!is.null(group_by) && !all(group_by %in% names(pairs))) group_by <- NULL
  group <- if (is.null(group_by)) rep("all", nrow(pairs)) else do.call(paste, c(lapply(pairs[group_by], as.character), sep = " / "))
  st <- classifyDifferential(pairs$beta_A, pairs$beta_B, pairs$sig_A, pairs$sig_B, pairs$sig_diff)
  opp <- sign(pairs$beta_A) * sign(pairs$beta_B) < 0
  n_edges <- ifelse(st == "none", 0L, ifelse(st == "reversed" | (st == "divergent" & opp %in% TRUE), 2L, 1L))
  data.frame(group = group, layer = if ("layer" %in% names(pairs)) as.character(pairs$layer) else "regulator",
             regulator = as.character(pairs$regulator), target = as.character(pairs$target),
             beta_A = pairs$beta_A, beta_B = pairs$beta_B, sig_A = pairs$sig_A, sig_B = pairs$sig_B,
             sig_diff = pairs$sig_diff, status = st, edges = n_edges, stringsAsFactors = FALSE)
}

# tells what is not drawn and why (only statuses that are in the table but not in `show`)
.explainStatuses <- function(status, show, sig_diff) {
  tab    <- table(factor(status, levels = .statusLevels))
  drawn  <- tab[names(tab) %in% show & tab > 0]
  hidden <- tab[!names(tab) %in% show & tab > 0]
  fmt    <- function(x) paste0(names(x), ": ", as.integer(x), collapse = ", ")
  out <- character()
  if (length(hidden)) {
    out <- c(out, "i" = paste0(sum(tab), " regulations in the table: ", sum(drawn), " drawn (", if (length(drawn)) fmt(drawn) else "none", ")."),
             " " = paste0("Not drawn, because their status is not in `show`: ", fmt(hidden), "."),
             " " = paste0("Add them with show = c(..., \"", names(hidden)[1], "\")."))
  }
  if ("conserved" %in% show && tab[["conserved"]] == 0) {
    out <- c(out, "!" = paste0("No conserved regulation: it needs both slopes detected with the same sign and a contrast that is NOT significant; ",
                               sum(!sig_diff %in% TRUE), " of ", length(status), " rows have a non-significant contrast",
                               if (all(sig_diff %in% TRUE)) " (the table only has significant contrasts, so none can be conserved)." else "."))
  }
  if ("divergent" %in% show && tab[["divergent"]] == 0) {
    out <- c(out, "!" = "No divergent regulation: it needs a significant contrast with neither slope detected on its own.")
  }
  out
}

# ---- step 1 of a contrast: what happened to each regulation ---------------------------------

#' Select the regulations of a contrast
#'
#' The counterpart of \code{\link{selectEdges}} for two groups. Every regulation is
#' classified with \code{\link{classifyDifferential}} and only the statuses in
#' \code{show} are kept; each becomes one edge, except a \strong{reversal}, which
#' becomes \strong{two}: one per group, each with the sign it has in that group.
#'
#' \itemize{
#'   \item \code{only_A} / \code{only_B}: the edge belongs to that group (solid line).
#'   \item \code{stronger_A} / \code{stronger_B}: the same regulation in both, but stronger
#'     in one; the edge belongs to that group (dashed line).
#'   \item \code{reversed}: two edges, one of the colour of A and one of the colour of B.
#'   \item \code{conserved}: grey, only if you ask for it in \code{show}.
#'   \item \code{divergent}: the difference is significant but neither slope is, by itself
#'     (e.g. negative in A and positive in B, both within the noise, yet clearly different).
#'     It is drawn dotted, to say that the evidence is in the difference and not in each
#'     group: with opposite signs, as \strong{two dotted arrows} (one per group, like a
#'     reversal); with the same sign, as one dotted edge of the group with the larger
#'     slope. Only if you ask for it in \code{show}.
#'   \item \code{suggestive} (dotted): significant in one group but the difference is not;
#'     inconclusive. Only if you ask for it.
#' }
#'
#' \code{effect_percentile} works on the size of the change, \eqn{|\beta_A - \beta_B|}, inside
#' each group (panel), like in \code{selectEdges()}.
#'
#' @param pairs Either the output of \code{contrastPatterns()} (a simple contrast; it is
#'   used as it is: A is the first compared level and B the second, a slope is "detected"
#'   when its \code{state} is \code{+} or \code{-}, and the panels are the \code{by} of the
#'   contrast), or a table with \code{target}, \code{regulator}, \code{beta_A},
#'   \code{beta_B}, \code{sig_A}, \code{sig_B}, \code{sig_diff} and, if present,
#'   \code{layer} and the column in \code{group_by} (see \code{\link{pairSlopes}}).
#' @param group_by,groups,effect_percentile,top_regulators,top_targets,top,regulators,targets
#'   As in \code{\link{selectEdges}}.
#' @param group_labels Names of A and B, for the legend. Default: the compared levels when
#'   \code{pairs} comes from \code{contrastPatterns()}, otherwise \code{c("A", "B")}.
#' @param contrast If \code{pairs} has several contrasts, the one to draw.
#' @param show Statuses to draw. Default: the ones that changed (\code{only_A}, \code{only_B},
#'   \code{stronger_A}, \code{stronger_B}, \code{reversed}). A message says how many regulations
#'   of the table are left out because their status is not in \code{show}; \code{\link{contrastStatus}}
#'   shows the status of each one.
#' @return A data.frame with one row per edge: \code{group}, \code{layer},
#'   \code{regulator}, \code{target}, \code{status}, \code{owner} (\code{A}, \code{B} or
#'   \code{shared}), \code{side} (A or B for the two arrows of a reversal), \code{sign},
#'   \code{estimate} (the slope shown), \code{beta_A}, \code{beta_B}, \code{weight}
#'   (width of the edge), \code{linetype} and \code{curved}.
#' @family scMore networks
#' @export
contrastEdges <- function(pairs, group_by = "group", groups = NULL, group_labels = NULL,
                          show = .changedStatus, effect_percentile = 0,
                          top_regulators = NULL, top_targets = NULL, top = NULL,
                          regulators = NULL, targets = NULL, contrast = NULL) {
  if (!is.data.frame(pairs)) cli::cli_abort("{.arg pairs} must be a data.frame.")
  if (.isPatterns(pairs)) {                    # the output of contrastPatterns()
    conv <- .fromPatterns(pairs, group_by, contrast)
    pairs <- conv$pairs; group_by <- conv$group_by
    if (is.null(group_labels)) group_labels <- conv$labels
    if (is.null(group_by)) group_by <- "group"        # no panels: one group called "all"
  }
  if (is.null(group_labels)) group_labels <- c("A", "B")
  need <- c("target", "regulator", "beta_A", "beta_B", "sig_A", "sig_B", "sig_diff")
  miss <- setdiff(need, names(pairs))
  if (length(miss)) cli::cli_abort(c("{.arg pairs} needs the columns {.val {miss}}.", "i" = "See {.fn pairSlopes} to build it from {.fn testRegulation}."))
  bad <- setdiff(show, .statusLevels)
  if (length(bad)) cli::cli_abort("Unknown status in {.arg show}: {.val {bad}}. Valid: {.val {.statusLevels}}.")
  if (!is.numeric(effect_percentile) || length(effect_percentile) != 1 || effect_percentile < 0 || effect_percentile >= 1) {
    cli::cli_abort("{.arg effect_percentile} must be a number in [0, 1).")
  }
  if (length(group_labels) != 2) cli::cli_abort("{.arg group_labels} must have two names.")
  if (!is.null(group_by) && !all(group_by %in% names(pairs))) {
    if (identical(group_by, "group")) group_by <- NULL
    else cli::cli_abort("Column(s) {.val {setdiff(group_by, names(pairs))}} of {.arg group_by} not found in {.arg pairs}.")
  }
  group <- if (is.null(group_by)) rep("all", nrow(pairs)) else do.call(paste, c(lapply(pairs[group_by], as.character), sep = " / "))
  status <- classifyDifferential(pairs$beta_A, pairs$beta_B, pairs$sig_A, pairs$sig_B, pairs$sig_diff)
  keep <- status %in% show
  if (!is.null(groups)) {
    badg <- setdiff(groups, unique(group))
    if (length(badg)) cli::cli_abort("Unknown {.arg groups}: {.val {badg}}. Available: {.val {unique(group)}}.")
    keep <- keep & group %in% groups
  }
  in_scope <- if (is.null(groups)) rep(TRUE, length(group)) else group %in% groups
  msg <- .explainStatuses(status[in_scope], show, pairs$sig_diff[in_scope])
  if (length(msg)) cli::cli_inform(msg)
  layer <- if ("layer" %in% names(pairs)) as.character(pairs$layer) else rep("regulator", nrow(pairs))
  layer[is.na(layer)] <- "regulator"
  bA <- ifelse(is.na(pairs$beta_A), 0, pairs$beta_A); bB <- ifelse(is.na(pairs$beta_B), 0, pairs$beta_B)
  d <- data.frame(group = group, layer = layer, regulator = as.character(pairs$regulator), target = as.character(pairs$target),
                  status = as.character(status), beta_A = bA, beta_B = bB, delta = abs(bA - bB),
                  sig_A = pairs$sig_A %in% TRUE, stringsAsFactors = FALSE)[keep, , drop = FALSE]
  if (!is.null(regulators)) d <- d[d$regulator %in% regulators, , drop = FALSE]
  if (!is.null(targets))    d <- d[d$target %in% targets, , drop = FALSE]
  if (nrow(d) == 0) cli::cli_abort("No edges left after filtering.")
  if (effect_percentile > 0) {
    cut <- tapply(d$delta, d$group, stats::quantile, probs = effect_percentile, names = FALSE)
    d <- d[d$delta >= cut[d$group], , drop = FALSE]
  }
  d$weight <- d$delta
  d <- .limitEdges(d, top_regulators, top_targets, top)

  # one edge per regulation, two for a reversal
  mk <- function(i, owner, beta, linetype, weight = d$delta[i], side = NA_character_) {
    if (!length(i)) return(NULL)
    data.frame(group = d$group[i], layer = d$layer[i], regulator = d$regulator[i], target = d$target[i],
               status = d$status[i], owner = owner[i], side = side, sign = ifelse(beta[i] > 0, "positive", "negative"),
               estimate = beta[i], beta_A = d$beta_A[i], beta_B = d$beta_B[i], weight = weight, linetype = linetype,
               curved = !is.na(side), stringsAsFactors = FALSE)
  }
  s <- d$status; n <- nrow(d)
  owner <- function(x) rep(x, n)
  sug_owner <- ifelse(d$sig_A, "A", "B")
  div_all <- which(s == "divergent")
  div_opp <- div_all[sign(d$beta_A[div_all]) * sign(d$beta_B[div_all]) < 0]
  div_oth <- setdiff(div_all, div_opp)
  parts <- list(
    mk(which(s == "only_A"),     owner("A"), d$beta_A, "solid"),
    mk(which(s == "only_B"),     owner("B"), d$beta_B, "solid"),
    mk(which(s == "stronger_A"), owner("A"), d$beta_A, "dashed"),
    mk(which(s == "stronger_B"), owner("B"), d$beta_B, "dashed"),
    mk(which(s == "reversed"),   owner("A"), d$beta_A, "solid", abs(d$beta_A[which(s == "reversed")]), "A"),
    mk(which(s == "reversed"),   owner("B"), d$beta_B, "solid", abs(d$beta_B[which(s == "reversed")]), "B"),
    mk(which(s == "conserved"),  owner("shared"), d$beta_A, "solid", (abs(d$beta_A) + abs(d$beta_B))[which(s == "conserved")] / 2),
    mk(which(s == "suggestive"), sug_owner, ifelse(d$sig_A, d$beta_A, d$beta_B), "dotted", abs(ifelse(d$sig_A, d$beta_A, d$beta_B))[which(s == "suggestive")]),
    mk(div_opp, owner("A"), d$beta_A, "dotted", abs(d$beta_A[div_opp]), "A"),      # divergent, opposite signs: two dotted arrows
    mk(div_opp, owner("B"), d$beta_B, "dotted", abs(d$beta_B[div_opp]), "B"),
    mk(div_oth, ifelse(abs(d$beta_A) >= abs(d$beta_B), "A", "B"), ifelse(abs(d$beta_A) >= abs(d$beta_B), d$beta_A, d$beta_B), "dotted")
  )
  e <- do.call(rbind, parts)
  e <- e[order(match(e$group, if (is.null(groups)) unique(d$group) else groups)), , drop = FALSE]
  rownames(e) <- NULL
  attr(e, "group_labels") <- group_labels
  e
}

# ---- the networks ---------------------------------------------------------------------------

#' Build the networks of a contrast
#'
#' The counterpart of \code{\link{regulationNetwork}} for two groups: runs
#' \code{\link{contrastEdges}}, \code{\link{buildNodes}} and \code{\link{styleNetwork}}
#' and returns one igraph per panel (e.g. per cell type) with the colours of the two
#' groups, dashed lines for "stronger in one", and two arrows for a reversal. The
#' arguments are the ones of those functions and of \code{regulationNetwork()}.
#'
#' @inheritParams contrastEdges
#' @inheritParams buildNodes
#' @inheritParams styleNetwork
#' @inheritParams regulationNetwork
#' @return A list of igraph objects (class \code{scMoreNetworks}) or one igraph.
#' @examples
#' \dontrun{
#' res  <- contrastRegulation(fit, compare = "condition", by = "cell_type")   # components = TRUE
#' x    <- contrastPatterns(res)
#' nets <- contrastNetwork(x, effect_percentile = 0.8)      # A and B are the two compared levels
#' plotNetwork(nets)
#' }
#' @family scMore networks
#' @export
contrastNetwork <- function(pairs, group_by = "group", groups = NULL, group_labels = NULL,
                            show = .changedStatus, effect_percentile = 0,
                            top_regulators = NULL, top_targets = NULL, top = NULL,
                            regulators = NULL, targets = NULL, contrast = NULL,
                            merge_regulator_target = TRUE, merge_layers = NULL,
                            directed = TRUE, label_top = 20, label_targets = 0,
                            edge_alpha = NULL, size_scale = NULL, shapes = NULL,
                            layout = c("fr", "stars", "stress", "kk", "circle"),
                            shared_layout = FALSE, seed = 1) {
  layout <- match.arg(layout)
  edges <- contrastEdges(pairs, group_by, groups, group_labels, show, effect_percentile,
                         top_regulators, top_targets, top, regulators, targets, contrast)
  .assemble(edges, groups, merge_regulator_target, merge_layers, directed, label_top, label_targets,
            edge_alpha, size_scale, shapes, layout, shared_layout, seed)
}

# ===========================================================================
# PART 3: layout, plots and files
#
# igraph and visNetwork drawings, the legends and networkSave()
# ===========================================================================

#' @export
print.scMoreNetworks <- function(x, ...) {
  cat("<scMoreNetworks> a list of ", length(x), " igraph object(s)\n", sep = "")
  for (p in names(x)) {
    g <- x[[p]]
    cat(sprintf("  $%-12s %3d nodes (%d both), %3d edges\n", paste0("`", p, "`"), igraph::vcount(g),
                sum(igraph::V(g)$role == "both"), igraph::ecount(g)))
  }
  invisible(x)
}

# ---- layout ------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# centre and scale to [-1, 1], keeping the aspect ratio
.normXY <- function(xy) {
  xy  <- as.matrix(xy)
  xy  <- sweep(xy, 2, (apply(xy, 2, max) + apply(xy, 2, min)) / 2, "-")
  r   <- max(abs(xy))
  if (!is.finite(r) || r == 0) r <- 1
  xy / r
}

# Regulator-centred "stars": regulators (the nodes that regulate something, including
# the TFs that are also targets) are laid out apart from each other, and the targets
# that only they regulate are spread in a disc around them; a target shared by several
# regulators sits between them. The result is open and readable, unlike a force layout
# of hundreds of weakly connected targets.
.layoutStars <- function(ids, pairs) {
  pairs <- pairs[pairs$from != pairs$to, , drop = FALSE]
  reg   <- intersect(ids, unique(pairs$from))
  if (length(reg) == 0) return(igraph::layout_in_circle(igraph::make_empty_graph(length(ids))))
  pure  <- setdiff(ids, reg)
  pt    <- pairs[pairs$to %in% pure, , drop = FALSE]
  n_pt  <- as.integer(table(factor(pt$from, levels = reg)))
  rho   <- 0.45 + 0.19 * sqrt(n_pt)                    # radius of the star of each regulator

  # layout of the regulators: links between them + shared targets
  he <- pairs[pairs$to %in% reg, c("from", "to"), drop = FALSE]
  sh <- NULL
  if (nrow(pt)) {
    by_t <- split(pt$from, pt$to); by_t <- by_t[lengths(by_t) > 1]
    if (length(by_t)) sh <- do.call(rbind, lapply(by_t, function(h) {
      h <- utils::head(unique(h)[order(-n_pt[match(unique(h), reg)])], 8)
      cm <- t(utils::combn(h, 2)); data.frame(from = cm[, 1], to = cm[, 2], stringsAsFactors = FALSE)
    }))
  }
  he <- rbind(he, sh)
  hg <- igraph::graph_from_data_frame(if (is.null(he) || nrow(he) == 0) data.frame(from = character(), to = character()) else he,
                                      directed = FALSE, vertices = data.frame(id = reg))
  hg <- igraph::simplify(hg, edge.attr.comb = "sum")
  hub <- function(gg) {
    if (igraph::vcount(gg) < 3 || !requireNamespace("graphlayouts", quietly = TRUE)) return(igraph::layout_with_fr(gg, niter = 800))
    graphlayouts::layout_with_stress(gg)
  }
  xy <- if (length(reg) == 1) matrix(0, 1, 2) else igraph::layout_components(hg, layout = hub)
  xy <- as.matrix(xy)
  # keep the stars apart
  for (it in seq_len(250)) {
    dx <- outer(xy[, 1], xy[, 1], "-"); dy <- outer(xy[, 2], xy[, 2], "-")
    d  <- sqrt(dx^2 + dy^2); diag(d) <- Inf
    ov <- pmax(outer(rho, rho, "+") + 0.5 - d, 0)
    if (max(ov) < 1e-3) break
    d[d < 1e-6] <- 1e-6
    xy[, 1] <- xy[, 1] + rowSums(ov * dx / d) * 0.4
    xy[, 2] <- xy[, 2] + rowSums(ov * dy / d) * 0.4
  }
  pos <- matrix(NA_real_, length(ids), 2); rownames(pos) <- ids
  pos[reg, ] <- xy
  golden <- 2.399963
  if (nrow(pt)) {
    owners <- tapply(pt$from, pt$to, function(h) sort(unique(h)))
    key    <- vapply(owners, paste, "", collapse = "\r")
    for (k in unique(key)) {
      tg <- sort(names(owners)[key == k]); hs <- owners[[tg[1]]]; n <- length(tg)
      i  <- seq_len(n)
      if (length(hs) == 1) {
        r <- rho[match(hs, reg)] * (0.55 + 0.45 * sqrt((i - 0.5) / n))
        a <- i * golden + (match(hs, reg) %% 7)
        pos[tg, ] <- cbind(xy[match(hs, reg), 1] + r * cos(a), xy[match(hs, reg), 2] + r * sin(a))
      } else {                                           # shared: between its regulators
        ctr <- colMeans(xy[match(hs, reg), , drop = FALSE])
        r <- (0.2 + 0.1 * sqrt(n)) * sqrt((i - 0.5) / n); a <- i * golden
        pos[tg, ] <- cbind(ctr[1] + r * cos(a), ctr[2] + r * sin(a))
      }
    }
  }
  pos[is.na(pos[, 1]), ] <- 0
  pos
}

.layoutGraph <- function(ids, pairs, method) {
  if (method == "stars") return(.layoutStars(ids, pairs))
  g <- igraph::simplify(igraph::graph_from_data_frame(pairs, directed = FALSE, vertices = data.frame(id = ids)))
  stress <- function(gg) {
    if (igraph::vcount(gg) < 3 || !requireNamespace("graphlayouts", quietly = TRUE)) return(igraph::layout_with_fr(gg))
    graphlayouts::layout_with_stress(gg)
  }
  switch(method,
         stress = igraph::layout_components(g, layout = stress),
         fr     = igraph::layout_with_fr(g),
         kk     = igraph::layout_components(g, layout = igraph::layout_with_kk),
         circle = igraph::layout_in_circle(g))
}

#' Compute node positions
#'
#' Positions are stored in each graph as the graph attribute \code{layout} (a
#' matrix with one row per vertex, in \eqn{[-1, 1]}).
#'
#' With \code{shared = FALSE} (default) each group gets its own layout: every graph
#' is compact and shows its own hubs and modules, but nodes move between groups.
#' With \code{shared = TRUE} the layout is computed once on the edges of all
#' groups, so a node keeps its position in every group and a specific node can be
#' followed across groups; with many groups each graph looks sparse.
#'
#' @param nets A list of networks (\code{regulationNetwork()}) or one igraph.
#' @param method \code{"fr"} (default): Fruchterman-Reingold, the layout of the MORE
#'   plots. \code{"stars"}: regulators apart from each other, with the targets that
#'   only one regulator controls in a disc around it. \code{"stress"} and
#'   \code{"kk"}: other force layouts. \code{"circle"}: all nodes on a circle. With
#'   hundreds of weakly connected targets the force layouts pile them in the middle;
#'   \code{"stars"} spreads them.
#' @param shared One layout for all groups (\code{TRUE}) or one per group (\code{FALSE}).
#' @param seed Random seed.
#' @return The same object with the layout stored.
#' @family scMore networks
#' @export
networkLayout <- function(nets, method = c("fr", "stars", "stress", "kk", "circle"), shared = FALSE, seed = 1) {
  method <- match.arg(method)
  single <- inherits(nets, "igraph")
  L <- if (single) stats::setNames(list(nets), if (is.null(nets$group)) "all" else nets$group) else nets
  pairs_of <- function(g) {
    d <- unique(igraph::as_data_frame(g, what = "edges")[, c("from", "to")])
    if (!igraph::is_directed(g)) {                       # undirected edges: regulator first
      role <- stats::setNames(igraph::V(g)$role, igraph::V(g)$name)
      sw <- role[d$from] == "target" & role[d$to] != "target"
      d[sw, c("from", "to")] <- d[sw, c("to", "from")]
    }
    d
  }
  set.seed(seed)
  if (shared) {
    ids <- unique(unlist(lapply(L, function(g) igraph::V(g)$name)))
    xy  <- .normXY(.layoutGraph(ids, unique(do.call(rbind, lapply(L, pairs_of))), method))
    for (p in names(L)) L[[p]]$layout <- xy[match(igraph::V(L[[p]])$name, ids), , drop = FALSE]
  } else {
    for (p in names(L)) {
      g <- L[[p]]
      set.seed(seed)                       # same seed for every group: a group looks the same whichever others are drawn
      L[[p]]$layout <- .normXY(.layoutGraph(igraph::V(g)$name, pairs_of(g), method))
    }
  }
  if (single) L[[1]] else L
}

# ---- plotting ----------------------------------------------------------------

# igraph has no bar arrowhead: draw a short bar at the end of the negative edges. The end of
# an edge is on the border of the target; a curved edge (the two arrows of a reversal) ends
# along the tangent of igraph's spline, whose control point sits at the middle of the edge
# pushed sideways by the curvature.
.addBars <- function(g, lay) {
  if (!igraph::is_directed(g)) return(invisible())
  d <- igraph::as_data_frame(g, what = "edges")
  neg <- which(d$sign == "negative" & d$arrow.mode == 0)
  if (!length(neg)) return(invisible())
  vn <- igraph::V(g)$name; rad <- igraph::V(g)$size / 200
  i <- match(d$from[neg], vn); j <- match(d$to[neg], vn)
  keep <- i != j; neg <- neg[keep]; i <- i[keep]; j <- j[keep]
  if (!length(neg)) return(invisible())
  dx <- lay[j, 1] - lay[i, 1]; dy <- lay[j, 2] - lay[i, 2]
  len <- sqrt(dx^2 + dy^2); len[len == 0] <- 1
  ux <- dx / len; uy <- dy / len
  x1 <- lay[i, 1] + ux * rad[i]; y1 <- lay[i, 2] + uy * rad[i]
  x2 <- lay[j, 1] - ux * rad[j]; y2 <- lay[j, 2] - uy * rad[j]
  lambda <- if ("curved" %in% names(d)) d$curved[neg] else 0
  spx <- (x1 + x2) / 2 - lambda / 2 * (y2 - y1); spy <- (y1 + y2) / 2 + lambda / 2 * (x2 - x1)
  tx <- x2 - spx; ty <- y2 - spy; tl <- sqrt(tx^2 + ty^2); tl[tl == 0] <- 1
  tx <- tx / tl; ty <- ty / tl
  half <- 0.035
  bx <- x2 - tx * 0.004; by <- y2 - ty * 0.004
  graphics::segments(bx - ty * half, by + tx * half, bx + ty * half, by - tx * half,
                     col = d$color[neg], lwd = d$width[neg] * 1.1, lend = "butt")
  invisible()
}

# legend order of the layers: the usual ones first (TF, ATAC, miRNA), then the rest
.orderLayers <- function(x) x[order(match(.layerKey(x), c("TF", "ATAC", "MIRNA")), x)]

.drawLegend <- function(L) {
  v <- lapply(L, function(g) igraph::V(g))
  layers <- .orderLayers(unique(unlist(lapply(v, function(x) unique(x$layer[x$role != "target"])))))
  has_both <- any(vapply(v, function(x) any(x$role == "both"), logical(1)))
  lcol <- L[[1]]$layer_colors; scol <- L[[1]]$sign_colors
  tcol <- if (is.null(L[[1]]$target_color)) "#DADADA" else L[[1]]$target_color
  directed <- igraph::is_directed(L[[1]])
  ink <- "#2B2B2B"; edge_grey <- "#6B6B6B"
  graphics::par(mar = c(0, 0, 0, 0)); graphics::plot.new()
  cex <- 0.8
  step <- graphics::strheight("Mg", cex = cex) * 2.1
  y <- 1 - step * 0.6
  heading <- function(txt) { graphics::text(0.02, y, txt, adj = c(0, 0.5), font = 2, cex = cex + 0.05, col = ink); y <<- y - step }
  # nodes: a symbol and its meaning
  heading("Node (layer)")
  item <- function(pch, bg, col, lwd, label) {
    graphics::points(0.09, y, pch = pch, bg = bg, col = col, lwd = lwd, cex = 1.7)
    graphics::text(0.2, y, label, adj = c(0, 0.5), cex = cex, col = ink); y <<- y - step
  }
  for (l in layers) item(if (.layerKey(l) == "ATAC") 22 else 21, lcol[[l]], edge_grey, 1, l)
  item(21, tcol, edge_grey, 1, "Target gene")
  if (has_both) item(if (.layerKey(layers[1]) == "ATAC") 22 else 21, lcol[[layers[1]]], tcol, 2.5, "Also a target")
  y <- y - step * 0.5
  key <- L[[1]]$edge_key
  # little samples of edges: a line ending in an arrowhead or in a bar, optionally curved
  sample_edge <- function(col, head = "arrow", lty = 1, bend = 0) {
    x0 <- 0.02; x1 <- 0.17
    if (bend == 0) { xs <- c(x0, x1); ys <- c(y, y) }
    else { pts <- graphics::xspline(c(x0, (x0 + x1) / 2, x1), c(y, y + bend, y), shape = c(0, 1, 0), open = TRUE, draw = FALSE)
    xs <- pts$x; ys <- pts$y }
    n <- length(xs)
    if (!directed) { graphics::lines(xs, ys, col = col, lwd = 2.6, lty = lty); return(invisible()) }
    m <- max(which(sqrt((xs - xs[n])^2 + (ys - ys[n])^2) > 0.012), 1)     # a point far enough to give the direction
    if (head == "arrow") {
      graphics::lines(xs, ys, col = col, lwd = 2.6, lty = lty)                       # the shaft, with its line type
      hl <- sqrt((xs[n] - xs[m])^2 + (ys[n] - ys[m])^2)
      hx <- xs[n] - (xs[n] - xs[m]) / hl * 0.03; hy <- ys[n] - (ys[n] - ys[m]) / hl * 0.03   # only the last bit
      graphics::arrows(hx, hy, xs[n], ys[n], length = 0.09, angle = 25, col = col, lwd = 2.6, lty = 1)   # the head
    } else {
      graphics::lines(xs, ys, col = col, lwd = 2.6, lty = lty)
      tx <- xs[n] - xs[m]; ty <- ys[n] - ys[m]; tl <- sqrt(tx^2 + ty^2); tx <- tx / tl; ty <- ty / tl
      h <- step * 0.32
      graphics::segments(xs[n] - ty * h * 0.6, ys[n] + tx * h * 0.6, xs[n] + ty * h * 0.6, ys[n] - tx * h * 0.6, col = col, lwd = 2.6, lend = "butt")
    }
    invisible()
  }
  entry <- function(label, ...) { sample_edge(...); graphics::text(0.2, y, label, adj = c(0, 0.5), cex = cex, col = ink); y <<- y - step }
  grey <- "#7F7F7F"
  if (identical(key$kind, "contrast")) {
    gl <- L[[1]]$group_labels; owners <- unique(unlist(lapply(L, function(g) g$edge_key$owners)))
    lts <- unique(unlist(lapply(L, function(g) g$edge_key$linetypes)))
    rev <- any(vapply(L, function(g) isTRUE(g$edge_key$reversed), logical(1)))
    ocol <- scMoreColors("owner")
    heading("Edge belongs to")
    if ("A" %in% owners) entry(gl[1], ocol[["A"]])
    if ("B" %in% owners) entry(gl[2], ocol[["B"]])
    if ("shared" %in% owners) entry("Conserved", ocol[["shared"]])
    y <- y - step * 0.4; heading("Line")
    entry("only in that group", grey)
    if ("dashed" %in% lts) entry("stronger in that group", grey, lty = 2)
    if ("dotted" %in% lts) entry("weak evidence", grey, lty = 3)
    y <- y - step * 0.4; heading("Sign")
    entry("positive slope", grey, "arrow"); entry("negative slope", grey, "bar")
    if (rev) {
      y <- y - step * 0.4; heading("Opposite signs: one arrow per group")
      yy <- y; sample_edge(ocol[["A"]], "arrow", bend = step * 0.22); y <- yy
      sample_edge(ocol[["B"]], "bar", bend = -step * 0.22)
      graphics::text(0.2, y, "in the two groups", adj = c(0, 0.5), cex = cex, col = ink); y <- y - step
    }
  } else {
    heading("Regulation")
    entry("Positive slope", scol[["positive"]], "arrow")
    entry("Negative slope", scol[["negative"]], "bar")
  }
  invisible()
}

#' Plot regulation networks
#'
#' Draws one network, or a list of them as panels, with igraph graphics and the
#' style stored in the graphs, plus a legend and a title per group. The graphs are
#' plain igraph objects: \code{plot(g, layout = g$layout)} also works, and what you
#' change in their attributes (\code{V(g)$color}, \code{E(g)$width}...) shows up
#' here.
#'
#' @param nets One igraph or the list returned by \code{\link{regulationNetwork}}.
#' @param engine \code{"igraph"} (static, default) or \code{"visNetwork"}
#'   (interactive; see \code{\link{as_visNetwork}}).
#' @param groups Names of the groups to draw; default all.
#' @param mfrow Number of rows and columns of the panels, \code{c(rows, cols)};
#'   default up to three columns (two with visNetwork).
#' @param legend Draw the legend (with visNetwork, only in the first panel).
#' @param bars Draw a bar at the end of the negative edges (igraph only; visNetwork has
#'   them built in). Without it, a negative edge has no arrowhead.
#' @param main Titles of the panels; default the group names. \code{""} for none.
#' @param ... Passed to \code{plot.igraph} (e.g. \code{vertex.label.cex}) or, with
#'   visNetwork, to \code{\link{as_visNetwork}} (e.g. \code{physics}, \code{height}).
#' @return With igraph, the networks, invisibly. With visNetwork, a widget (one
#'   graph) or a browsable page with the widgets of all groups.
#' @family scMore networks
#' @export
plotNetwork <- function(nets, engine = c("igraph", "visNetwork"), groups = NULL, mfrow = NULL,
                        legend = TRUE, main = NULL, bars = TRUE, ...) {
  engine <- match.arg(engine)
  single <- inherits(nets, "igraph")
  if (single) nets <- stats::setNames(list(nets), if (is.null(nets$group)) "all" else nets$group)
  if (!is.null(groups)) {
    if (length(bad <- setdiff(groups, names(nets)))) cli::cli_abort("Unknown {.arg groups}: {.val {bad}}.")
    nets <- nets[groups]
  }
  if (engine == "visNetwork") {
    titles <- if (is.null(main)) names(nets) else rep_len(main, length(nets))
    ws <- lapply(seq_along(nets), function(i) {
      as_visNetwork(nets[[i]], main = titles[i], legend = legend && i == 1, ...)
    })
    if (length(ws) == 1) return(ws[[1]])
    return(.visPanels(ws, ncol = if (is.null(mfrow)) min(length(ws), 2) else mfrow[2]))
  }
  for (p in names(nets)) if (is.null(nets[[p]]$layout)) nets[[p]] <- networkLayout(nets[[p]])
  n <- length(nets)
  if (is.null(mfrow)) { nc <- min(n, 3); mfrow <- c(ceiling(n / nc), nc) }
  k <- prod(mfrow) - n                       # empty cells, drawn after the panels
  slots <- matrix(c(seq_len(n), if (k > 0) n + seq_len(k)), mfrow[1], mfrow[2], byrow = TRUE)
  mat <- if (legend) cbind(slots, n + k + 1) else slots
  op <- graphics::par(no.readonly = TRUE)
  on.exit({ graphics::layout(1); suppressWarnings(graphics::par(op)) }, add = TRUE)
  graphics::layout(mat, widths = if (legend) c(rep(1, mfrow[2]), 0.55) else rep(1, mfrow[2]))
  titles <- if (is.null(main)) names(nets) else rep_len(main, n)
  lim <- c(-1.08, 1.08)
  for (i in seq_len(n)) {
    graphics::par(mar = c(0.2, 0.2, 1.9, 0.2), bg = NA)
    g <- nets[[i]]
    graphics::plot(g, layout = g$layout, rescale = FALSE, xlim = lim, ylim = lim, asp = 1, main = "", ...)
    if (bars) .addBars(g, g$layout)
    if (nzchar(titles[i])) graphics::title(main = titles[i], adj = 0, font.main = 2, cex.main = 1, line = 0.4, col.main = "#2B2B2B")
  }
  for (i in seq_len(k)) graphics::plot.new()
  if (legend) .drawLegend(nets)
  invisible(nets)
}

#' Save network plots with a standard size and a transparent background
#'
#' Opens a png, pdf or svg device (chosen by the extension of \code{filename}),
#' calls \code{\link{plotNetwork}} and closes it. With \code{.html}, saves the
#' interactive visNetwork version as a self-contained file.
#'
#' @param nets Networks (see \code{\link{plotNetwork}}).
#' @param filename File name.
#' @param size Size preset of \code{\link{scMoreFigureSize}}; ignored if
#'   \code{width} and \code{height} (mm) are given.
#' @param width,height Size in millimetres.
#' @param dpi Resolution for png.
#' @param transparent Transparent background (\code{FALSE}: white).
#' @param ... Passed to \code{\link{plotNetwork}}.
#' @return The file name, invisibly.
#' @family scMore networks
#' @export
networkSave <- function(nets, filename, size = "double", width = NULL, height = NULL,
                        dpi = 300, transparent = TRUE, ...) {
  if (is.null(width) || is.null(height)) {
    d <- scMoreFigureSize(size); width <- d[["width"]]; height <- d[["height"]]
  }
  ext <- tolower(tools::file_ext(filename))
  bg  <- if (transparent) "transparent" else "white"
  if (ext == "html") {
    if (!requireNamespace("htmlwidgets", quietly = TRUE)) cli::cli_abort("Saving html needs the {.pkg htmlwidgets} package.")
    .saveSelfContained(plotNetwork(nets, engine = "visNetwork", ...), filename)
    return(invisible(filename))
  }
  switch(ext,
         png = grDevices::png(filename, width = width, height = height, units = "mm", res = dpi, bg = bg,
                              type = if (capabilities("cairo")) "cairo" else "Xlib"),
         pdf = grDevices::pdf(filename, width = width / 25.4, height = height / 25.4, bg = bg),
         svg = grDevices::svg(filename, width = width / 25.4, height = height / 25.4, bg = bg),
         cli::cli_abort("Unsupported format {.val {ext}}: use png, pdf, svg or html."))
  on.exit(grDevices::dev.off(), add = TRUE)
  plotNetwork(nets, ...)
  invisible(filename)
}

# ---- visNetwork ---------------------------------------------------------------

# "#RRGGBBAA" -> list(colour "#RRGGBB", opacity 0-1)
.splitAlpha <- function(col) {
  op <- ifelse(nchar(col) == 9, strtoi(substr(col, 8, 9), 16L) / 255, 1)
  list(color = substr(col, 1, 7), opacity = round(op, 2))
}

#' Convert one network to an interactive visNetwork widget
#'
#' Reads the style stored in the igraph (colours, shapes, sizes, labels, widths,
#' signs) and the layout, so the widget looks like the static plot: nodes coloured
#' by layer (squares for ATAC), thick border for nodes that are regulator and
#' target, red/blue edges with an arrowhead for a positive slope and a bar for a
#' negative one. Hovering a node or edge gives its details; clicking a node
#' highlights its neighbours; the box on top selects a node by name.
#'
#' The igraph layout gives the starting positions. With physics on (default) the
#' nodes settle with a force simulation (ForceAtlas2) that keeps the stars open and
#' avoids overlaps, and they react when you drag one. \code{physics = "freeze"}
#' settles them and then stops the simulation; \code{FALSE} keeps the igraph
#' positions fixed.
#'
#' @param g One igraph from \code{\link{regulationNetwork}}.
#' @param main Title.
#' @param legend Add the legend.
#' @param physics \code{TRUE} (live force simulation), \code{"freeze"} (settle, then
#'   stop) or \code{FALSE} (fixed positions).
#' @param height,width Size of the widget (CSS).
#' @param shapes Change the shapes at drawing time, without rebuilding the networks:
#'   a named vector whose names are layers (\code{gene} for the targets), e.g.
#'   \code{shapes = c(TF = "star", gene = "square")}. Valid shapes: dot (circle),
#'   square, triangle, triangleDown, diamond, hexagon, star. The size is adjusted so
#'   that every shape covers a similar area.
#' @param ... Not used.
#' @return A \code{visNetwork} widget; add steps to it as usual
#'   (\code{visNetwork::visOptions(...)}).
#' @family scMore networks
#' @export
as_visNetwork <- function(g, main = NULL, legend = TRUE, physics = TRUE,
                          height = "560px", width = "100%", shapes = NULL, ...) {
  if (!requireNamespace("visNetwork", quietly = TRUE)) cli::cli_abort("This needs the {.pkg visNetwork} package.")
  if (is.null(g$layout)) g <- networkLayout(g)
  vd <- igraph::as_data_frame(g, what = "vertices")
  ed <- igraph::as_data_frame(g, what = "edges")
  directed <- igraph::is_directed(g)
  contrast <- identical(g$edge_key$kind, "contrast")
  gl <- if (is.null(g$group_labels)) c("A", "B") else g$group_labels
  lab   <- function(id) vd$symbol[match(id, vd$name)]
  din <- igraph::degree(g, mode = "in"); dout <- igraph::degree(g, mode = "out")
  lay <- g$layout
  # shapes: the stored ones, with the changes asked for here
  vs <- if (!is.null(g$vis_shapes)) g$vis_shapes else stats::setNames(unname(vd$vis.shape[match(unique(vd$layer), vd$layer)]), unique(vd$layer))
  if (!is.null(shapes)) vs <- .visShapes(names(vs), c(vs, shapes))
  node_shape <- unname(vs[vd$layer])
  nominal <- if ("vis.size" %in% names(vd)) vd$vis.size else 2 + 1.4 * vd$size
  node_size  <- nominal * unname(.visAreaFactor[node_shape])
  nodes <- data.frame(
    id = vd$name, label = vd$label,
    title = paste0("<b>", vd$symbol, "</b><br>layer: ", vd$layer, "<br>role: ", vd$role,
                   "<br>targets: ", dout[vd$name], " | regulators: ", din[vd$name]),
    shape = node_shape, size = node_size,
    x = lay[, 1] * 450, y = -lay[, 2] * 450,
    color.background = vd$color, color.border = vd$frame.color, borderWidth = vd$frame.width,
    font.size = 15, font.color = "#2B2B2B", font.strokeWidth = 4, font.strokeColor = "#FFFFFF",
    stringsAsFactors = FALSE)
  col <- .splitAlpha(ed$color)
  slope_txt <- paste0("slope: ", signif(ed$estimate, 3))
  if (contrast) slope_txt <- ifelse(!is.na(ed$side), paste0("slope in ", ifelse(ed$side == "A", gl[1], gl[2]), ": ", signif(ed$estimate, 3)), slope_txt)
  tip <- paste0(lab(ed$from), " &rarr; ", lab(ed$to), "<br>", slope_txt, "<br>", ifelse(ed$sign == "positive", "positive", "negative"))
  if (contrast) {
    nice <- c(only_A = "only in %s", only_B = "only in %s", stronger_A = "stronger in %s", stronger_B = "stronger in %s",
              reversed = "reversed", conserved = "conserved", suggestive = "inconclusive (the difference is not significant)",
              divergent = "only the difference is significant")
    who  <- ifelse(ed$owner == "A", gl[1], ifelse(ed$owner == "B", gl[2], ""))
    st   <- ifelse(grepl("%s", nice[ed$status], fixed = TRUE), sprintf(nice[ed$status], who), nice[ed$status])
    tip  <- paste0(tip, "<br>", st, "<br>", gl[1], ": ", signif(ed$beta_A, 3), " | ", gl[2], ": ", signif(ed$beta_B, 3))
  }
  edges <- data.frame(
    from = ed$from, to = ed$to, title = tip,
    color.color = col$color, color.opacity = col$opacity, width = ed$width,
    stringsAsFactors = FALSE)
  if ("lty" %in% names(ed)) {      # dashed = stronger in one group, dotted = inconclusive
    edges$dashes <- I(lapply(ed$lty, function(l) if (l == 2) c(9, 6) else if (l == 3) c(2, 5) else FALSE))
  }
  if ("curved" %in% names(ed) && any(ed$curved != 0)) {   # the two arrows of a reversal bend to opposite sides
    edges$smooth <- I(lapply(ed$curved, function(cv) if (cv == 0) list(enabled = FALSE)
                             else list(enabled = TRUE, type = if (cv > 0) "curvedCW" else "curvedCCW", roundness = 0.3)))
  }
  if (directed) {
    edges$arrows <- I(lapply(ed$sign, function(s)
      list(to = list(enabled = TRUE, type = if (s == "positive") "arrow" else "bar", scaleFactor = 0.55))))
  }
  w <- visNetwork::visNetwork(nodes, edges, width = width, height = height,
                              main = if (!is.null(main) && nzchar(main)) list(text = main,
                                                                              style = "font-family:sans-serif;font-weight:bold;font-size:16px;text-align:left;color:#2B2B2B;") else NULL)
  w <- visNetwork::visEdges(w, smooth = FALSE, selectionWidth = 2, hoverWidth = 1.4)
  mode <- if (isFALSE(physics)) "off" else if (identical(physics, "freeze")) "freeze" else "on"
  if (mode == "off") {
    w <- visNetwork::visPhysics(w, enabled = FALSE)
  } else {
    w <- visNetwork::visPhysics(w, enabled = TRUE, solver = "forceAtlas2Based",
                                forceAtlas2Based = list(gravitationalConstant = -45, centralGravity = 0.012, springLength = 40,
                                                        springConstant = 0.08, damping = 0.75, avoidOverlap = 0.5),
                                stabilization = list(enabled = TRUE, iterations = 200, updateInterval = 50), minVelocity = 0.8)
    if (mode == "freeze") w <- visNetwork::visEvents(w, stabilizationIterationsDone = "function () { this.setOptions({ physics: false }); }")
  }
  w <- visNetwork::visInteraction(w, hover = TRUE, tooltipDelay = 60, navigationButtons = TRUE,
                                  hideEdgesOnDrag = nrow(ed) > 600)
  w <- visNetwork::visOptions(w, highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE),
                              nodesIdSelection = list(enabled = TRUE, useLabels = FALSE, main = "Select a node"))
  # the legend is HTML + SVG under the graph: visNetwork's own legend is a second network whose
  # nodes drift with the physics and that gets clipped when it has many rows
  if (legend) {
    w <- htmlwidgets::appendContent(w, .visLegendHtml(g, vd, ed, vs, contrast, gl, directed))
  }
  w
}

# ---- legend of the visNetwork widget, as HTML + SVG -------------------------------------

.svgShape <- function(shape, fill, stroke = "#6B6B6B", sw = 1, s = 18) {
  c0 <- s / 2; r <- s / 2 - sw - 0.5
  pts <- function(x, y) paste(sprintf("%.1f,%.1f", x, y), collapse = " ")
  poly <- function(x, y) sprintf('<polygon points="%s" fill="%s" stroke="%s" stroke-width="%g" stroke-linejoin="round"/>', pts(x, y), fill, stroke, sw)
  inner <- switch(shape,
                  dot = sprintf('<circle cx="%g" cy="%g" r="%g" fill="%s" stroke="%s" stroke-width="%g"/>', c0, c0, r, fill, stroke, sw),
                  square = sprintf('<rect x="%g" y="%g" width="%g" height="%g" fill="%s" stroke="%s" stroke-width="%g"/>', c0 - r * 0.9, c0 - r * 0.9, r * 1.8, r * 1.8, fill, stroke, sw),
                  diamond = poly(c(c0, c0 + r, c0, c0 - r), c(c0 - r, c0, c0 + r, c0)),
                  triangle = poly(c(c0, c0 + r, c0 - r), c(c0 - r, c0 + r * 0.8, c0 + r * 0.8)),
                  triangleDown = poly(c(c0, c0 + r, c0 - r), c(c0 + r, c0 - r * 0.8, c0 - r * 0.8)),
                  hexagon = poly(c0 + r * cos(seq(0, 5) * pi / 3), c0 + r * sin(seq(0, 5) * pi / 3)),
                  star = { a <- seq(0, 9) * pi / 5 - pi / 2; rr <- rep(c(r, r * 0.45), 5); poly(c0 + rr * cos(a), c0 + rr * sin(a)) },
                  sprintf('<circle cx="%g" cy="%g" r="%g" fill="%s" stroke="%s" stroke-width="%g"/>', c0, c0, r, fill, stroke, sw))
  sprintf('<svg width="%d" height="%d" viewBox="0 0 %d %d" style="vertical-align:middle">%s</svg>', s, s, s, s, inner)
}

# a sample edge: a line (straight, dashed or curved) ending in an arrowhead or a bar
.svgEdge <- function(col, head = "arrow", dash = NULL, curve = 0, w = 46, h = 22) {
  y <- h / 2; x0 <- 2; x1 <- w - 7; cy <- y + curve
  d <- if (curve == 0) sprintf("M%g %g L%g %g", x0, y, x1, y) else sprintf("M%g %g Q%g %g %g %g", x0, y, (x0 + x1) / 2, cy, x1, y)
  da <- if (is.null(dash)) "" else sprintf(' stroke-dasharray="%s"', dash)
  ang <- atan2(y - cy, x1 - (x0 + x1) / 2) * 180 / pi
  head_svg <- if (head == "arrow") sprintf('<polygon points="%g,%g %g,%g %g,%g" fill="%s" transform="rotate(%g %g %g)"/>', x1 + 5, y, x1 - 3, y - 4.5, x1 - 3, y + 4.5, col, ang, x1, y)
  else sprintf('<line x1="%g" y1="%g" x2="%g" y2="%g" stroke="%s" stroke-width="2.6" transform="rotate(%g %g %g)"/>', x1 + 1, y - 6, x1 + 1, y + 6, col, ang, x1, y)
  sprintf('<svg width="%d" height="%d" viewBox="0 0 %d %d" style="vertical-align:middle"><path d="%s" fill="none" stroke="%s" stroke-width="2.6" stroke-linecap="round"%s/>%s</svg>', w, h, w, h, d, col, da, head_svg)
}

.visLegendHtml <- function(g, vd, ed, vs, contrast, gl, directed) {
  esc <- htmltools::htmlEscape
  layers <- .orderLayers(unique(vd$layer[vd$role != "target"]))
  lcol <- g$layer_colors; scol <- g$sign_colors
  tcol <- if (is.null(g$target_color)) "#DADADA" else g$target_color
  item <- function(svg, label) sprintf('<span style="display:inline-flex;align-items:center;gap:6px;margin:2px 14px 2px 0;white-space:nowrap">%s<span>%s</span></span>', svg, esc(label))
  block <- function(title, items) sprintf('<div style="margin:4px 26px 4px 0"><div style="font-weight:600;margin-bottom:2px">%s</div><div>%s</div></div>', esc(title), paste(items, collapse = ""))
  node_items <- c(vapply(layers, function(l) item(.svgShape(vs[[l]], lcol[[l]]), l), ""), item(.svgShape(vs[["gene"]], tcol), "Target gene"))
  if (any(vd$role == "both")) node_items <- c(node_items, item(.svgShape(vs[[layers[1]]], lcol[[layers[1]]], stroke = tcol, sw = 2, s = 20), "Also a target"))
  blocks <- list(block("Node (layer)", node_items))
  grey <- "#7F7F7F"
  if (contrast) {
    ocol <- scMoreColors("owner"); owners <- unique(ed$owner)
    own <- c(if ("A" %in% owners) item(.svgEdge(ocol[["A"]]), gl[1]), if ("B" %in% owners) item(.svgEdge(ocol[["B"]]), gl[2]),
             if ("shared" %in% owners) item(.svgEdge(ocol[["shared"]]), "Conserved"))
    blocks <- c(blocks, list(block("Edge belongs to", own)))
    ln <- c(item(.svgEdge(grey), "only in that group"), if (any(ed$lty == 2)) item(.svgEdge(grey, dash = "7 5"), "stronger in that group"),
            if (any(ed$lty == 3)) item(.svgEdge(grey, dash = "2 4"), "weak evidence"))
    blocks <- c(blocks, list(block("Line", ln)))
    blocks <- c(blocks, list(block("Sign", c(item(.svgEdge(grey, "arrow"), "positive slope"), item(.svgEdge(grey, "bar"), "negative slope")))))
    if (any(!is.na(ed$side))) {
      both <- sprintf('<span style="display:inline-block;position:relative;width:46px;height:30px;vertical-align:middle"><span style="position:absolute;left:0;top:-2px">%s</span><span style="position:absolute;left:0;top:10px">%s</span></span>',
                      .svgEdge(ocol[["A"]], "arrow", curve = -7), .svgEdge(ocol[["B"]], "bar", curve = 7))
      blocks <- c(blocks, list(block("Opposite signs: one arrow per group", item(both, "in the two groups"))))
    }
  } else {
    blocks <- c(blocks, list(block("Regulation", c(item(.svgEdge(scol[["positive"]], "arrow"), "positive slope"),
                                                   item(.svgEdge(scol[["negative"]], "bar"), "negative slope")))))
  }
  htmltools::HTML(sprintf('<div style="display:flex;flex-wrap:wrap;font:12px/1.3 sans-serif;color:#2B2B2B;padding:6px 4px 2px 4px">%s</div>', paste(unlist(blocks), collapse = "")))
}

# size factor that gives every visNetwork shape the area of a circle of the same nominal size
.visAreaFactor <- c(dot = 1, diamond = 1.25, square = 0.886, triangle = 1.555, triangleDown = 1.555, hexagon = 1.1, star = 1.3)

# several widgets in a grid
.visPanels <- function(ws, ncol = 2) {
  cells <- lapply(ws, function(w) htmltools::tags$div(style = "min-width:0;", w))
  grid <- htmltools::tags$div(
    style = sprintf("display:grid;grid-template-columns:repeat(%d,minmax(0,1fr));gap:14px;background:#fff;padding:10px;border-radius:8px;", ncol),
    cells)
  htmltools::browsable(grid)
}

# save a widget or a page of widgets as ONE html file: the scripts and styles of the
# htmlwidgets dependencies are embedded, so it opens anywhere (no pandoc needed)
.saveSelfContained <- function(x, file, title = "scMore network") {
  rt   <- htmltools::renderTags(htmltools::tags$div(style = "background:#fff;padding:10px;", x))
  deps <- htmltools::resolveDependencies(rt$dependencies)
  tmp  <- tempfile("scMore_deps_"); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  rd <- function(base, f) paste(readLines(file.path(base, f), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  inl <- vapply(deps, function(d) {
    d <- htmltools::copyDependencyToDir(d, tmp, mustWork = TRUE)
    base <- d$src$file
    paste(c(vapply(d$stylesheet, function(f) paste0("<style>\n", rd(base, f), "\n</style>"), ""),
            vapply(d$script, function(f) paste0("<script>\n", rd(base, f), "\n</script>"), "")), collapse = "\n")
  }, character(1))
  page <- paste0("<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n",
                 "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n",
                 "<title>", htmltools::htmlEscape(title), "</title>\n",
                 "<style>body{margin:0;font-family:sans-serif;}</style>\n",
                 paste(rt$head, collapse = "\n"), "\n", paste(inl, collapse = "\n"), "\n</head>\n<body>\n",
                 rt$html, "\n</body>\n</html>\n")
  writeLines(enc2utf8(page), file, useBytes = TRUE)
  invisible(file)
}
