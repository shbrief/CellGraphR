#' Build a spatial igraph from a nuclei list
#'
#' Constructs one of five spatial graph types (Delaunay, k-NN, MST, Gabriel, or
#' RNG) from a list of nucleus objects and returns an \code{igraph} with
#' \code{x}, \code{y}, and \code{type} vertex attributes.  The resulting object
#' can be passed directly to \code{\link{plotGraph}}.
#'
#' @param nuclei_list A list of nucleus objects. Each element must contain:
#'   \describe{
#'     \item{centroid}{Numeric vector of length 2 (x, y coordinates).}
#'     \item{type}{Character scalar — cell-type label.}
#'   }
#' @param type Character(1). Graph type to build. One of \code{"delaunay"},
#'   \code{"knn"}, \code{"mst"}, \code{"gabriel"}, or \code{"rng"}.
#' @param k Integer(1). Number of nearest neighbours for \code{type = "knn"}.
#'   Ignored for all other types. Default \code{5L}.
#'
#' @return An undirected \code{igraph} with vertex attributes \code{x},
#'   \code{y} (centroid coordinates) and \code{type} (cell-type label).
#'
#' @examples
#' set.seed(42)
#' nuclei <- lapply(seq_len(30), function(i)
#'     list(centroid = runif(2, 0, 500),
#'          type     = sample(c("neoplastic", "stromal", "immune"), 1)))
#'
#' g_del <- buildSpatialGraph(nuclei, type = "delaunay")
#' g_knn <- buildSpatialGraph(nuclei, type = "knn", k = 5)
#' g_mst <- buildSpatialGraph(nuclei, type = "mst")
#'
#' @importFrom igraph graph_from_data_frame
#' @export
buildSpatialGraph <- function(nuclei_list,
                               type = c("delaunay", "knn", "mst", "gabriel", "rng"),
                               k    = 5L) {

    type      <- match.arg(type)
    centroids <- t(vapply(nuclei_list, function(n) n$centroid, numeric(2L)))
    types     <- vapply(nuclei_list, function(n) n$type, character(1L))
    n         <- nrow(centroids)

    if (n < 3L)
        stop("nuclei_list must contain at least 3 nuclei to build a spatial graph.")

    # Build edge data.frame
    if (type == "mst") {
        mst   <- .build_mst_edges(centroids)
        edges <- data.frame(from = mst$edges[, 1L], to = mst$edges[, 2L])
    } else {
        adj <- switch(type,
            delaunay = .build_delaunay_adj(centroids),
            knn      = .build_knn_adj(centroids, k),
            gabriel  = .build_gabriel_adj(centroids),
            rng      = .build_rng_adj(centroids)
        )
        idx   <- which(adj & upper.tri(adj), arr.ind = TRUE)
        edges <- data.frame(from = idx[, 1L], to = idx[, 2L])
    }

    verts <- data.frame(
        name = as.character(seq_len(n)),
        x    = centroids[, 1L],
        y    = centroids[, 2L],
        type = types,
        stringsAsFactors = FALSE
    )

    igraph::graph_from_data_frame(edges, directed = FALSE, vertices = verts)
}


#' Visualize a cell spatial graph
#'
#' Renders an \code{igraph} object as a \pkg{ggplot2} figure.  Nodes are
#' positioned at their spatial coordinates (auto-detected from vertex
#' attributes) or, if no coordinate columns are found, using a
#' Fruchterman-Reingold layout.  Nodes can be coloured by any vertex attribute.
#'
#' @param g An undirected \code{igraph} object, typically produced by
#'   \code{\link{importGraphH5}} or \code{\link{buildSpatialGraph}}.
#' @param x_col Character(1) or \code{NULL}. Vertex attribute name for x
#'   coordinates.  When \code{NULL} (default) the function probes common names
#'   \code{"x"}, \code{"cx"}, \code{"centroid_x"}, \code{"X"} in that order
#'   and falls back to a force-directed layout if none are found.
#' @param y_col Character(1) or \code{NULL}. Same as \code{x_col} for y
#'   coordinates (\code{"y"}, \code{"cy"}, \code{"centroid_y"}, \code{"Y"}).
#' @param color_by Character(1) or \code{NULL}. Vertex attribute used to colour
#'   nodes.  When \code{NULL} (default) the function probes
#'   \code{"cell_type"}, \code{"type"}, \code{"label"}.
#' @param title Character(1) or \code{NULL}. Plot title. Defaults to the
#'   \code{patient_id} graph attribute when available.
#' @param node_size Numeric(1). Point size passed to \code{geom_point}.
#'   Default \code{1}.
#' @param edge_alpha Numeric(1) in \eqn{[0,1]}. Transparency of edges.
#'   Default \code{0.3}.
#' @param image Character(1) file path to a source image (PNG, JPEG, SVS, NDPI,
#'   TIFF) or a pre-loaded raster array (e.g. from \code{png::readPNG}).
#'   When provided the image is rendered as a background layer beneath the
#'   graph.  Default \code{NULL} (no image).
#' @param image_scale Numeric(1). Scale factor applied to vertex coordinates
#'   before overlaying on the image.  Use this when the graph coordinates and
#'   image pixels are at different resolutions (e.g. \code{0.5} if the image is
#'   half the resolution of the coordinate space).  Default \code{1} (no
#'   scaling).
#' @param wsi_dim Numeric(2) or \code{NULL}. Width and height of the original
#'   whole-slide image in pixels, \code{c(width, height)}.  When provided, the
#'   exact scale factors \code{image_width / wsi_width} and
#'   \code{image_height / wsi_height} are used to map graph coordinates
#'   (in WSI pixel space) to the thumbnail image.  This gives pixel-accurate
#'   alignment and is recommended over \code{image_scale} when the WSI
#'   dimensions are known.  Default \code{NULL} (auto-detect).
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' set.seed(42)
#' nuclei <- lapply(seq_len(40), function(i)
#'     list(centroid = runif(2, 0, 500),
#'          type     = sample(c("neoplastic", "stromal", "immune"), 1)))
#'
#' g <- buildSpatialGraph(nuclei, type = "delaunay")
#' plotGraph(g)
#' plotGraph(g, color_by = "type", title = "Delaunay graph")
#'
#' @importFrom igraph vertex_attr_names vertex_attr graph_attr as_data_frame
#'   layout_with_fr V
#' @importFrom ggplot2 ggplot aes geom_segment geom_point scale_colour_discrete
#'   theme_void labs theme element_text annotation_raster scale_y_reverse
#'   coord_cartesian
#' @export
plotGraph <- function(g,
                      x_col       = NULL,
                      y_col       = NULL,
                      color_by    = NULL,
                      title       = NULL,
                      node_size   = 1,
                      edge_alpha  = 0.3,
                      image       = NULL,
                      image_scale = 1,
                      wsi_dim     = NULL) {

    # ── 1. Resolve spatial layout ─────────────────────────────────────────────
    coords <- .resolve_coord_attrs(g, x_col, y_col, required = FALSE)
    x_col  <- coords$x_col
    y_col  <- coords$y_col

    use_spatial <- !is.null(x_col) && !is.null(y_col)

    if (use_spatial) {
        coords <- cbind(
            as.numeric(igraph::vertex_attr(g, x_col)),
            as.numeric(igraph::vertex_attr(g, y_col))
        )
    } else {
        message("No spatial coordinate columns found; using Fruchterman-Reingold layout.")
        coords <- igraph::layout_with_fr(g)
    }

    # ── 1b. Load image and apply scale ──────────────────────────────────────
    has_image <- !is.null(image)
    img_raster <- NULL
    if (has_image) {
        img_raster <- .load_image_raster(image)
        img_h <- nrow(img_raster)
        img_w <- ncol(img_raster)

        if (!is.null(wsi_dim)) {
            # Exact scaling from known WSI dimensions
            coords[, 1L] <- coords[, 1L] * (img_w / wsi_dim[1L])
            coords[, 2L] <- coords[, 2L] * (img_h / wsi_dim[2L])
        } else if (image_scale != 1) {
            coords[, 1L] <- coords[, 1L] * image_scale
            coords[, 2L] <- coords[, 2L] * image_scale
        } else {
            # Auto-scale: map graph coordinate extent to image pixels.
            max_x <- max(coords[, 1L], na.rm = TRUE)
            max_y <- max(coords[, 2L], na.rm = TRUE)
            if (max_x > img_w || max_y > img_h) {
                auto_scale <- min(img_w / max_x, img_h / max_y)
                coords[, 1L] <- coords[, 1L] * auto_scale
                coords[, 2L] <- coords[, 2L] * auto_scale
                message(sprintf(
                    "Auto-scaled coords (factor: %.4f). ",
                    auto_scale),
                    "Set 'wsi_dim' for exact alignment.")
            }
        }
    }

    # ── 2. Node data ──────────────────────────────────────────────────────────
    nodes <- data.frame(x = coords[, 1L], y = coords[, 2L],
                        stringsAsFactors = FALSE)

    if (is.null(color_by))
        color_by <- .resolve_celltype_attr(g, required = FALSE)

    has_color <- !is.null(color_by) &&
        color_by %in% igraph::vertex_attr_names(g)
    if (has_color)
        nodes$color <- as.character(igraph::vertex_attr(g, color_by))

    # ── 3. Edge data ──────────────────────────────────────────────────────────
    el      <- igraph::as_data_frame(g, what = "edges")
    vnames  <- igraph::V(g)$name
    if (!is.null(vnames)) {
        idx_from <- match(el$from, vnames)
        idx_to   <- match(el$to,   vnames)
    } else {
        idx_from <- as.integer(el$from)
        idx_to   <- as.integer(el$to)
    }
    edges <- data.frame(
        x    = coords[idx_from, 1L],
        y    = coords[idx_from, 2L],
        xend = coords[idx_to,   1L],
        yend = coords[idx_to,   2L]
    )

    # ── 4. Resolve title ──────────────────────────────────────────────────────
    if (is.null(title)) {
        pid   <- igraph::graph_attr(g, "patient_id")
        title <- if (!is.null(pid) && !is.na(pid)) pid else "Cell Graph"
    }

    # ── 5. Assemble ggplot ────────────────────────────────────────────────────
    p <- ggplot2::ggplot()

    if (has_image) {
        # Image raster: origin at top-left, y increases downward.
        # annotation_raster maps ymin→bottom, ymax→top, so we pass
        # ymin=img_h and ymax=0 to flip the image correctly, then
        # use scale_y_reverse so y=0 is at the top of the plot.
        p <- p +
            ggplot2::annotation_raster(
                img_raster,
                xmin = 0, xmax = img_w,
                ymin = -img_h, ymax = 0
            )
        # Negate y coords so row 0 (image top) = plot top
        nodes$y    <- -nodes$y
        edges$y    <- -edges$y
        edges$yend <- -edges$yend
    }

    p <- p +
        ggplot2::geom_segment(
            data = edges,
            ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
            colour = if (has_image) "white" else "grey60",
            alpha  = edge_alpha,
            linewidth = 0.3
        )

    if (has_color) {
        p <- p +
            ggplot2::geom_point(
                data = nodes,
                ggplot2::aes(x = x, y = y, colour = color),
                size = node_size
            ) +
            ggplot2::scale_colour_discrete(name = color_by)
    } else {
        p <- p +
            ggplot2::geom_point(
                data = nodes,
                ggplot2::aes(x = x, y = y),
                size = node_size, colour = if (has_image) "yellow" else "steelblue"
            )
    }

    p <- p +
        ggplot2::theme_void() +
        ggplot2::labs(title = title) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

    if (has_image)
        p <- p + ggplot2::coord_cartesian(
            xlim = c(0, img_w), ylim = c(-img_h, 0), expand = FALSE
        )

    p
}
