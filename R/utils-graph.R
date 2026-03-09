# Internal helpers for igraph vertex-attribute resolution
# These functions provide DRY auto-detection of coordinate and cell-type
# attributes across all functions that accept igraph objects.

#' Resolve a vertex attribute name from candidate names
#'
#' @param g An \code{igraph} object.
#' @param candidates Character vector of candidate attribute names, checked
#'   in order.
#' @param label Character(1). Human-readable label for error messages.
#' @param required Logical(1). If \code{TRUE}, \code{stop()} when no
#'   candidate is found.  If \code{FALSE}, return \code{NULL}.
#' @return Character(1) — the first matching attribute name, or \code{NULL}.
#' @keywords internal
.resolve_vertex_attr <- function(g, candidates, label = "attribute",
                                  required = TRUE) {
    vatts <- igraph::vertex_attr_names(g)
    found <- candidates[candidates %in% vatts]
    if (length(found)) return(found[1L])
    if (required)
        stop("Could not find a ", label, " vertex attribute. ",
             "Tried: ", paste(candidates, collapse = ", "), ". ",
             "Available: ", paste(vatts, collapse = ", "), ".",
             call. = FALSE)
    NULL
}


#' Resolve x and y coordinate attribute names
#'
#' @param g An \code{igraph} object.
#' @param x_col Character(1) or \code{NULL}. Explicit override.
#' @param y_col Character(1) or \code{NULL}. Explicit override.
#' @param required Logical(1). Whether coordinates must be found.
#' @return A named list with elements \code{x_col} and \code{y_col}
#'   (character(1) each, or \code{NULL} when not required and not found).
#' @keywords internal
.resolve_coord_attrs <- function(g, x_col = NULL, y_col = NULL,
                                  required = TRUE) {
    if (is.null(x_col))
        x_col <- .resolve_vertex_attr(
            g, c("x", "cx", "centroid_x", "X"),
            label = "x-coordinate", required = required
        )
    if (is.null(y_col))
        y_col <- .resolve_vertex_attr(
            g, c("y", "cy", "centroid_y", "Y"),
            label = "y-coordinate", required = required
        )
    list(x_col = x_col, y_col = y_col)
}


#' Resolve cell-type attribute name
#'
#' @param g An \code{igraph} object.
#' @param cell_type_col Character(1) or \code{NULL}. Explicit override.
#' @param required Logical(1). Whether the attribute must be found.
#' @return Character(1) — the resolved attribute name, or \code{NULL}.
#' @keywords internal
.resolve_celltype_attr <- function(g, cell_type_col = NULL,
                                    required = TRUE) {
    if (!is.null(cell_type_col)) return(cell_type_col)
    .resolve_vertex_attr(
        g, c("cell_type", "type", "label"),
        label = "cell-type", required = required
    )
}


#' Assert that an object is an igraph
#'
#' @param g Object to check.
#' @param caller Character(1). Name of calling function for error messages.
#' @keywords internal
.assert_igraph <- function(g, caller = "") {
    if (!igraph::is.igraph(g))
        stop(if (nzchar(caller)) paste0(caller, ": "),
             "'g' must be an igraph object.", call. = FALSE)
}


#' Load an image file as an RGB raster array
#'
#' Reads PNG, JPEG, or whole-slide image formats (.svs, .ndpi, .tiff) and
#' returns a 3-D numeric array (height x width x channels) suitable for
#' \code{ggplot2::annotation_raster}.
#'
#' @param image Character(1) file path, or a numeric array already loaded.
#' @return A numeric array with dimensions (height, width, 3 or 4).
#' @keywords internal
.load_image_raster <- function(image) {
    if (is.array(image) || is.matrix(image))
        return(image)

    if (!is.character(image) || length(image) != 1L)
        stop("'image' must be a file path or a raster array.", call. = FALSE)

    if (!file.exists(image))
        stop("Image file not found: ", image, call. = FALSE)

    ext <- tolower(tools::file_ext(image))

    if (ext == "png") {
        if (!requireNamespace("png", quietly = TRUE))
            stop("Package 'png' is required to read PNG images. ",
                 "Install it with install.packages(\"png\").", call. = FALSE)
        return(png::readPNG(image))
    }

    if (ext %in% c("jpg", "jpeg")) {
        if (!requireNamespace("jpeg", quietly = TRUE))
            stop("Package 'jpeg' is required to read JPEG images. ",
                 "Install it with install.packages(\"jpeg\").", call. = FALSE)
        return(jpeg::readJPEG(image))
    }

    if (ext %in% c("svs", "ndpi", "tiff", "tif")) {
        openslide <- tryCatch(
            reticulate::import("openslide", delay_load = FALSE),
            error = function(e)
                stop("Python package 'openslide-python' is required to read ",
                     "whole-slide images. Install it with: ",
                     "pip install openslide-python", call. = FALSE)
        )
        np  <- reticulate::import("numpy", delay_load = FALSE)
        wsi <- openslide$OpenSlide(normalizePath(image))
        dim <- wsi$dimensions
        # Use highest-level thumbnail that is at least 2000px wide
        best_level <- as.integer(wsi$level_count - 1L)
        for (lv in seq(0L, wsi$level_count - 1L)) {
            lv_dim <- wsi$level_dimensions[[lv + 1L]]
            if (lv_dim[[1L]] <= 4000L) {
                best_level <- as.integer(lv)
                break
            }
        }
        region <- wsi$read_region(
            reticulate::tuple(0L, 0L),
            best_level,
            wsi$level_dimensions[[best_level + 1L]]
        )
        img_arr <- np$array(region$convert("RGB"))
        img_arr <- reticulate::py_to_r(img_arr)
        return(img_arr / 255)
    }

    stop("Unsupported image format: .", ext,
         ". Supported: png, jpg, jpeg, svs, ndpi, tiff, tif.", call. = FALSE)
}
