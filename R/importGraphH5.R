#' Import cell-graph topology from a PyTables HDF5 file
#'
#' Reads the \code{nodes} and \code{edges} datasets from a PyTables-formatted
#' \code{.h5} file produced by the cell-graph construction pipeline and
#' optionally converts the result to an \code{igraph} object.
#'
#' @param fname Character(1). Path to a \code{.h5} graph file. The file must
#'   contain two PyTables datasets: \code{nodes} (column \code{node_id}) and
#'   \code{edges} (columns \code{source} and \code{target}).
#' @param as_igraph Logical(1). If \code{TRUE} (default), return an
#'   \code{igraph} object. If \code{FALSE}, return a list with \code{nodes}
#'   and \code{edges} data frames.
#' @param patient_id Character(1) or \code{NULL}. Patient identifier stored as
#'   a graph-level attribute. Inferred from the TCGA barcode in \code{fname}
#'   when \code{NULL}.
#'
#' @return An undirected \code{igraph} object (or list when
#'   \code{as_igraph = FALSE}) with the following graph-level attributes:
#'   \describe{
#'     \item{patient_id}{TCGA patient barcode.}
#'     \item{source_file}{Absolute path of the source file.}
#'     \item{n_nodes}{Number of cell nodes.}
#'     \item{n_edges}{Number of edges.}
#'   }
#'
#' @details
#' The \code{.h5} files use the PyTables HDF5 format written by
#' \code{pandas.DataFrame.to_hdf()}. This function reads them via
#' \pkg{reticulate} using \code{h5py} to access the raw HDF5 datasets and
#' \code{pandas} to construct data frames — the \code{tables} (PyTables)
#' back-end is \emph{not} required. Ensure that a Python environment with
#' \code{h5py >= 3.0} and \code{pandas >= 1.3} is available before calling
#' this function (e.g. \code{pip install h5py pandas}).
#'
#' @examples
#' \dontrun{
#' h5_file <- system.file("extdata", "example_graph.h5", package = "CellGraphR")
#' g <- importGraphH5(h5_file)
#' igraph::vcount(g)   # number of cell nodes
#' igraph::ecount(g)   # number of cell-cell edges
#' igraph::graph_attr(g, "patient_id")
#' }
#'
#' @importFrom reticulate import py_to_r
#' @importFrom igraph graph_from_data_frame graph_attr vcount ecount
#' @export
importGraphH5 <- function(fname, as_igraph = TRUE, patient_id = NULL) {

    fname <- normalizePath(fname, mustWork = TRUE)

    if (is.null(patient_id)) {
        patient_id <- regmatches(
            basename(fname),
            regexpr("TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}", basename(fname))
        )
        if (length(patient_id) == 0L) patient_id <- NA_character_
    }

    h5py <- reticulate::import("h5py", delay_load = FALSE)
    pd   <- reticulate::import("pandas", delay_load = FALSE)
    np   <- reticulate::import("numpy", delay_load = FALSE)

    hf <- h5py$File(fname, "r")
    on.exit(hf$close(), add = TRUE)

    nodes <- reticulate::py_to_r(pd$DataFrame(np$array(hf[["nodes/table"]])))
    edges <- reticulate::py_to_r(pd$DataFrame(np$array(hf[["edges/table"]])))

    if (!as_igraph) {
        return(list(nodes = nodes, edges = edges, patient_id = patient_id))
    }

    g <- igraph::graph_from_data_frame(
        d        = edges[, c("source", "target")],
        directed = FALSE,
        vertices = nodes
    )

    igraph::graph_attr(g, "patient_id")  <- patient_id
    igraph::graph_attr(g, "source_file") <- fname
    igraph::graph_attr(g, "n_nodes")     <- igraph::vcount(g)
    igraph::graph_attr(g, "n_edges")     <- igraph::ecount(g)

    g
}


#' Batch import of cell-graph topology files
#'
#' Applies \code{\link{importGraphH5}} to a vector of file paths and returns a
#' named list of \code{igraph} objects (or raw node/edge lists).
#'
#' @param fnames Character vector of paths to \code{.h5} graph files.
#' @param as_igraph Logical(1). Passed to \code{\link{importGraphH5}}.
#' @param verbose Logical(1). Print progress messages.
#'
#' @return A named list (names = TCGA patient barcodes) of \code{igraph}
#'   objects or raw lists, depending on \code{as_igraph}.
#'
#' @examples
#' \dontrun{
#' h5_files <- list.files(
#'     system.file("extdata", package = "CellGraphR"),
#'     pattern = "_graph\\.h5$",
#'     full.names = TRUE
#' )
#' graphs <- importAllGraphH5(h5_files)
#' sapply(graphs, igraph::vcount)   # cells per slide
#' }
#'
#' @importFrom reticulate import
#' @export
importAllGraphH5 <- function(fnames, as_igraph = TRUE, verbose = TRUE) {

    patient_ids <- regmatches(
        basename(fnames),
        regexpr("TCGA-[0-9A-Z]{2}-[0-9A-Z]{4}", basename(fnames))
    )

    dup_idx <- which(duplicated(patient_ids))
    if (length(dup_idx) > 0L) {
        message("Duplicate patient IDs detected; keeping first slide per patient: ",
                paste(patient_ids[dup_idx], collapse = ", "))
        fnames      <- fnames[-dup_idx]
        patient_ids <- patient_ids[-dup_idx]
    }

    result        <- vector("list", length(fnames))
    names(result) <- patient_ids

    for (i in seq_along(fnames)) {
        if (verbose) message(i, " / ", length(fnames), ": importing ", basename(fnames[i]))
        result[[i]] <- importGraphH5(fnames[i],
                                     as_igraph  = as_igraph,
                                     patient_id = patient_ids[i])
    }

    result
}
