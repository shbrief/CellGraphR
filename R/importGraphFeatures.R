#' Import patient-level graph features as a SummarizedExperiment
#'
#' Reads the pre-computed patient-level cell-graph feature matrix (CSV) and
#' returns it as a \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#' ready for multi-modal downstream analysis.
#'
#' @param fname Character(1). Path to the CSV file of patient-level graph
#'   features. Expected columns: \code{patient_id} plus one column per
#'   feature.
#' @param colData A \code{data.frame} or \code{DataFrame} with additional
#'   sample-level annotations (rows must be patient IDs). \code{NULL} to skip.
#' @param feature_groups Logical(1). If \code{TRUE} (default), annotate each
#'   feature with its thematic group in \code{rowData}.
#'
#' @return A \code{SummarizedExperiment} object with:
#'   \describe{
#'     \item{assay \dQuote{graph_features}}{numeric matrix of dimension
#'       \emph{features} \eqn{\times} \emph{patients}.}
#'     \item{colData}{patient ID plus any additional annotations from
#'       \code{colData}.}
#'     \item{rowData}{feature name and, when \code{feature_groups = TRUE},
#'       a \code{feature_group} column classifying each feature into one of
#'       14 thematic groups (see \code{\link{.classify_graph_features}}).}
#'   }
#'
#' @examples
#' \dontrun{
#' csv_file <- system.file(
#'     "extdata", "example_graph_features.csv",
#'     package = "CellGraphR"
#' )
#' graph_se <- importGraphFeatures(csv_file)
#' dim(graph_se)       # features x patients
#' rowData(graph_se)   # feature annotations with group labels
#' }
#'
#' @importFrom utils read.csv
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @export
importGraphFeatures <- function(fname, colData = NULL, feature_groups = TRUE) {

    df <- utils::read.csv(fname, check.names = FALSE)
    df <- df[!is.na(df$patient_id) & nzchar(as.character(df$patient_id)), ]

    if (nrow(df) == 0L)
        stop("No valid rows found in '", fname,
             "'. Check that 'patient_id' column is present and non-empty.")

    rownames(df) <- df$patient_id
    df$patient_id <- NULL

    mat <- t(as.matrix(df))
    storage.mode(mat) <- "numeric"

    ## ---- column (sample) metadata ----------------------------------------
    col_df <- S4Vectors::DataFrame(patient_id = colnames(mat),
                                   row.names   = colnames(mat))

    if (!is.null(colData)) {
        colData <- as.data.frame(colData)
        shared  <- intersect(colnames(mat), rownames(colData))
        if (length(shared) == 0L)
            warning("No matching patients between 'colData' and feature matrix.")
        extra <- colData[shared, , drop = FALSE]
        col_df[shared, colnames(extra)] <- extra
    }

    ## ---- row (feature) metadata ------------------------------------------
    row_df <- S4Vectors::DataFrame(feature_name = rownames(mat),
                                   row.names     = rownames(mat))

    if (feature_groups) {
        row_df$feature_group <- .classify_graph_features(rownames(mat))
    }

    SummarizedExperiment::SummarizedExperiment(
        assays  = list(graph_features = mat),
        colData = col_df,
        rowData = row_df
    )
}


#' Add graph features to a MultiAssayExperiment
#'
#' Integrates a \code{SummarizedExperiment} of patient-level graph features
#' (produced by \code{\link{importGraphFeatures}}) into an existing
#' \code{\link[MultiAssayExperiment]{MultiAssayExperiment}}, aligning by
#' patient ID.
#'
#' @param mae A \code{MultiAssayExperiment}.
#' @param graph_se A \code{SummarizedExperiment} produced by
#'   \code{\link{importGraphFeatures}}, with patient IDs as column names.
#' @param experiment_name Character(1). Name of the new experiment slot in the
#'   returned \code{MultiAssayExperiment}. Default \code{"graph_features"}.
#'
#' @return A \code{MultiAssayExperiment} with the graph features added as an
#'   additional experiment.
#'
#' @details
#' Only patients present in \code{colData(mae)} are retained in
#' \code{graph_se}. The \code{sampleMap} is constructed automatically,
#' assuming that column names of \code{graph_se} correspond directly to the
#' primary sample identifiers in \code{mae}.
#'
#' @examples
#' \dontrun{
#' library(MultiAssayExperiment)
#' # mae <- curatedTCGAData::curatedTCGAData("PRAD", "RNASeq2GeneNorm*", ...)
#' # graph_se <- importGraphFeatures("path/to/features.csv")
#' mae2 <- addGraphToMAE(mae, graph_se)
#' experiments(mae2)
#' }
#'
#' @importFrom MultiAssayExperiment ExperimentList colData
#' @importFrom S4Vectors DataFrame
#' @export
addGraphToMAE <- function(mae, graph_se, experiment_name = "graph_features") {

    primary_ids <- rownames(MultiAssayExperiment::colData(mae))
    shared      <- intersect(colnames(graph_se), primary_ids)

    if (length(shared) == 0L)
        stop("No shared patient IDs between 'graph_se' and 'mae'. ",
             "Check that column names of 'graph_se' match primary identifiers in 'mae'.")

    if (length(shared) < ncol(graph_se))
        message(ncol(graph_se) - length(shared),
                " sample(s) in 'graph_se' not found in 'mae' and will be dropped.")

    graph_se_sub <- graph_se[, shared]

    smap <- S4Vectors::DataFrame(
        assay   = experiment_name,
        primary = shared,
        colname = shared
    )

    new_el <- MultiAssayExperiment::ExperimentList(
        setNames(list(graph_se_sub), experiment_name)
    )

    suppressWarnings(c(mae, new_el, sampleMap = smap))
}
