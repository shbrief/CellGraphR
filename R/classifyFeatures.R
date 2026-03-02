#' Classify graph feature names into thematic groups
#'
#' Assigns each feature name to one of 14 biologically motivated thematic
#' groups using regex pattern matching. Called automatically by
#' \code{\link{importGraphFeatures}} when \code{feature_groups = TRUE}.
#'
#' @param feature_names Character vector of feature names.
#'
#' @return Character vector of group labels, the same length as
#'   \code{feature_names}. Features that do not match any pattern are labelled
#'   \code{"other"}.
#'
#' @details
#' The 14 groups and the naming patterns that trigger them:
#' \describe{
#'   \item{global_topology}{n_cells, n_edges, mean/std/max/skewness/kurtosis of degree.}
#'   \item{cell_type_composition}{frac_* — cell-type fractions.}
#'   \item{cell_type_mixing}{edge_frac_* — fraction of edges between cell types.}
#'   \item{neighborhood_radius}{mean_deg_*_rN, mean_mix_*_rN, n_cells_*_rN — multi-scale (r20/40/80/100) features.}
#'   \item{per_celltype_graph}{mean_deg_*, mean_mix_*, std_deg_*, std_mix_* without radius suffix.}
#'   \item{cell_type_indicator}{is_* — binary presence flags.}
#'   \item{local_node}{degree, mixing_index, local_frac_*.}
#'   \item{voronoi}{voronoi_* — Voronoi tessellation metrics.}
#'   \item{nearest_neighbor_dist}{nn_dist_* — inter-cell-type nearest-neighbor distances.}
#'   \item{clustering}{clustering_coeff_* — local clustering coefficients.}
#'   \item{spatial_autocorrelation}{morans_I_* — Moran's I statistics.}
#'   \item{tumor_architecture}{tumor_*, invasion_front, interior_tumor, neoplastic_ecc/hull.}
#'   \item{immune_spatial}{immune_*, TIL_* — immune infiltration metrics.}
#'   \item{network_community}{modularity_*, n_connected_*, celltype_entropy, local_entropy.}
#'   \item{network_motif}{total_triangles, homotypic_triangles, heterotypic_triangles, homotypic_triangle_frac.}
#'   \item{spatial_distribution}{features ending in _dispersion or _ratio (not otherwise classified).}
#' }
#'
#' @keywords internal
.classify_graph_features <- function(feature_names) {

    groups <- rep("other", length(feature_names))

    groups[grepl("^(n_cells|n_edges|mean_degree|std_degree|max_degree|degree_skewness|degree_kurtosis)",
                 feature_names)] <- "global_topology"
    groups[grepl("^frac_", feature_names)]      <- "cell_type_composition"
    groups[grepl("^edge_frac_", feature_names)] <- "cell_type_mixing"
    groups[grepl("^mean_deg_.*_r[0-9]+$|^mean_mix_.*_r[0-9]+$|^n_cells_.*_r[0-9]+",
                 feature_names)] <- "neighborhood_radius"
    groups[grepl("^mean_deg_|^mean_mix_|^std_deg_|^std_mix_", feature_names) &
               !grepl("_r[0-9]+$", feature_names)] <- "per_celltype_graph"
    groups[grepl("^(is_)", feature_names)]          <- "cell_type_indicator"
    groups[grepl("^(degree$|mixing_index$|local_frac_)", feature_names)] <- "local_node"
    groups[grepl("^voronoi_", feature_names)]        <- "voronoi"
    groups[grepl("^nn_dist_", feature_names)]        <- "nearest_neighbor_dist"
    groups[grepl("^clustering_coeff_", feature_names)] <- "clustering"
    groups[grepl("^morans_I_", feature_names)]       <- "spatial_autocorrelation"
    groups[grepl("^tumor_|invasion_front|interior_tumor|neoplastic_ecc|neoplastic_hull",
                 feature_names)] <- "tumor_architecture"
    groups[grepl("^immune_|^TIL_", feature_names)]  <- "immune_spatial"
    groups[grepl("^modularity_|^n_connected_|^celltype_entropy|^local_entropy",
                 feature_names)] <- "network_community"
    groups[grepl("^(total_triangles|homotypic_triangles|heterotypic_triangles|homotypic_triangle_frac)",
                 feature_names)] <- "network_motif"
    groups[grepl("_(dispersion|ratio)$", feature_names) &
               groups == "other"] <- "spatial_distribution"

    groups
}
