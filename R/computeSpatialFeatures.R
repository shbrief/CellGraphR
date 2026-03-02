#' Extract advanced spatial relationship features
#'
#' Computes a comprehensive set of spatial point-pattern statistics from a list
#' of nucleus objects, including spatial autocorrelation (Moran's I),
#' multi-scale local density and diversity, neighborhood composition,
#' Hopkins clustering tendency, Ripley's K/L functions, and nearest-neighbor
#' statistics.
#'
#' @param nuclei_list A list of nucleus objects. Each element must have at
#'   minimum: \code{centroid} (numeric vector of length 2, x/y coordinates)
#'   and \code{type} (character, cell-type label).
#' @param image_dims Numeric vector \code{c(width, height)} specifying the
#'   observation window in pixels.
#' @param verbose Logical. Print progress messages.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{spatial_autocorr}{Moran's I and interpretation.}
#'     \item{multiscale}{Local density and diversity at scales 25, 50, 100,
#'       200 pixels.}
#'     \item{neighborhood}{k-NN homogeneity and diversity for k = 3, 5, 10,
#'       15.}
#'     \item{clustering}{Hopkins statistic and clustering interpretation.}
#'     \item{ripley_k}{Ripley K/L function deviation statistics.}
#'     \item{nearest_neighbor}{Clark-Evans index and k-NN distance summaries.}
#'   }
#'
#' @examples
#' \dontrun{
#' feats <- extract_advanced_spatial_features(nuclei_list, image_dims = c(1024, 1024))
#' feats$spatial_autocorr$morans_i
#' feats$nearest_neighbor$clark_evans_index
#' }
#'
#' @importFrom spatstat.geom owin ppp marks area nndist
#' @importFrom spatstat.explore Kest Lest
#' @export
extract_advanced_spatial_features <- function(nuclei_list, image_dims, verbose = FALSE) {

    if (verbose) cat("  Processing spatial features for", length(nuclei_list), "nuclei...\n")

    centroids <- t(sapply(nuclei_list, function(n) n$centroid))
    types     <- sapply(nuclei_list, function(n) n$type)

    W  <- spatstat.geom::owin(c(0, image_dims[1L]), c(0, image_dims[2L]))
    pp <- spatstat.geom::ppp(centroids[, 1L], centroids[, 2L],
                              window = W,
                              marks  = factor(types))

    if (verbose) cat("    Computing spatial autocorrelation...\n")
    spatial_autocorr <- compute_spatial_autocorrelation(pp)

    if (verbose) cat("    Multi-scale point pattern analysis...\n")
    multiscale <- compute_multiscale_analysis(pp)

    if (verbose) cat("    Neighborhood structure analysis...\n")
    neighborhood <- compute_neighborhood_features(nuclei_list, centroids, types)

    if (verbose) cat("    Spatial clustering analysis...\n")
    clustering <- compute_spatial_clustering(pp)

    if (verbose) cat("    Ripley K/L and nearest-neighbor analysis...\n")
    ripley_k        <- compute_ripley_k_analysis(pp)
    nearest_neighbor <- compute_nearest_neighbor_analysis(pp)

    if (verbose) cat("  Spatial feature extraction completed.\n")

    list(
        spatial_autocorr  = spatial_autocorr,
        multiscale        = multiscale,
        neighborhood      = neighborhood,
        clustering        = clustering,
        ripley_k          = ripley_k,
        nearest_neighbor  = nearest_neighbor
    )
}


# ---- spatial analysis helpers ------------------------------------------------

#' Compute spatial autocorrelation via Moran's I
#'
#' Uses inverse-distance weighting (with a 30th-percentile distance cutoff)
#' to construct the spatial weights matrix.
#'
#' @param pp A \code{spatstat.geom::ppp} marked point pattern.
#' @return Named list: \code{morans_i}, \code{expected_i}, \code{p_value}
#'   (always \code{NA}; requires permutation test), \code{interpretation}.
#' @keywords internal
compute_spatial_autocorrelation <- function(pp) {

    type_numeric <- as.numeric(spatstat.geom::marks(pp))

    if (length(unique(type_numeric)) < 2L || pp$n < 3L) {
        return(list(morans_i = 0, p_value = 1,
                    interpretation = "insufficient_data"))
    }

    coords     <- cbind(pp$x, pp$y)
    distances  <- as.matrix(dist(coords))
    cutoff     <- quantile(distances[upper.tri(distances)], 0.3)
    weights    <- 1 / (distances + 1)
    weights[distances > cutoff] <- 0
    diag(weights) <- 0

    n         <- length(type_numeric)
    mean_t    <- mean(type_numeric)
    deviations <- type_numeric - mean_t

    numerator         <- sum(weights * outer(deviations, deviations))
    sum_weights        <- sum(weights)
    sum_sq_deviations  <- sum(deviations^2)

    if (sum_weights == 0 || sum_sq_deviations == 0) {
        return(list(morans_i = 0, p_value = 1,
                    interpretation = "no_spatial_correlation"))
    }

    morans_i   <- (n / sum_weights) * (numerator / sum_sq_deviations)
    expected_i <- -1 / (n - 1L)

    interpretation <- if (morans_i > expected_i + 0.1) {
        "positive_autocorrelation"
    } else if (morans_i < expected_i - 0.1) {
        "negative_autocorrelation"
    } else {
        "random_distribution"
    }

    list(morans_i = morans_i, expected_i = expected_i,
         p_value = NA_real_, interpretation = interpretation)
}


#' Compute multi-scale local density and diversity
#'
#' For each spatial scale (radius in pixels), computes the mean and variance
#' of local nuclear density and Shannon diversity of cell types in the
#' neighbourhood.
#'
#' @param pp A \code{spatstat.geom::ppp} marked point pattern.
#' @return Named list with one entry per scale (\code{scale_25}, \code{scale_50},
#'   \code{scale_100}, \code{scale_200}), each containing
#'   \code{mean_density}, \code{var_density}, \code{max_density},
#'   \code{mean_diversity}, \code{var_diversity}.
#' @keywords internal
compute_multiscale_analysis <- function(pp) {

    scales <- c(25, 50, 100, 200)
    result <- list()

    for (r in scales) {
        local_density   <- numeric(pp$n)
        local_diversity <- numeric(pp$n)

        for (i in seq_len(pp$n)) {
            d      <- sqrt((pp$x - pp$x[i])^2 + (pp$y - pp$y[i])^2)
            nb_idx <- which(d <= r & d > 0)
            local_density[i] <- length(nb_idx) / (pi * r^2)

            if (length(nb_idx) > 0L) {
                tab   <- table(spatstat.geom::marks(pp)[nb_idx])
                props <- tab / sum(tab)
                local_diversity[i] <- -sum(props * log(props + 1e-10))
            }
        }

        result[[paste0("scale_", r)]] <- list(
            mean_density  = mean(local_density),
            var_density   = var(local_density),
            max_density   = max(local_density),
            mean_diversity = mean(local_diversity),
            var_diversity  = var(local_diversity)
        )
    }
    result
}


#' Compute k-NN neighbourhood homogeneity and diversity
#'
#' @param nuclei_list List of nucleus objects.
#' @param centroids Numeric matrix of centroid coordinates.
#' @param types Character vector of cell-type labels.
#' @return Named list with one entry per k (\code{k3}, \code{k5}, \code{k10},
#'   \code{k15}), each containing \code{mean_homogeneity},
#'   \code{var_homogeneity}, \code{mean_diversity}, \code{var_diversity},
#'   \code{max_homogeneity}, \code{min_diversity}.
#' @keywords internal
compute_neighborhood_features <- function(nuclei_list, centroids, types) {

    n_nuclei  <- nrow(centroids)
    distances <- as.matrix(dist(centroids))
    result    <- list()

    for (k in c(3L, 5L, 10L, 15L)) {
        if (k >= n_nuclei) next

        homogeneity <- numeric(n_nuclei)
        diversity   <- numeric(n_nuclei)

        for (i in seq_len(n_nuclei)) {
            nb_idx      <- order(distances[i, ])[seq(2L, k + 1L)]
            nb_types    <- types[nb_idx]
            tab         <- table(nb_types)
            props       <- tab / sum(tab)
            homogeneity[i] <- max(props)
            diversity[i]   <- -sum(props * log(props + 1e-10))
        }

        result[[paste0("k", k)]] <- list(
            mean_homogeneity = mean(homogeneity),
            var_homogeneity  = var(homogeneity),
            mean_diversity   = mean(diversity),
            var_diversity    = var(diversity),
            max_homogeneity  = max(homogeneity),
            min_diversity    = min(diversity)
        )
    }
    result
}


#' Compute Hopkins statistic for spatial clustering tendency
#'
#' @param pp A \code{spatstat.geom::ppp} point pattern.
#' @return Named list: \code{hopkins_statistic} and
#'   \code{clustering_interpretation} (\code{"highly_clustered"},
#'   \code{"regular_pattern"}, or \code{"random_pattern"}).
#' @keywords internal
compute_spatial_clustering <- function(pp) {

    n_sample <- min(30L, pp$n)

    if (n_sample < 10L) {
        return(list(hopkins_statistic = NA_real_,
                    clustering_interpretation = "insufficient_data"))
    }

    sample_idx    <- sample(pp$n, n_sample)
    sample_pts    <- cbind(pp$x[sample_idx], pp$y[sample_idx])
    random_pts    <- cbind(runif(n_sample, pp$window$xrange[1L], pp$window$xrange[2L]),
                           runif(n_sample, pp$window$yrange[1L], pp$window$yrange[2L]))

    real_nn   <- numeric(n_sample)
    random_nn <- numeric(n_sample)

    for (i in seq_len(n_sample)) {
        d_real            <- sqrt((pp$x - sample_pts[i, 1L])^2 + (pp$y - sample_pts[i, 2L])^2)
        d_real[sample_idx[i]] <- Inf
        real_nn[i]        <- min(d_real)

        d_rand            <- sqrt((pp$x - random_pts[i, 1L])^2 + (pp$y - random_pts[i, 2L])^2)
        random_nn[i]      <- min(d_rand)
    }

    U      <- sum(random_nn)
    W      <- sum(real_nn)
    hopkins <- U / (U + W)

    interpretation <- if (hopkins > 0.75) {
        "highly_clustered"
    } else if (hopkins < 0.25) {
        "regular_pattern"
    } else {
        "random_pattern"
    }

    list(hopkins_statistic = hopkins, clustering_interpretation = interpretation)
}


#' Compute Ripley's K and L function statistics
#'
#' @param pp A \code{spatstat.geom::ppp} point pattern.
#' @return Named list: \code{k_max_deviation}, \code{k_clustering_strength},
#'   \code{k_pattern}, \code{l_max_deviation}.
#' @importFrom spatstat.explore Kest Lest
#' @keywords internal
compute_ripley_k_analysis <- function(pp) {

    if (pp$n < 10L) {
        return(list(k_max_deviation = NA_real_, k_clustering_strength = NA_real_,
                    k_pattern = "insufficient_data", l_max_deviation = NA_real_))
    }

    tryCatch({
        K <- spatstat.explore::Kest(pp, correction = "isotropic")
        L <- spatstat.explore::Lest(pp, correction = "isotropic")

        dev   <- K$iso - K$theo
        k_max <- max(dev, na.rm = TRUE)
        list(
            k_max_deviation      = k_max,
            k_clustering_strength = mean(dev[dev > 0], na.rm = TRUE),
            k_pattern             = if (k_max > 0) "clustered" else "regular_or_random",
            l_max_deviation       = max(abs(L$iso - L$theo), na.rm = TRUE)
        )
    }, error = function(e) {
        list(k_max_deviation = NA_real_, k_clustering_strength = NA_real_,
             k_pattern = "computation_failed", l_max_deviation = NA_real_)
    })
}


#' Compute comprehensive nearest-neighbor statistics
#'
#' @param pp A \code{spatstat.geom::ppp} point pattern.
#' @return Named list including \code{nn_mean}, \code{nn_median},
#'   \code{nn_std}, \code{nn_min}, \code{nn_max}, \code{nn_range},
#'   \code{clark_evans_index}, \code{ce_interpretation}, \code{nn_cv}, and
#'   k-th nearest-neighbor means/stds for k = 2 to min(5, n - 1).
#' @importFrom spatstat.geom nndist area
#' @keywords internal
compute_nearest_neighbor_analysis <- function(pp) {

    if (pp$n < 2L) {
        empty <- stats::setNames(rep(NA_real_, 8L),
                                 c("nn_mean", "nn_median", "nn_std", "nn_min",
                                   "nn_max", "nn_range", "clark_evans_index", "nn_cv"))
        return(c(as.list(empty), list(ce_interpretation = "single_point")))
    }

    nn  <- spatstat.geom::nndist(pp)
    dens <- pp$n / spatstat.geom::area(pp$window)

    ce_index <- mean(nn) / (1 / (2 * sqrt(dens)))
    ce_interp <- if (ce_index > 1.15) "dispersed" else if (ce_index < 0.85) "clustered" else "random"

    result <- list(
        nn_mean            = mean(nn),
        nn_median          = median(nn),
        nn_std             = stats::sd(nn),
        nn_min             = min(nn),
        nn_max             = max(nn),
        nn_range           = diff(range(nn)),
        clark_evans_index  = ce_index,
        ce_interpretation  = ce_interp,
        nn_cv              = stats::sd(nn) / mean(nn)
    )

    for (k in seq(2L, min(5L, pp$n - 1L))) {
        knn <- spatstat.geom::nndist(pp, k = k)
        result[[paste0("nn", k, "_mean")]] <- mean(knn)
        result[[paste0("nn", k, "_std")]]  <- stats::sd(knn)
    }

    result
}
