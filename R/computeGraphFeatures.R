#' Extract graph-based tissue organization features
#'
#' Computes multiple graph representations of cell spatial organization from a
#' list of nucleus objects and returns a nested list of graph-topology metrics.
#' Five graph types are computed: Delaunay triangulation, k-nearest-neighbor
#' (k = 3, 5, 7, 10), minimum spanning tree, Gabriel graph, and relative
#' neighborhood graph.
#'
#' @param nuclei_list A list of nucleus objects. Each element must have at
#'   minimum: \code{centroid} (numeric vector of length 2, x/y coordinates)
#'   and \code{type} (character, cell-type label).
#' @param image_dims Numeric vector \code{c(width, height)} specifying the
#'   image dimensions in pixels. Currently unused by the graph functions but
#'   included for API consistency with
#'   \code{\link{extract_advanced_spatial_features}}.
#' @param verbose Logical. Print progress messages.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{delaunay}{Graph metrics from the Delaunay triangulation.}
#'     \item{knn}{Named sub-list with one entry per k (\code{k3}, \code{k5},
#'       \code{k7}, \code{k10}), each containing k-NN graph metrics.}
#'     \item{mst}{Minimum spanning tree metrics.}
#'     \item{gabriel}{Gabriel graph metrics.}
#'     \item{rng}{Relative neighborhood graph metrics.}
#'   }
#'
#' @examples
#' \dontrun{
#' # nuclei_list is the 'nuc' element of a HoVerNet JSON parsed with jsonToSe
#' feats <- extract_graph_features(nuclei_list, image_dims = c(1024, 1024))
#' feats$delaunay$transitivity
#' feats$knn$k5$clustering_coefficient
#' }
#'
#' @export
extract_graph_features <- function(nuclei_list, image_dims, verbose = FALSE) {

    if (verbose) cat("  Processing graph features for", length(nuclei_list), "nuclei...\n")

    centroids <- t(sapply(nuclei_list, function(n) n$centroid))
    types     <- sapply(nuclei_list, function(n) n$type)

    graph_features <- list()

    if (verbose) cat("    Computing Delaunay triangulation...\n")
    graph_features$delaunay <- compute_delaunay_graph(centroids, types)

    if (verbose) cat("    Computing k-NN graphs...\n")
    knn_features <- list()
    for (k in c(3, 5, 7, 10)) {
        if (k < nrow(centroids)) {
            knn_features[[paste0("k", k)]] <- compute_knn_graph(centroids, types, k)
        }
    }
    graph_features$knn <- knn_features

    if (verbose) cat("    Computing minimum spanning tree...\n")
    graph_features$mst <- compute_mst_features(centroids, types)

    if (verbose) cat("    Computing Gabriel graph...\n")
    graph_features$gabriel <- compute_gabriel_graph(centroids, types)

    if (verbose) cat("    Computing relative neighborhood graph...\n")
    graph_features$rng <- compute_relative_neighborhood_graph(centroids, types)

    if (verbose) cat("  Graph feature extraction completed.\n")

    graph_features
}


# ---- graph constructors ------------------------------------------------------

# ---- adjacency builders (internal) -------------------------------------------

.build_delaunay_adj <- function(centroids) {
    n         <- nrow(centroids)
    distances <- as.matrix(dist(centroids))
    threshold <- quantile(distances[upper.tri(distances)], 0.2)
    adj       <- (distances <= threshold) & (distances > 0)
    for (i in seq_len(n)) {
        if (sum(adj[i, ]) < 2L) {
            nn         <- order(distances[i, ])[seq(2L, 3L)]
            adj[i, nn] <- TRUE
            adj[nn, i] <- TRUE
        }
    }
    adj | t(adj)
}

.build_knn_adj <- function(centroids, k) {
    n         <- nrow(centroids)
    distances <- as.matrix(dist(centroids))
    adj       <- matrix(FALSE, n, n)
    for (i in seq_len(n))
        adj[i, order(distances[i, ])[seq(2L, k + 1L)]] <- TRUE
    adj & t(adj)
}

.build_mst_edges <- function(centroids) {
    n           <- nrow(centroids)
    distances   <- as.matrix(dist(centroids))
    mst_edges   <- matrix(0L, n - 1L, 2L)
    mst_weights <- numeric(n - 1L)
    visited     <- logical(n)
    visited[1L] <- TRUE
    for (i in seq_len(n - 1L)) {
        min_dist <- Inf
        min_edge <- c(0L, 0L)
        for (v1 in which(visited))
            for (v2 in which(!visited))
                if (distances[v1, v2] < min_dist) {
                    min_dist <- distances[v1, v2]
                    min_edge <- c(v1, v2)
                }
        mst_edges[i, ]  <- min_edge
        mst_weights[i]  <- min_dist
        visited[min_edge[2L]] <- TRUE
    }
    list(edges = mst_edges, weights = mst_weights)
}

.build_gabriel_adj <- function(centroids) {
    n   <- nrow(centroids)
    adj <- matrix(FALSE, n, n)
    for (i in seq_len(n - 1L)) {
        for (j in seq(i + 1L, n)) {
            cx     <- (centroids[i, 1L] + centroids[j, 1L]) / 2
            cy     <- (centroids[i, 2L] + centroids[j, 2L]) / 2
            radius <- sqrt((centroids[i, 1L] - centroids[j, 1L])^2 +
                               (centroids[i, 2L] - centroids[j, 2L])^2) / 2
            inside <- any(vapply(seq_len(n), function(k) {
                if (k == i || k == j) return(FALSE)
                sqrt((centroids[k, 1L] - cx)^2 + (centroids[k, 2L] - cy)^2) < radius
            }, logical(1L)))
            if (!inside) adj[i, j] <- adj[j, i] <- TRUE
        }
    }
    adj
}

.build_rng_adj <- function(centroids) {
    n   <- nrow(centroids)
    adj <- matrix(FALSE, n, n)
    for (i in seq_len(n - 1L)) {
        for (j in seq(i + 1L, n)) {
            edge_len   <- sqrt((centroids[i, 1L] - centroids[j, 1L])^2 +
                                   (centroids[i, 2L] - centroids[j, 2L])^2)
            lune_empty <- !any(vapply(seq_len(n), function(k) {
                if (k == i || k == j) return(FALSE)
                di <- sqrt((centroids[k, 1L] - centroids[i, 1L])^2 +
                               (centroids[k, 2L] - centroids[i, 2L])^2)
                dj <- sqrt((centroids[k, 1L] - centroids[j, 1L])^2 +
                               (centroids[k, 2L] - centroids[j, 2L])^2)
                max(di, dj) < edge_len
            }, logical(1L)))
            if (lune_empty) adj[i, j] <- adj[j, i] <- TRUE
        }
    }
    adj
}

# ------------------------------------------------------------------------------

#' Compute Delaunay triangulation graph features
#'
#' Approximates a Delaunay triangulation via adaptive distance thresholding
#' and returns graph-topology metrics.
#'
#' @param centroids Numeric matrix (\emph{n} \eqn{\times} 2) of cell centroid
#'   coordinates.
#' @param types Character vector of length \emph{n} with cell-type labels.
#'
#' @return Named list: \code{n_edges}, \code{density}, \code{mean_degree},
#'   \code{max_degree}, \code{degree_variance}, \code{transitivity},
#'   \code{assortativity}.
#'
#' @keywords internal
compute_delaunay_graph <- function(centroids, types) {

    n_nuclei <- nrow(centroids)

    if (n_nuclei < 3) {
        return(list(n_edges = 0, density = 0, mean_degree = 0, max_degree = 0,
                    degree_variance = 0, transitivity = 0, assortativity = NA))
    }

    tryCatch({
        adjacency <- .build_delaunay_adj(centroids)

        n_edges  <- sum(adjacency) / 2
        density  <- n_edges / (n_nuclei * (n_nuclei - 1) / 2)
        degrees  <- rowSums(adjacency)

        list(
            n_edges        = n_edges,
            density        = density,
            mean_degree    = mean(degrees),
            max_degree     = max(degrees),
            degree_variance = var(degrees),
            transitivity   = compute_transitivity(adjacency),
            assortativity  = compute_assortativity(adjacency, types)
        )
    }, error = function(e) {
        list(n_edges = 0, density = 0, mean_degree = 0, max_degree = 0,
             degree_variance = 0, transitivity = 0, assortativity = NA)
    })
}


#' Compute k-nearest-neighbor graph features
#'
#' Builds a mutual k-NN graph (an edge exists only if both endpoints list each
#' other as a k-nearest neighbor) and returns topology metrics.
#'
#' @param centroids Numeric matrix of centroid coordinates.
#' @param types Character vector of cell-type labels.
#' @param k Integer. Number of nearest neighbors.
#'
#' @return Named list: \code{n_edges}, \code{mean_degree},
#'   \code{clustering_coefficient}, \code{degree_variance},
#'   \code{assortativity}, \code{path_length}.
#'
#' @keywords internal
compute_knn_graph <- function(centroids, types, k) {

    n_nuclei <- nrow(centroids)

    if (k >= n_nuclei || n_nuclei < 3L) {
        return(list(n_edges = 0, mean_degree = 0, clustering_coefficient = 0,
                    degree_variance = 0, assortativity = NA, path_length = Inf))
    }

    adjacency <- .build_knn_adj(centroids, k)

    n_edges <- sum(adjacency) / 2
    degrees <- rowSums(adjacency)

    # Local clustering coefficients
    cc <- numeric(n_nuclei)
    for (i in seq_len(n_nuclei)) {
        nb <- which(adjacency[i, ])
        if (length(nb) >= 2L) {
            triangles <- 0L
            for (j in seq_len(length(nb) - 1L)) {
                for (l in seq(j + 1L, length(nb))) {
                    if (adjacency[nb[j], nb[l]]) triangles <- triangles + 1L
                }
            }
            possible <- length(nb) * (length(nb) - 1L) / 2
            cc[i] <- triangles / possible
        }
    }

    list(
        n_edges               = n_edges,
        mean_degree           = mean(degrees),
        clustering_coefficient = mean(cc, na.rm = TRUE),
        degree_variance       = var(degrees),
        assortativity         = compute_assortativity(adjacency, types),
        path_length           = compute_average_path_length(adjacency)
    )
}


#' Compute minimum spanning tree features
#'
#' Uses Prim's algorithm to construct the MST and returns edge-length and
#' type-mixing statistics.
#'
#' @param centroids Numeric matrix of centroid coordinates.
#' @param types Character vector of cell-type labels.
#'
#' @return Named list: \code{total_length}, \code{mean_edge_length},
#'   \code{max_edge_length}, \code{edge_length_variance}, \code{type_mixing},
#'   \code{mst_diameter}.
#'
#' @keywords internal
compute_mst_features <- function(centroids, types) {

    n_nuclei <- nrow(centroids)

    if (n_nuclei < 2L) {
        return(list(total_length = 0, mean_edge_length = 0, max_edge_length = 0,
                    edge_length_variance = 0, type_mixing = 0, mst_diameter = 0))
    }

    mst         <- .build_mst_edges(centroids)
    mst_edges   <- mst$edges
    mst_weights <- mst$weights

    type_mixing <- mean(types[mst_edges[, 1L]] != types[mst_edges[, 2L]])

    list(
        total_length        = sum(mst_weights),
        mean_edge_length    = mean(mst_weights),
        max_edge_length     = max(mst_weights),
        edge_length_variance = var(mst_weights),
        type_mixing         = type_mixing,
        mst_diameter        = compute_mst_diameter(mst_edges, mst_weights)
    )
}


#' Compute Gabriel graph features
#'
#' A Gabriel graph edge exists between two points if the open disc with that
#' edge as its diameter contains no other points.
#'
#' @param centroids Numeric matrix of centroid coordinates.
#' @param types Character vector of cell-type labels (currently unused but
#'   retained for API consistency).
#'
#' @return Named list: \code{n_edges}, \code{density}, \code{mean_degree},
#'   \code{transitivity}.
#'
#' @keywords internal
compute_gabriel_graph <- function(centroids, types) {

    n_nuclei <- nrow(centroids)

    if (n_nuclei < 3L) {
        return(list(n_edges = 0, density = 0, mean_degree = 0, transitivity = 0))
    }

    adjacency <- .build_gabriel_adj(centroids)

    degrees <- rowSums(adjacency)
    list(
        n_edges     = sum(adjacency) / 2,
        density     = sum(adjacency) / 2 / (n_nuclei * (n_nuclei - 1L) / 2),
        mean_degree = mean(degrees),
        transitivity = compute_transitivity(adjacency)
    )
}


#' Compute relative neighborhood graph features
#'
#' An RNG edge exists between two points if no third point lies within the
#' "lune" region (intersection of the two discs centred at the endpoints with
#' radius equal to the edge length).
#'
#' @param centroids Numeric matrix of centroid coordinates.
#' @param types Character vector of cell-type labels (currently unused).
#'
#' @return Named list: \code{n_edges}, \code{density}, \code{mean_degree},
#'   \code{degree_variance}.
#'
#' @keywords internal
compute_relative_neighborhood_graph <- function(centroids, types) {

    n_nuclei <- nrow(centroids)

    if (n_nuclei < 3L) {
        return(list(n_edges = 0, density = 0, mean_degree = 0, degree_variance = 0))
    }

    adjacency <- .build_rng_adj(centroids)

    degrees <- rowSums(adjacency)
    list(
        n_edges        = sum(adjacency) / 2,
        density        = sum(adjacency) / 2 / (n_nuclei * (n_nuclei - 1L) / 2),
        mean_degree    = mean(degrees),
        degree_variance = var(degrees)
    )
}


# ---- graph analysis helpers --------------------------------------------------

#' Compute global transitivity (clustering coefficient)
#'
#' @param adjacency Logical or binary adjacency matrix.
#' @return Numeric scalar in \eqn{[0, 1]}.
#' @keywords internal
compute_transitivity <- function(adjacency) {
    adjacency <- adjacency != 0  # coerce numeric 0/1 to logical
    n <- nrow(adjacency)
    total_triplets  <- 0L
    closed_triplets <- 0L
    for (i in seq_len(n)) {
        nb <- which(adjacency[i, ])
        if (length(nb) >= 2L) {
            total_triplets <- total_triplets + length(nb) * (length(nb) - 1L) / 2L
            for (j in seq_len(length(nb) - 1L)) {
                for (k in seq(j + 1L, length(nb))) {
                    if (adjacency[nb[j], nb[k]]) closed_triplets <- closed_triplets + 1L
                }
            }
        }
    }
    if (total_triplets == 0L) return(0)
    closed_triplets / total_triplets
}


#' Compute type assortativity (homophily)
#'
#' @param adjacency Adjacency matrix.
#' @param types Character (or numeric) vector of node type labels.
#' @return Numeric assortativity coefficient, or \code{NA} if no edges.
#' @keywords internal
compute_assortativity <- function(adjacency, types) {
    if (!is.numeric(types)) {
        u <- unique(types)
        types <- setNames(seq_along(u), u)[types]
    }
    n           <- length(types)
    total_edges <- sum(adjacency) / 2
    if (total_edges == 0) return(NA_real_)

    sum_xy <- sum_x <- sum_y <- sum_x2 <- sum_y2 <- 0
    for (i in seq_len(n - 1L)) {
        for (j in seq(i + 1L, n)) {
            if (adjacency[i, j]) {
                xi <- types[i]; yj <- types[j]
                sum_xy <- sum_xy + 2 * xi * yj
                sum_x  <- sum_x  + xi + yj
                sum_y  <- sum_y  + yj + xi
                sum_x2 <- sum_x2 + xi^2 + yj^2
                sum_y2 <- sum_y2 + yj^2 + xi^2
            }
        }
    }
    m   <- total_edges * 2
    mxy <- sum_xy / m; mx <- sum_x / m; my <- sum_y / m
    mx2 <- sum_x2 / m; my2 <- sum_y2 / m
    denom <- sqrt((mx2 - mx^2) * (my2 - my^2))
    if (denom == 0) return(0)
    (mxy - mx * my) / denom
}


#' Compute average shortest path length (Floyd-Warshall)
#'
#' @param adjacency Adjacency matrix.
#' @return Mean finite pairwise shortest-path length, or \code{Inf} if the
#'   graph is disconnected with no finite paths.
#' @keywords internal
compute_average_path_length <- function(adjacency) {
    n    <- nrow(adjacency)
    dist <- matrix(Inf, n, n)
    diag(dist) <- 0
    dist[adjacency > 0] <- 1
    for (k in seq_len(n))
        for (i in seq_len(n))
            for (j in seq_len(n))
                if (dist[i, k] + dist[k, j] < dist[i, j])
                    dist[i, j] <- dist[i, k] + dist[k, j]
    finite_d <- dist[is.finite(dist) & dist > 0]
    if (length(finite_d) == 0L) return(Inf)
    mean(finite_d)
}


#' Compute MST diameter via double-DFS
#'
#' @param mst_edges Integer matrix (\emph{(n-1)} \eqn{\times} 2) of MST edges.
#' @param mst_weights Numeric vector of corresponding edge weights.
#' @return Numeric MST diameter (longest weighted path).
#' @keywords internal
compute_mst_diameter <- function(mst_edges, mst_weights) {
    if (nrow(mst_edges) == 0L) return(0)
    n <- max(mst_edges)
    adj <- vector("list", n)
    for (i in seq_len(nrow(mst_edges))) {
        u <- mst_edges[i, 1L]; v <- mst_edges[i, 2L]; w <- mst_weights[i]
        adj[[u]] <- rbind(adj[[u]], c(v, w))
        adj[[v]] <- rbind(adj[[v]], c(u, w))
    }
    find_farthest_node(adj, find_farthest_node(adj, 1L)$node)$distance
}


#' DFS to find the farthest node from a starting node
#'
#' @param adj_list Adjacency list (list of 2-column matrices: neighbor, weight).
#' @param start Integer start node index.
#' @return List with \code{node} (integer) and \code{distance} (numeric).
#' @keywords internal
find_farthest_node <- function(adj_list, start) {
    n       <- length(adj_list)
    visited <- logical(n)
    dists   <- numeric(n)

    dfs <- function(node, d) {
        visited[node] <<- TRUE
        dists[node]   <<- d
        if (!is.null(adj_list[[node]])) {
            for (i in seq_len(nrow(adj_list[[node]]))) {
                nb <- adj_list[[node]][i, 1L]
                w  <- adj_list[[node]][i, 2L]
                if (!visited[nb]) dfs(nb, d + w)
            }
        }
    }
    dfs(start, 0)
    list(node = which.max(dists), distance = max(dists))
}


# ---- advanced graph features -------------------------------------------------

#' Compute graph spectral features from Laplacian eigenvalues
#'
#' @param adjacency Adjacency matrix.
#' @return Named list: \code{algebraic_connectivity} (second-smallest
#'   eigenvalue of the Laplacian), \code{spectral_gap}, \code{trace},
#'   \code{spectral_radius}.
#' @keywords internal
compute_spectral_features <- function(adjacency) {
    tryCatch({
        L    <- diag(rowSums(adjacency)) - adjacency
        evs  <- sort(eigen(L, only.values = TRUE)$values)
        list(
            algebraic_connectivity = evs[2L],
            spectral_gap           = evs[length(evs)] - evs[length(evs) - 1L],
            trace                  = sum(evs),
            spectral_radius        = max(abs(evs))
        )
    }, error = function(e) {
        list(algebraic_connectivity = 0, spectral_gap = 0, trace = 0, spectral_radius = 0)
    })
}


#' Count graph motifs (triangles, 3-stars, open 3-paths)
#'
#' @param adjacency Adjacency matrix.
#' @return Named list: \code{triangles}, \code{stars}, \code{paths}.
#' @keywords internal
compute_motif_counts <- function(adjacency) {
    adjacency <- adjacency != 0  # coerce numeric 0/1 to logical
    n <- nrow(adjacency)
    if (n < 3L) return(list(triangles = 0L, stars = 0L, paths = 0L))

    triangles <- 0L
    for (i in seq_len(n - 2L))
        for (j in seq(i + 1L, n - 1L))
            for (k in seq(j + 1L, n))
                if (adjacency[i, j] && adjacency[j, k] && adjacency[i, k])
                    triangles <- triangles + 1L

    stars <- sum(sapply(seq_len(n), function(i) {
        d <- sum(adjacency[i, ])
        if (d >= 3L) choose(d, 3L) else 0L
    }))

    paths <- 0L
    for (i in seq_len(n)) {
        nb <- which(adjacency[i, ])
        if (length(nb) >= 2L)
            for (j in seq_len(length(nb) - 1L))
                for (k in seq(j + 1L, length(nb)))
                    if (!adjacency[nb[j], nb[k]]) paths <- paths + 1L
    }

    list(triangles = triangles, stars = stars, paths = paths)
}
