test_that("extract_graph_features returns expected structure", {
    set.seed(42L)
    n <- 20L
    nuclei_list <- lapply(seq_len(n), function(i) {
        list(centroid = runif(2L, 0, 100),
             type     = sample(c("A", "B", "C"), 1L))
    })

    feats <- extract_graph_features(nuclei_list, image_dims = c(200, 200))

    expect_type(feats, "list")
    expect_named(feats, c("delaunay", "knn", "mst", "gabriel", "rng"))

    # Delaunay metrics
    expect_true(all(c("n_edges", "density", "mean_degree", "transitivity") %in%
                        names(feats$delaunay)))

    # k-NN sub-list
    expect_named(feats$knn, c("k3", "k5", "k7", "k10"))
    expect_true("clustering_coefficient" %in% names(feats$knn$k5))
})

test_that("compute_transitivity returns value in [0, 1]", {
    adj <- matrix(c(0,1,1, 1,0,1, 1,1,0), 3L, 3L)
    t   <- CellGraphR:::compute_transitivity(adj)
    expect_gte(t, 0)
    expect_lte(t, 1)
})

test_that("compute_mst_features handles single-point input", {
    feats <- CellGraphR:::compute_mst_features(matrix(c(1,2), 1L, 2L), "A")
    expect_equal(feats$total_length, 0)
})
