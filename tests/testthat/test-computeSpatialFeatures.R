test_that("extract_advanced_spatial_features returns expected structure", {
    set.seed(1L)
    n <- 30L
    nuclei_list <- lapply(seq_len(n), function(i) {
        list(centroid = runif(2L, 0, 100),
             type     = sample(c("neoplastic", "stromal", "immune"), 1L))
    })

    feats <- extract_advanced_spatial_features(nuclei_list, image_dims = c(100, 100))

    expect_type(feats, "list")
    expect_named(feats,
                 c("spatial_autocorr", "multiscale", "neighborhood",
                   "clustering", "ripley_k", "nearest_neighbor"))
})

test_that("Moran's I is within plausible range", {
    set.seed(7L)
    n <- 50L
    nuclei_list <- lapply(seq_len(n), function(i) {
        list(centroid = runif(2L, 0, 500),
             type     = sample(c("A", "B"), 1L))
    })
    feats <- extract_advanced_spatial_features(nuclei_list, image_dims = c(500, 500))
    mi    <- feats$spatial_autocorr$morans_i
    expect_true(is.numeric(mi))
    expect_true(mi >= -1.5 && mi <= 1.5)
})

test_that("nearest_neighbor_analysis returns clark_evans_index", {
    set.seed(3L)
    n <- 20L
    nuclei_list <- lapply(seq_len(n), function(i) {
        list(centroid = runif(2L, 0, 200),
             type     = "A")
    })
    feats <- extract_advanced_spatial_features(nuclei_list, image_dims = c(200, 200))
    expect_true("clark_evans_index" %in% names(feats$nearest_neighbor))
    expect_true(is.numeric(feats$nearest_neighbor$clark_evans_index))
})
