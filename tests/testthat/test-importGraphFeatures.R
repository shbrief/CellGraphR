test_that("importGraphFeatures returns a SummarizedExperiment", {
    csv <- system.file("extdata", "example_graph_features.csv",
                       package = "CellGraphR")
    skip_if(!file.exists(csv), "example_graph_features.csv not found")

    se <- importGraphFeatures(csv)

    expect_s4_class(se, "SummarizedExperiment")
    expect_true("graph_features" %in% SummarizedExperiment::assayNames(se))
    expect_true("patient_id" %in% colnames(SummarizedExperiment::colData(se)))
    expect_true("feature_group" %in% colnames(SummarizedExperiment::rowData(se)))
})

test_that("importGraphFeatures dimensions are features x patients", {
    csv <- system.file("extdata", "example_graph_features.csv",
                       package = "CellGraphR")
    skip_if(!file.exists(csv), "example_graph_features.csv not found")

    se <- importGraphFeatures(csv)
    mat <- SummarizedExperiment::assay(se, "graph_features")

    # rows = features, cols = patients
    expect_gt(nrow(mat), 1L)
    expect_gt(ncol(mat), 1L)
    expect_true(is.numeric(mat))
})

test_that(".classify_graph_features returns correct group labels", {
    feats <- c("n_cells", "frac_neoplastic", "edge_frac_A__B",
               "mean_deg_neoplastic_r20", "morans_I_neoplastic",
               "TIL_at_margin_frac", "voronoi_area_mean",
               "total_triangles", "some_unknown_feature_dispersion")
    grps <- CellGraphR:::.classify_graph_features(feats)

    expect_equal(grps[1L], "global_topology")
    expect_equal(grps[2L], "cell_type_composition")
    expect_equal(grps[3L], "cell_type_mixing")
    expect_equal(grps[4L], "neighborhood_radius")
    expect_equal(grps[5L], "spatial_autocorrelation")
    expect_equal(grps[6L], "immune_spatial")
    expect_equal(grps[7L], "voronoi")
    expect_equal(grps[8L], "network_motif")
    expect_equal(grps[9L], "spatial_distribution")
})
