# CellGraphR <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**CellGraphR** is an R/Bioconductor package for importing, analysing, and
visualising cell-level spatial graphs derived from whole-slide histopathology
images. It bridges a Python cell-graph construction pipeline with the
Bioconductor analysis ecosystem, providing a unified workflow from raw graph
data to clinical association testing.

## Features

- **Import** pre-computed PyTables HDF5 cell-graph files into `igraph` objects
- **Construct** spatial graphs (Delaunay, k-NN, Gabriel, relative neighbourhood,
  MST) from cell coordinates
- **Extract** graph-topology features (transitivity, path length, motifs,
  spectral metrics) and spatial statistics (Moran's I, Ripley's K/L,
  neighbourhood composition, Hopkins clustering)
- **Store** patient-level feature matrices as `SummarizedExperiment` and
  integrate with `MultiAssayExperiment`
- **Visualise** spatial graphs overlaid on tissue images with `ggplot2`
- **Test** associations between MOFA latent factors and clinical variables

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("shbrief/CellGraphR")
```

### System requirements

CellGraphR uses **reticulate** to read HDF5 graph files produced by a Python
pipeline. The following Python packages are required:

- Python >= 3.8
- pandas >= 1.3
- tables >= 3.6 (PyTables)

## Quick start

```r
library(CellGraphR)

# --- 1. Import a cell graph from an HDF5 file ---
h5_path <- system.file("extdata", "TCGA-EJ-5494_graph.h5", package = "CellGraphR")
g <- importGraphH5(h5_path)

# --- 2. Build a spatial graph from nuclei coordinates ---
sg <- buildSpatialGraph(nuclei_list, graph_type = "delaunay")

# --- 3. Visualise the graph ---
plotGraph(g)

# --- 4. Extract graph features ---
topo_feats <- extract_graph_features(nuclei_list)
spatial_feats <- extract_advanced_spatial_features(nuclei_list)

# --- 5. Import features into a SummarizedExperiment ---
csv_path <- system.file("extdata", "example_graph_features.csv", package = "CellGraphR")
se <- importGraphFeatures(csv_path)

# --- 6. Test factor–clinical associations ---
results <- testFactorAssociations(factor_scores, clinical_df)
```

## Vignettes

Detailed walkthroughs are available at the
[package website](https://shbrief.github.io/CellGraphR):

| Vignette | Description |
|---|---|
| [Graph Visualization](https://shbrief.github.io/CellGraphR/articles/graph_visualization.html) | Plotting spatial graphs on tissue images |
| [Graph Analysis with igraph](https://shbrief.github.io/CellGraphR/articles/graph_igraph_analysis.html) | Extracting topology and spatial features |
| [Multimodal Analysis](https://shbrief.github.io/CellGraphR/articles/graph_multimodal_analysis.html) | Integrating graph features with multi-omics data |
| [EDA — PRAD Graphs](https://shbrief.github.io/CellGraphR/articles/EDA_prad_graph.html) | Exploratory analysis of prostate adenocarcinoma graphs |

## Contributing

Bug reports and feature requests are welcome on the
[GitHub issue tracker](https://github.com/shbrief/CellGraphR/issues).

## License

Artistic-2.0
