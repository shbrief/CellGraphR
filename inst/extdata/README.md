# Graph Features

This directory contains cell-level graph data and patient-level graph-derived
features extracted from TCGA whole-slide histopathology images.

## Files

### HDF5 Graph Files (`*_graph.h5`)

Per-slide cell graphs stored in HDF5 format with two tables:

- **`nodes`** — Cell-level node attributes (cell type, coordinates, etc.)
- **`edges`** — Edge list (`source`, `target`) representing spatial adjacency between cells

Example slides:
- `TCGA-EJ-5494-01Z-00-DX1.…_graph.h5` (~50 MB, ~1M cells)
- `TCGA-KK-A59Y-01Z-00-DX1.…_graph.h5` (~47 MB, ~1M cells)

### AnnData File (`*.h5ad.gz`)

- `TCGA-KK-A59Y-01Z-00-DX1.….h5ad.gz` — Compressed AnnData object for the KK-A59Y slide, suitable for use with `importGraphH5()`.

### Slide Image (`*.png`)

- `TCGA-KK-A59Y-01Z-00-DX1.….png` — Thumbnail of the KK-A59Y whole-slide image for visualization.

### Per-Sample Feature CSVs (`TCGA-HI-7168-*`)

- `TCGA-HI-7168-01Z-00-DX1.…_graph_features_summary.csv` — Cell-level graph feature summary for one slide.
- `TCGA-HI-7168-01Z-00-DX1.…_graph_smooth_embedding.csv` — Smooth embedding coordinates for cells in the HI-7168 slide.

### Cohort-Level Feature CSVs

| File | Description |
|---|---|
| `all_fullcohort_combined_expanded_patient_level_features_2.csv` | Patient-level summary features aggregated from cell graphs. Each row is one patient. |
| `example_graph_features.csv` | Small example graph feature table for package examples and tests. |
| `tcga_multiscale_graph_features_per_slide.csv` | Per-slide multiscale graph features across the TCGA cohort. |

#### Feature categories in the cohort CSVs

| Category | Examples |
|---|---|
| **Cell-type fractions** | `frac_neoplastic`, `frac_stromal`, `frac_inflammatory`, … |
| **Degree statistics** | `mean_degree`, `std_degree`, `max_degree`, per-cell-type mean degree |
| **Mixing indices** | `mean_mixing`, per-cell-type mixing index (homotypic vs. heterotypic neighbors) |
| **Edge-type fractions** | `edge_frac_neoplastic__stromal`, `edge_frac_inflammatory__neoplastic`, … |
| **Multi-scale features** | Degree/mixing at radii r20, r40, r80, r100 |
| **Spatial statistics** | Moran's I, clustering coefficients, Voronoi area, nearest-neighbor distances |
| **Tumor architecture** | Island count/size, fragmentation, invasion front ratio, eccentricity |
| **Immune microenvironment** | TIL margin fraction, immune exclusion score, immune clustering |
| **Graph topology** | Modularity, entropy, triangle counts, connected components |
