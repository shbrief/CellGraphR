# Graph Features

This directory contains cell-level graph data and patient-level graph-derived features extracted from TCGA whole-slide histopathology images.

## Files

### HDF5 Graph Files (`*_graph.h5`)

Per-slide cell graphs stored in HDF5 format with two tables:

- **`nodes`** — Cell-level node attributes (cell type, coordinates, etc.)
- **`edges`** — Edge list (`source`, `target`) representing spatial adjacency between cells

These can be loaded into a NetworkX graph as shown in `load_graph_sample.py`.

Example slides:
- `TCGA-EJ-5494-01Z-00-DX1.…_graph.h5` (~50 MB, ~1M cells)
- `TCGA-KK-A59Y-01Z-00-DX1.…_graph.h5` (~47 MB, ~1M cells)

### `all_fullcohort_combined_expanded_patient_level_features_2.csv`

Patient-level summary features aggregated from the cell graphs. Each row is one patient. Key feature groups include:

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

### `load_graph_sample.py`

Helper script demonstrating how to load an HDF5 graph file into a pandas DataFrame and NetworkX graph.
