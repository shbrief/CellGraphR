import pandas as pd
import networkx as nx

graph = '~/Packages/CellGraphR/inst/extdata/TCGA-EJ-5494-01Z-00-DX1.f466c17c-4671-45e0-aca9-0ceeec2c0280_graph.h5'
nodes_df = pd.read_hdf(graph, key='nodes')
edges_df = pd.read_hdf(graph, key='edges')
G = nx.from_pandas_edgelist(edges_df, source='source', target='target', create_using=nx.Graph)

import anndata as ad

file_path = "/Users/sehyunoh/Packages/CellGraphR/inst/extdata/TCGA-KK-A59Y-01Z-00-DX1.0DB9AF4D-4879-4295-A321-CBEE0FE6B2A0.h5ad"
adata = ad.read_h5ad(file_path)
coords = adata.obsm["spatial"]
print(coords)
