import pandas as pd
import networkx as nx

graph = '~/Packages/HistoImageR/data/graph_features/TCGA-EJ-5494-01Z-00-DX1.f466c17c-4671-45e0-aca9-0ceeec2c0280_graph.h5'
nodes_df = pd.read_hdf(graph, key='nodes')
edges_df = pd.read_hdf(graph, key='edges')
G = nx.from_pandas_edgelist(edges_df, source='source', target='target', create_using=nx.Graph)
