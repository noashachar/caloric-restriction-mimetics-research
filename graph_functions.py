import os
from parse_desc import parse_desc_file
from parse_genes import k_most_expressed_genes
from glob import glob
import networkx as nx
from networkx.algorithms import bipartite


def build__k_most_expressed__homogeneous_graph(dirpath, k):    
    category = os.path.basename(dirpath)  # "CRMs" or "Fasting" or "CR"
    
    B = nx.Graph(name=f'Gene-{category}')
    
    data_paths_pairs = zip(
        sorted(glob(f'{dirpath}/*.desc.txt')), 
        sorted(glob(f'{dirpath}/*.genes.txt'))
    )
    
    for desc_path, genes_path in data_paths_pairs:
        desc_dict = parse_desc_file(desc_path)
        
        B.add_node(desc_dict['Name'], category=category, **desc_dict)
        
        genes_arr = k_most_expressed_genes(genes_path, k)
        
        for _, (gene_symbol, logFC) in genes_arr.iterrows():
            B.add_node(gene_symbol, category='gene',  Name=gene_symbol, Organism=desc_dict['Organism'])
            B.add_edge(desc_dict['Name'], gene_symbol, logFC=logFC)
    
    assert bipartite.is_bipartite(B)
    return B


def without_lone_nodes(G):
    """
    Returns a subgraph with the vertices whose degree is at least 2.
    """
    while True:
        partialG_nodes = [node for node, degree in dict(G.degree()).items() if degree >= 2]
        if len(partialG_nodes) == len(G.nodes()):
            return G
        G = G.subgraph(partialG_nodes).copy()
        
        
def compose_graphs(g1, g2, g3):
    res = nx.compose(g1, g2)
    res = nx.compose(res, g3)
    return res
