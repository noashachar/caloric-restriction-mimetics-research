import re
import pickle
import mygene
import os.path
import pandas as pd
from glob import glob
from research_data import paths_pairs
from parse_desc import parse_desc_file


def remove_symbol_version(gene: str):
    """E.g. ENSG00000140455.11 -> ENSG00000140455"""
    return re.sub(r"\.\d+$", "", gene)


def aggregate_gene_symbols_and_query_mygene():
    print("Pre-Processing Gene Symbols...")

    genes_of = {"human": set(), "mouse": set()}

    for desc_path, genes_path in paths_pairs:
        df = pd.read_csv(genes_path, sep="\t")
        organism = parse_desc_file(desc_path)["Organism"]
        gene_symbols = {remove_symbol_version(symbol) for symbol in df.iloc[:, 0]}
        genes_of[organism].update(gene_symbols)

    mg = mygene.MyGeneInfo()

    def _normalize_gene_symbols(gene_symbols, species):
        scopes = "symbol,entrezgene,ensemblgene,refseq,mgi,hgnc,name,alias,uniprot,clone,other"
        query_results = mg.querymany(
            gene_symbols, scopes=scopes, fields="symbol", species=species
        )
        return {res["query"]: res["symbol"] for res in query_results if "symbol" in res}

    print("> Querying MyGene for Mouse genes...")
    normalized_symbols_mouse = _normalize_gene_symbols(genes_of["mouse"], "mouse")

    print("> Querying MyGene for Human genes...")
    normalized_symbols_human = _normalize_gene_symbols(genes_of["human"], "human")

    # Convert all symbols to lower case to minimize human/mouse query results disagreements

    normalized_symbols_mouse = {
        key.lower(): value.lower() for key, value in normalized_symbols_mouse.items()
    }

    normalized_symbols_human = {
        key.lower(): value.lower() for key, value in normalized_symbols_human.items()
    }

    disagreements = [
        human_k
        for human_k, human_v in normalized_symbols_human.items()
        if human_v != normalized_symbols_mouse.get(human_k, human_v)
    ]

    print("> Human & Mouse Symbols Map Disagreements:", len(disagreements))

    unifed_gene_symbol_map = {}

    for human_k, human_v in normalized_symbols_human.items():
        # if human and mouse agree on symbol name
        if human_v == normalized_symbols_mouse.get(human_k, human_v):
            unifed_gene_symbol_map[human_k] = human_v
        # if human and mouse don't agree on symbol name
        else:
            mouse_v = normalized_symbols_mouse[human_k]
            # pick the shorter one (e.g. fbxl9p vs. fbxl9)
            unifed_gene_symbol_map[human_k] = min([human_v, mouse_v], key=len)

    # add missing entries from mouse dict
    for mouse_k, mouse_v in normalized_symbols_mouse.items():
        if mouse_k not in unifed_gene_symbol_map:
            unifed_gene_symbol_map[mouse_k] = mouse_v

    print("> All Disagreements have been Settled.")

    return unifed_gene_symbol_map


def get_unifed_gene_symbol_map(cache_fname):
    if os.path.isfile(cache_fname):
        with open(cache_fname, "rb") as file:
            return pickle.load(file)

    unifed_gene_symbol_map = aggregate_gene_symbols_and_query_mygene()

    with open(cache_fname, "wb") as file:
        pickle.dump(unifed_gene_symbol_map, file)

    return unifed_gene_symbol_map


"""
a global singleton object
"""
unifed_gene_symbol_map = get_unifed_gene_symbol_map(cache_fname="gene_symbol_map.pkl")

known_gene_symbols = list(unifed_gene_symbol_map.values())


def normalize_symbol(gene: str):
    """
    # 1. remove gene symbol version number
    # 2. if it's in the dict, map to shorter common symbol
    """
    gene = gene.lower()
    gene = remove_symbol_version(gene)
    return unifed_gene_symbol_map.get(gene, gene)


def all_genes_in_file(genes_path: str):
    df = pd.read_csv(genes_path, sep="\t")

    # only keep rows that are significant
    df = df[df["P.Value"] <= 0.05]

    # normalize the gene symbols
    df["gene_symbol"] = df["gene_symbol"].apply(normalize_symbol)

    return df[["gene_symbol", "logFC"]]


def k_most_expressed_genes(genes_path: str, k: int):
    df = all_genes_in_file(genes_path)
    # the k rows with most extreme logFC value
    rows = df["logFC"].abs().nlargest(k).index
    return df.loc[rows].reset_index(drop=True)
