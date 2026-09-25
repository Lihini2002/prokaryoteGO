import gzip
from pathlib import Path
import networkx as nx
import obonet
import pandas as pd


# 1. Load GO Graph from OBO file
def load_go_graph(obo_file_path: str) -> nx.MultiDiGraph:
    """Parses go.obo into a NetworkX graph."""
    with open(obo_file_path, "r") as f:
        graph = obonet.read_obo(f)
    print(f"Loaded GO Graph with {len(graph)} nodes.")
    return graph


# 2. Extract Terms by Aspect (BP, MF, CC)
def split_go_aspects(graph: nx.MultiDiGraph):
    """Splits terms into Biological Process (BP), Molecular Function (MF), and

    Cellular Component (CC).
    """
    aspects = {"BP": set(), "MF": set(), "CC": set()}

    for node, data in graph.nodes(data=True):
        namespace = data.get("namespace")
        if namespace == "biological_process":
            aspects["BP"].add(node)
        elif namespace == "molecular_function":
            aspects["MF"].add(node)
        elif namespace == "cellular_component":
            aspects["CC"].add(node)

    print(
        f"Terms per aspect -> BP: {len(aspects['BP'])}, MF: {len(aspects['MF'])}, CC: {len(aspects['CC'])}"
    )
    return aspects


# 3. Propagate Annotations (True Path Rule)
def propagate_go_terms(
    direct_go_terms: set, graph: nx.MultiDiGraph
) -> set[str]:
    """Given a set of direct GO term annotations, returns the set including all

    ancestor GO terms up to the root.
    """
    propagated_terms = set(direct_go_terms)

    for term in direct_go_terms:
        if term in graph:
            # nx.ancestors finds all nodes reachable by following directed edges upwards
            ancestors = nx.ancestors(graph, term)
            propagated_terms.update(ancestors)

    return propagated_terms


# 4. Load SIFTS GO Mapping & Build Target Label Dict
def process_sifts_annotations(
    sifts_file: str, graph: nx.MultiDiGraph, target_uniprot_ids: set
) -> dict[str, dict[str, set[str]]]:
    """Processes SIFTS mapping (e.g., pdb2uniprot / pdb_chain_go.csv), extracts

    direct annotations, and applies True Path Rule propagation.
    """
    # SIFTS CSV/TSV format usually has: PDB, CHAIN, SP_PRIMARY (UniProt), GO_ID
    df = pd.read_csv(sifts_file, comment="#", low_memory=False)

    # Standardize columns (adjust column names based on your SIFTS file version)
    uniprot_col = "SP_PRIMARY" if "SP_PRIMARY" in df.columns else "UniProt"
    go_col = "GO_ID" if "GO_ID" in df.columns else "GO"

    df = df.dropna(subset=[uniprot_col, go_col])
    df = df[df[uniprot_col].isin(target_uniprot_ids)]

    # Group direct annotations by protein
    direct_annotations = df.groupby(uniprot_col)[go_col].apply(set).to_dict()

    # Apply True Path Rule per protein
    propagated_annotations = {}
    aspects_keys = split_go_aspects(graph)

    for uniprot_id, direct_terms in direct_annotations.items():
        all_terms = propagate_go_terms(direct_terms, graph)

        # Separate propagated terms into BP, MF, and CC sets
        propagated_annotations[uniprot_id] = {
            "BP": all_terms.intersection(aspects_keys["BP"]),
            "MF": all_terms.intersection(aspects_keys["MF"]),
            "CC": all_terms.intersection(aspects_keys["CC"]),
        }

    return propagated_annotations


# --- Example Pipeline Usage ---
if __name__ == "__main__":
    OBO_FILE = "go.obo"  # Download from http://current.geneontology.org/ontology/go.obo
    SIFTS_FILE = "pdb_chain_go.csv"

    # 1. Load DAG
    go_graph = load_go_graph(OBO_FILE)

    # 2. Example protein set
    my_proteins = {"P0A6F5", "P0A6L2", "P0A725"}

    # 3. Get fully propagated labels for all 3 GO types
    protein_labels = process_sifts_annotations(
        SIFTS_FILE, go_graph, my_proteins
    )

    # Inspect propagated GO terms for one protein
    sample_id = "P0A6F5"
    if sample_id in protein_labels:
        print(
            f"\nSample output for {sample_id}:"
        )
        print(f"BP terms count: {len(protein_labels[sample_id]['BP'])}")
        print(f"MF terms count: {len(protein_labels[sample_id]['MF'])}")
        print(f"CC terms count: {len(protein_labels[sample_id]['CC'])}")