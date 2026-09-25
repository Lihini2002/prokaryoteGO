
"""
This module contains functions to obtain the pdb ids of bacterial proteins and their associated metadata (resolution, R-free, and method) from the RCSB PDB API. 
The metadata is streamed to a file in batches to avoid memory issues with large datasets. 
Filter based on the resolution  
get the proper sequences frmo pdbseqres file 
run CD-HIT to remove redundancy and store the non-redundant sequences in a fasta file. 

After this , do chain wise analysis. 
A complex signature of the protein is generated with cluster ids of all the sub chains
if the complex signature is similar two proteins are considered similar.Only one of them needs to be retained. 
If else both are retained.
"""

# refactor this with data classes. 


# The bacterialquality id are found in data/raw/bacterial_pdb_ids.txt
# bacterial quality metadata is found in data/raw/bacterial_pdb_metadata.tsv


#filter based on resolution and R-free thresholds. 
    #RESOLUTION_THRESHOLD = 2.5  # Angstroms
    #R_FREE_THRESHOLD = 0.25  # R-free value


#FUNCTION to get the sequences of each of these from pdb seqres file and store them in a fasta file. 
#function to run CD-HIT to remove redundancy and store the non-redundant sequences in a fasta file. 

#make a complex signature of the protein with cluster ids of all the sub chains 
    #make complex signature of the proteins 
    #taking the pdb id from the non-redundant fasta file and 
    #match with the pdb id from pdbchaintaxonomy file. and match it 
        #getting the cluster ids of all the sub chains and storing them in a list. 

        #if two entries have the same complex signature then they are considered similar and only one of them is retained.
        #if the complex signature is different then both are retained. 


import os
import subprocess
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from pathlib import Path

# Locate src/ relative to this script's location, then get the project root (parent directory)
PREPROCESSING_DIR = Path(__file__).resolve().parent  # project/preprocessing/
PROJECT_ROOT = PREPROCESSING_DIR.parent

# --- Configuration Paths & Thresholds ---
DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
DATAFILES_DIR = PROJECT_ROOT / "datafiles"

# Ensure output directories exist
DATA_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Input Files
PDB_IDS_FILE = DATA_RAW_DIR / "bacterial_pdb_ids.txt"
METADATA_FILE = DATA_RAW_DIR / "bacterial_pdb_metadata.tsv"
PDB_SEQRES_FILE = DATAFILES_DIR / "pdb_seqres.txt"
PDB_TO_UNIPROT = DATAFILES_DIR / "pdb_chain_uniprot.tsv"

# Output Files
QUALITY_FASTA = DATA_PROCESSED_DIR / "quality_bacterial.fasta"
CDHIT_OUTPUT_FASTA = DATA_PROCESSED_DIR / "nr_bacterial.fasta"
CDHIT_CLSTR_FILE = DATA_PROCESSED_DIR / "nr_bacterial.fasta.clstr"
DEREPLICATED_FASTA = DATA_PROCESSED_DIR / "final_unique_complexes.fasta"
SELECTED_COMPLEX_PDB_CHAIN_SUMMARY_FILE = DATA_PROCESSED_DIR / "list_of_selected_complex_chain_summary.tsv"


RESOLUTION_THRESHOLD = 2.5  # Angstroms
R_FREE_THRESHOLD = 0.25     # Max acceptable R-free value

import pandas as pd




# filter the pdb ids based on resolution and R-free thresholds. 
# if you think there is insufficient data , skip this function. 

def filter_pdb_by_quality(ids_path: str, metadata_path: str) -> set:
    """Reads PDB IDs and filters them based on Resolution, R-free, and experimental method."""
    
    with open(ids_path, "r") as f:
        target_ids = {line.strip().upper() for line in f if line.strip()}

    df = pd.read_csv(metadata_path, sep="\t")
    df.columns = df.columns.str.lower().str.strip()
    df["pdb_id"] = df["pdb_id"].str.upper().str.strip()

    df = df[df["pdb_id"].isin(target_ids)].copy()

    # Convert numeric columns safely (turns 'N/A' or invalid values into NaN)
    df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
    df["r_free"] = pd.to_numeric(df["r_free"], errors="coerce")

    xray_pass = (
        (df["method"].str.contains("X-RAY", na=False, case=False)) & 
        (df["resolution"] <= 2.5) & 
        (df["r_free"].isna() | (df["r_free"] <= 0.25))
    )

    cryo_pass = (
        (df["method"].str.contains("ELECTRON MICROSCOPY", na=False, case=False)) & 
        (df["resolution"] <= 2.5)
    )

    nmr_pass = df["method"].str.contains("NMR|NUCLEAR MAGNETIC RESONANCE", na=False, case=False)

    filtered_df = df[xray_pass | cryo_pass | nmr_pass]

    return set(filtered_df["pdb_id"].unique())

# 2. Sequence Extraction
def extract_sequences(filtered_pdb_ids: set, seqres_path: str, output_fasta: str):
    """Parses local pdb_seqres file and saves sequences for filtered PDB IDs."""

    records_to_save = []

    # get the pdb id from the pdbseqres and getits sequence.
    # Iterate through PDB SEQRES file
    for record in SeqIO.parse(seqres_path, "fasta"):
        # Header format typically: 1abc_A mol:protein length:123 ...
       
        pdb_chain = record.id.split()[0]

        # extract the four letter pdb id from the chain id and convert to uppercase
        pdb_id = pdb_chain.split("_")[0].upper()


    # match it with the pdb id from the filtered pdb ids. 
        if pdb_id in filtered_pdb_ids:
            record.id = pdb_chain  # Set standard PDB_CHAIN identifier
            records_to_save.append(record)

    # writes all extracted pdbs into an output fasta file. 
    os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
    SeqIO.write(records_to_save, output_fasta, "fasta")
    print(f"Extracted {len(records_to_save)} chain sequences to {output_fasta}")


# 3. CD-HIT Execution
def run_cd_hit(input_fasta: str, output_fasta: str, sequence_identity: float = 0.9, kmer_size: int = 5):
    """Executes CD-HIT to cluster sequences at a given identity threshold."""
    cmd = [
        "cd-hit",
        "-i", input_fasta,
        "-o", output_fasta,
        "-c", str(sequence_identity),
        "-n", str(kmer_size),  # Word size for identity
        "-d", "0"
    ]
    subprocess.run(cmd, check=True)
    print(f"CD-HIT clustering complete. Output saved to {output_fasta}")


# 4. Complex Signature Generation & Dereplication
# this function gives individual chain ids and their cluster ids.
def parse_cdhit_clusters(clstr_path: str) -> dict:
    """Parses .clstr file to map each chain ID (e.g., '1ABC_A') to its CD-HIT cluster ID."""
    chain_to_cluster = {}
    current_cluster = None

    with open(clstr_path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                current_cluster = line.split()[1]  # Extract cluster index
            else:
                # Line format: 0 123aa, >1ABC_A... at 98.5%
                chain_id = line.split(">")[1].split("...")[0]
                chain_to_cluster[chain_id] = current_cluster

    return chain_to_cluster


# this function thinks that the pdb_ids have chains attached to them. 
def build_complex_signatures(chain_to_cluster: dict) -> dict:
    """Generates a sorted complex signature tuple for each PDB entry."""
    pdb_to_chains = defaultdict(list)
    for chain_id, cluster_id in chain_to_cluster.items():
        pdb_id = chain_id.split("_")[0].upper()
        pdb_to_chains[pdb_id].append(cluster_id)

    # Sort cluster IDs to make signature invariant to chain order
    complex_signatures = {
        pdb_id: tuple(sorted(clusters))
        for pdb_id, clusters in pdb_to_chains.items()
    }
    return complex_signatures

from collections import defaultdict
import pandas as pd


def dereplicate_complexes(
    complex_signatures: dict,
    chain_to_cluster: dict,
    sifts_path: str,
) -> set:
    """Retains one representative PDB ID per distinct complex signature, prioritizing

    candidates whose chains have valid UniProt mappings in SIFTS.
    """
    # 1. Read SIFTS file and extract valid (mapped) full chain IDs
    sifts_df = pd.read_csv(
        sifts_path,
        sep="\t",
        comment="#",
        usecols=["PDB", "CHAIN", "SP_PRIMARY"],
        dtype=str,
    ).dropna(subset=["SP_PRIMARY"])

    # Create lookup set of mapped chains (e.g., '1ABC_A')
    # Note: PDB is uppercase, CHAIN preserves case sensitivity
    # Fix: Convert both PDB and CHAIN to uppercase
    mapped_chains = set(
        sifts_df["PDB"].str.upper().str.strip()
        + "_"
        + sifts_df["CHAIN"].str.upper().str.strip()  # <--- Fixed
    )

    # 2. Map each PDB ID to its constituent chains in this dataset
    pdb_to_chains = defaultdict(list)
    for chain_id in chain_to_cluster.keys():
        pdb_id = chain_id.split("_")[0].upper()
        pdb_to_chains[pdb_id].append(chain_id)

    # 3. Group candidate PDB IDs by their complex signature
    signature_to_pdbs = defaultdict(list)
    for pdb_id, signature in complex_signatures.items():
        signature_to_pdbs[signature].append(pdb_id)

    retained_pdb_ids = set()

    # 4. For each signature group, choose the candidate PDB with the best UniProt coverage
    for signature, candidate_pdbs in signature_to_pdbs.items():
        best_pdb = None
        best_coverage = -1.0

        for pdb_id in candidate_pdbs:
            chains = pdb_to_chains.get(pdb_id, [])
            if not chains:
                continue

            # Calculate fraction of chains with valid UniProt IDs
            mapped_count = sum(1 for chain in chains if chain in mapped_chains)
            coverage = mapped_count / len(chains)

            # Retain the candidate with higher UniProt coverage
            if coverage > best_coverage:
                best_coverage = coverage
                best_pdb = pdb_id

            # Early exit if we find a candidate where 100% of chains are mapped
            if best_coverage == 1.0:
                break

        if best_pdb:
            retained_pdb_ids.add(best_pdb)

    print(
        f"Total signature groups: {len(signature_to_pdbs)} | Unique complexes retained: {len(retained_pdb_ids)}"
    )
    return retained_pdb_ids

# def dereplicate_complexes(complex_signatures: dict) -> set:
#     """Retains only one representative PDB ID per distinct complex signature."""
#     seen_signatures = set()
#     retained_pdb_ids = set()

#     for pdb_id, signature in complex_signatures.items():
#         if signature not in seen_signatures:
#             seen_signatures.add(signature)
#             retained_pdb_ids.add(pdb_id)

#     print(f"Total complexes: {len(complex_signatures)} | Unique complexes: {len(retained_pdb_ids)}")
#     return retained_pdb_ids

# def extract_pdb_chains(
#     fasta_path: str = DEREPLICATED_FASTA,
#     output_path: str = SELECTED_PDB_IDS_FILE,
# ):
#     ...
#     """Parses a FASTA file and extracts 'PDBID_CHAINID' identifiers into a text file."""
#     pdb_chains = []

#     with open(fasta_path, "r", encoding="utf-8") as infile:
#         for line in infile:
#             if line.startswith(">"):
#                 # Extract header ID (e.g., '10ac_A' from '>10ac_A mol:protein...')
#                 header_id = line.strip()[1:].split()[0]
#                 pdb_chains.append(header_id)

#     # Write the extracted IDs to the output file
#     with open(output_path, "w", encoding="utf-8") as outfile:
#         for chain_id in pdb_chains:
#             outfile.write(f"{chain_id}\n")

#     print(
#         f"Successfully extracted {len(pdb_chains)} PDB chain IDs to '{output_path}'."
#     )

from collections import defaultdict
import pandas as pd


def extract_pdb_chains_with_metadata(
    fasta_path: str = DEREPLICATED_FASTA,
    metadata_path: str = METADATA_FILE,
    output_path: str = SELECTED_COMPLEX_PDB_CHAIN_SUMMARY_FILE,
):
    """Parses dereplicated FASTA and matches with metadata to export:

    - Parent PDB ID
    - All original physical chains (from metadata)
    - Chains retained in the clustering FASTA
    """
    # 1. Extract clustered chains per parent PDB from FASTA
    clustered_chains = defaultdict(set)

    with open(fasta_path, "r", encoding="utf-8") as infile:
        for line in infile:
            if line.startswith(">"):
                full_id = line.strip()[1:].split()[0]  # e.g., '9ZDI_A'
                if "_" in full_id:
                    pdb_id, chain_id = full_id.split("_", 1)
                    clustered_chains[pdb_id.upper()].add(chain_id)

    # 2. Read metadata TSV file
    meta_df = pd.read_csv(metadata_path, sep="\t")
    meta_df.columns = meta_df.columns.str.lower().str.strip()
    meta_df["pdb_id"] = meta_df["pdb_id"].str.upper().str.strip()

    # Create fast lookup dict: { '9ZDI': 'A,B' }
    metadata_chains_map = meta_df.set_index("pdb_id")["chains"].to_dict()

    # 3. Build summary records
    records = []
    for pdb_id, retained_set in clustered_chains.items():
        retained_chains_str = ",".join(sorted(retained_set))

        # Retrieve original chains from metadata
        raw_meta_chains = metadata_chains_map.get(pdb_id, "")
        if pd.isna(raw_meta_chains) or not str(raw_meta_chains).strip():
            all_meta_chains_str = "N/A"
        else:
            all_meta_chains_str = ",".join(
                sorted(
                    [
                        c.strip()
                        for c in str(raw_meta_chains).split(",")
                        if c.strip()
                    ]
                )
            )

        records.append({
            "PDB_ID": pdb_id,
            "METADATA_CHAINS": all_meta_chains_str,
            "CLUSTERED_CHAINS": retained_chains_str,
            "METADATA_CHAIN_COUNT": len(all_meta_chains_str.split(","))
            if all_meta_chains_str != "N/A"
            else len(retained_set),
            "CLUSTERED_CHAIN_COUNT": len(retained_set),
        })

    # 4. Save to TSV file
    output_df = pd.DataFrame(records)
    output_df.to_csv(output_path, sep="\t", index=False)

    print(
        f"Successfully written complex chain summary for {len(output_df)} PDBs to '{output_path}'."
    )


# 5. Pipeline Execution
def run_pipeline():
    # Step 1: Filter by Resolution & R-free
    filtered_ids = filter_pdb_by_quality(PDB_IDS_FILE, METADATA_FILE)
    
    # Step 2: Extract sequences from seqres
    extract_sequences(filtered_ids, PDB_SEQRES_FILE, QUALITY_FASTA)
    
    # Step 3: Run CD-HIT sequence clustering
    run_cd_hit(QUALITY_FASTA, CDHIT_OUTPUT_FASTA, sequence_identity=0.90)
    
    # Step 4: Complex analysis based on sub-chain signatures
    chain_to_cluster = parse_cdhit_clusters(CDHIT_CLSTR_FILE)
    complex_signatures = build_complex_signatures(chain_to_cluster)
    unique_pdb_ids = dereplicate_complexes(complex_signatures ,chain_to_cluster,  PDB_TO_UNIPROT)

    # Step 5: Filter CD-HIT output FASTA to retain only unique complexes
    # Step 5 (Fixed): Filter the FULL sequence file for retained parent PDB IDs
    final_records = [
        rec for rec in SeqIO.parse(QUALITY_FASTA, "fasta")
        if rec.id.split("_")[0].upper() in unique_pdb_ids
    ]
    SeqIO.write(final_records, DEREPLICATED_FASTA, "fasta")
    print(f"Final dataset exported to {DEREPLICATED_FASTA}")


# write a test case to check if its both the same. 
    extract_pdb_chains_with_metadata(
        fasta_path=DEREPLICATED_FASTA,
        metadata_path=METADATA_FILE,
        output_path=SELECTED_COMPLEX_PDB_CHAIN_SUMMARY_FILE
    )

    


if __name__ == "__main__":
    run_pipeline()