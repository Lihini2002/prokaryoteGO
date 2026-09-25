from pathlib import Path
import pandas as pd

# Locate src/ relative to script location
PREPROCESSING_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = PREPROCESSING_DIR.parent

# Configuration Paths
DATA_RAW_DIR = PROJECT_ROOT / "data" / "raw"
DATA_PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
DATAFILES_DIR = PROJECT_ROOT / "datafiles"

DATA_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# Input Files (Fixed: Removed trailing comma)
SELECTED_COMPLEX_PDB_CHAIN_SUMMARY_FILE = (
    DATA_PROCESSED_DIR / "list_of_selected_complex_chain_summary.tsv"
)
pdb_to_uniprot_file = DATAFILES_DIR / "pdb_chain_uniprot.tsv"

# Output Files
SELECTED_QUATERNARY_UNIPROT_MAPPING_FILE = (
    DATA_PROCESSED_DIR / "selected_quaternary_uniprot_mapping.tsv"
)
def map_chains_to_uniprot(
    selected_chains_summary_path: Path | str,
    sifts_path: Path | str,
    output_path: Path | str,
):
    summary_df = pd.read_csv(selected_chains_summary_path, sep="\t")

    chain_records = []
    for _, row in summary_df.iterrows():
        pdb_id = str(row["PDB_ID"]).strip().upper()
        clustered_chains = str(row["CLUSTERED_CHAINS"]).split(",")

        for chain_id in clustered_chains:
            chain_id = chain_id.strip()
            if chain_id and chain_id != "nan":
                chain_records.append({
                    "PARENT_PDB_ID": pdb_id,
                    "CHAIN_ID": chain_id,
                    "MATCH_CHAIN_ID": chain_id.upper(),  # Standardized for matching
                    "FULL_CHAIN_ID": f"{pdb_id}_{chain_id}",
                })

    chains_df = pd.DataFrame(chain_records)

    # 2. Read SIFTS file and normalize casing
    sifts_df = pd.read_csv(
        sifts_path,
        sep="\t",
        comment="#",
        usecols=["PDB", "CHAIN", "SP_PRIMARY"],
        dtype=str,
    ).dropna(subset=["SP_PRIMARY"])

    sifts_df["PARENT_PDB_ID"] = sifts_df["PDB"].str.upper().str.strip()
    sifts_df["MATCH_CHAIN_ID"] = sifts_df["CHAIN"].str.upper().str.strip()
    sifts_df.rename(columns={"SP_PRIMARY": "UNIPROT_ID"}, inplace=True)

    # 3. Merge on uppercase-normalized PDB and CHAIN IDs
    merged_df = pd.merge(
        chains_df,
        sifts_df[["PARENT_PDB_ID", "MATCH_CHAIN_ID", "UNIPROT_ID"]],
        on=["PARENT_PDB_ID", "MATCH_CHAIN_ID"],
        how="left",
    )

    # Clean up temporary matching columns
    final_df = merged_df.drop_duplicates(subset=["FULL_CHAIN_ID", "UNIPROT_ID"])
    final_df = final_df[
        ["PARENT_PDB_ID", "CHAIN_ID", "FULL_CHAIN_ID", "UNIPROT_ID"]
    ]

    # Save outputs
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    final_df.to_csv(output_path, sep="\t", index=False)
    mapped_count = final_df["UNIPROT_ID"].notna().sum()
    print(
        f"Mapped {mapped_count}/{len(final_df)} chains to UniProt IDs. Saved to '{output_path}'."
    )

    # Find chains that failed to map
   # Filter unmapped chains
    unmapped = final_df[final_df["UNIPROT_ID"].isna()]

    total_chains = len(final_df)
    unmapped_count = len(unmapped)
    unmapped_pct = (unmapped_count / total_chains * 100) if total_chains > 0 else 0
    unique_pdbs_affected = unmapped["PARENT_PDB_ID"].nunique()

    # Print summary statistics
    print(f"\n--- Unmapped Chains Summary ---")
    print(f"Total unmapped chains : {unmapped_count} / {total_chains} ({unmapped_pct:.2f}%)")
    print(f"Unique PDBs affected  : {unique_pdbs_affected}")
    print(f"\nUnmapped chains sample (first {min(100, unmapped_count)}):\n")
    print(unmapped[["PARENT_PDB_ID", "CHAIN_ID", "FULL_CHAIN_ID"]].head(100))

    

    return final_df

# def map_chains_to_uniprot(
#     selected_chains_summary_path: Path | str,
#     sifts_path: Path | str,
#     output_path: Path | str,
# ):
#     """Parses complex chain summary TSV, expands clustered chains, maps to UniProt via SIFTS,

#     and exports to TSV.
#     """
#     # 1. Read summary TSV and expand CLUSTERED_CHAINS (e.g., 'A,B' -> 'A', 'B')
#     summary_df = pd.read_csv(selected_chains_summary_path, sep="\t")

#     chain_records = []
#     for _, row in summary_df.iterrows():
#         pdb_id = str(row["PDB_ID"]).strip().upper()
#         clustered_chains = str(row["CLUSTERED_CHAINS"]).split(",")

#         for chain_id in clustered_chains:
#             chain_id = chain_id.strip()
#             if chain_id:
#                 chain_records.append({
#                     "PARENT_PDB_ID": pdb_id,
#                     "CHAIN_ID": chain_id,
#                     "FULL_CHAIN_ID": f"{pdb_id}_{chain_id}",
#                 })

#     chains_df = pd.DataFrame(chain_records)

#     # 2. Read SIFTS mapping file
#     sifts_df = pd.read_csv(
#         sifts_path,
#         sep="\t",
#         comment="#",
#         usecols=["PDB", "CHAIN", "SP_PRIMARY"],
#         dtype=str,
#     )

#     sifts_df["PDB"] = sifts_df["PDB"].str.upper().str.strip()
#     sifts_df["CHAIN"] = sifts_df["CHAIN"].str.strip()
#     sifts_df.rename(
#         columns={
#             "PDB": "PARENT_PDB_ID",
#             "CHAIN": "CHAIN_ID",
#             "SP_PRIMARY": "UNIPROT_ID",
#         },
#         inplace=True,
#     )

#     # 3. Left join target chains with SIFTS mappings
#     merged_df = pd.merge(
#         chains_df, sifts_df, on=["PARENT_PDB_ID", "CHAIN_ID"], how="left"
#     )

#     # Remove duplicates while preserving unmapped chains
#     final_df = merged_df.drop_duplicates(subset=["FULL_CHAIN_ID", "UNIPROT_ID"])

#     # Arrange final columns
#     final_df = final_df[
#         ["PARENT_PDB_ID", "CHAIN_ID", "FULL_CHAIN_ID", "UNIPROT_ID"]
#     ]

#     # 4. Save to output TSV
#     Path(output_path).parent.mkdir(parents=True, exist_ok=True)
#     final_df.to_csv(output_path, sep="\t", index=False)

#     mapped_count = final_df["UNIPROT_ID"].notna().sum()
#     print(
#         f"Mapped {mapped_count}/{len(final_df)} chains to UniProt IDs. Saved to '{output_path}'."
#     )

#     # Find chains that failed to map
#     unmapped = final_df[final_df["UNIPROT_ID"].isna()]
#     print(f"Unmapped chains sample:\n{unmapped.head(100)}")
    
#     return final_df


if __name__ == "__main__":
    map_chains_to_uniprot(
        selected_chains_summary_path=SELECTED_COMPLEX_PDB_CHAIN_SUMMARY_FILE,
        sifts_path=pdb_to_uniprot_file,
        output_path=SELECTED_QUATERNARY_UNIPROT_MAPPING_FILE,
    )


