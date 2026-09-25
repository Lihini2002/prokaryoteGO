import os
import time
import requests
from pathlib import Path

# --- Configuration & Constants ---
SRC_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SRC_DIR.parent

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RAW_DATA_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "data", "raw")

IDS_OUTPUT_PATH = os.path.join(RAW_DATA_DIR, "bacterial_pdb_ids.txt")
METADATA_OUTPUT_PATH = os.path.join(RAW_DATA_DIR, "bacterial_pdb_metadata.tsv")

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"


def fetch_bacterial_pdb_ids(page_size=10000) -> list[str]:
    """Queries RCSB Search API to fetch all bacterial protein PDB IDs via pagination."""
    query_payload = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.taxonomy_lineage.id",
                        "operator": "exact_match",
                        "value": "2"  # Bacteria superkingdom
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "entity_poly.rcsb_entity_polymer_type",
                        "operator": "exact_match",
                        "value": "Protein"
                    }
                }
            ]
        },
        "return_type": "entry",
        "request_options": {
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]
        }
    }

    all_ids = []
    start = 0
    print("Fetching bacterial PDB IDs...")

    while True:
        query_payload["request_options"]["paginate"] = {"start": start, "rows": page_size}
        response = requests.post(SEARCH_URL, json=query_payload)
        response.raise_for_status()
        
        batch = [hit["identifier"] for hit in response.json().get("result_set", [])]
        if not batch:
            break
            
        all_ids.extend(batch)
        print(f"  Retrieved {len(all_ids)} IDs...")
        start += page_size

    return all_ids


def fetch_metadata_batch(entry_ids: list[str]) -> list[dict]:
    """Sends a GraphQL request to fetch resolution, R-free, method, and chain IDs for a batch of PDB IDs."""
    query = """
    query ($ids: [String!]!) {
        entries(entry_ids: $ids) {
            rcsb_id
            exptl { method }
            refine { ls_R_factor_R_free }
            rcsb_entry_info { resolution_combined }
            polymer_entities {
                entity_poly {
                    rcsb_entity_polymer_type
                }
                rcsb_polymer_entity_container_identifiers {
                    auth_asym_ids
                }
            }
        }
    }
    """
    resp = requests.post(GRAPHQL_URL, json={"query": query, "variables": {"ids": entry_ids}})
    resp.raise_for_status()
    return resp.json().get("data", {}).get("entries", [])


def stream_metadata_to_file(pdb_ids: list[str], output_path: str, batch_size=500):
    """Fetches metadata in batches and streams results directly to a TSV file on disk."""
    print(f"\nFetching metadata and chain IDs in batches of {batch_size}...")
    
    with open(output_path, "w") as f:
        # Write TSV Header with chains and full_chain_ids
        f.write("pdb_id\tchains\tfull_chain_ids\tmethod\tresolution\tr_free\n")

        for i in range(0, len(pdb_ids), batch_size):
            batch = pdb_ids[i:i + batch_size]
            
            # Retry mechanism
            entries = []
            for attempt in range(2):
                try:
                    entries = fetch_metadata_batch(batch)
                    break
                except requests.exceptions.RequestException as e:
                    if attempt == 1:
                        print(f"  Batch {i}-{i+batch_size} failed. Skipping: {e}")
                    else:
                        time.sleep(2)

            # Write batch directly to file
            for entry in entries:
                if not entry:
                    continue
                
                rcsb_id = entry.get("rcsb_id")
                methods = entry.get("exptl") or []
                method = methods[0]["method"] if methods else "N/A"
                
                refine = entry.get("refine") or []
                r_free = refine[0].get("ls_R_factor_R_free") if refine else "N/A"
                
                entry_info = entry.get("rcsb_entry_info") or {}
                res = entry_info.get("resolution_combined")
                resolution = res[0] if isinstance(res, list) and res else res or "N/A"

                # Extract individual protein chain IDs
                chains = []
                polymer_entities = entry.get("polymer_entities") or []
                for entity in polymer_entities:
                    poly_type = (entity.get("entity_poly") or {}).get("rcsb_entity_polymer_type")
                    if poly_type == "Protein":
                        identifiers = entity.get("rcsb_polymer_entity_container_identifiers") or {}
                        auth_asym_ids = identifiers.get("auth_asym_ids") or []
                        chains.extend(auth_asym_ids)

                chains_sorted = sorted(set(chains))
                chains_str = ",".join(chains_sorted) if chains_sorted else "N/A"
                full_chains_str = ",".join([f"{rcsb_id}_{c}" for c in chains_sorted]) if chains_sorted else "N/A"

                f.write(f"{rcsb_id}\t{chains_str}\t{full_chains_str}\t{method}\t{resolution}\t{r_free}\n")

            f.flush()
            print(f"  Processed {min(i + batch_size, len(pdb_ids))} / {len(pdb_ids)} entries")
            time.sleep(0.1)


def main():
    os.makedirs(RAW_DATA_DIR, exist_ok=True)

    # 1. Fetch & Save PDB IDs
    pdb_ids = fetch_bacterial_pdb_ids()
    with open(IDS_OUTPUT_PATH, "w") as f:
        f.write("\n".join(pdb_ids))
    print(f"Saved {len(pdb_ids)} PDB IDs to: {IDS_OUTPUT_PATH}")

    # 2. Fetch & Stream Metadata
    stream_metadata_to_file(pdb_ids, METADATA_OUTPUT_PATH)
    print(f"\nPipeline complete! Metadata saved to: {METADATA_OUTPUT_PATH}")


if __name__ == "__main__":
    main()