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
SELECTED_PDB_TO_UNIPROT_MAPPING_FILE = (
    DATA_PROCESSED_DIR / "selected_quaternary_uniprot_mapping.tsv"
)

# Output Files
SELECTED_PDB_UNIPROT_STING_FILE = (
    DATA_PROCESSED_DIR / "selected_string_networks.tsv"
)


import io
from pathlib import Path
import time
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util import Retry


def create_robust_session() -> requests.Session:
    """Creates a requests session with exponential backoff retry logic."""
    session = requests.Session()
    retries = Retry(
        total=5,
        backoff_factor=1,
        status_forcelist=[429, 500, 502, 503, 504],
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retries)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


def map_uniprot_to_string(
    input_file: str | Path,
    output_file: str | Path,
    batch_size: int = 100,  # Reduced from 500 to prevent server timeouts
    request_timeout: int = 60,  # Increased timeout allowance
    delay_between_batches: float = 0.2,  # Rate limiting safety delay
) -> pd.DataFrame:
    """Reads a PDB-to-UniProt TSV file, filters out rows lacking a UniProt ID,

    queries STRING DB REST API to map UniProt IDs to STRING IDs, and saves the
    result.
    """
    # 1. Read input file
    df = pd.read_csv(input_file, sep=r"\s+|\t+", engine="python")
    df.columns = df.columns.str.strip().str.upper()

    if "UNIPROT_ID" not in df.columns:
        raise KeyError(
            f"Column 'UNIPROT_ID' not found in input. Found columns: {list(df.columns)}"
        )

    # 2. Filter out invalid/empty UniProt IDs ("only necessary ones")
    valid_mask = (
        df["UNIPROT_ID"].notna()
        & (df["UNIPROT_ID"].astype(str).str.strip() != "")
        & (df["UNIPROT_ID"].astype(str).str.upper() != "NAN")
    )
    filtered_df = df[valid_mask].copy()
    filtered_df["UNIPROT_ID"] = (
        filtered_df["UNIPROT_ID"].astype(str).str.strip()
    )

    print(
        f"Filtered dataset: Retained {len(filtered_df)}/{len(df)} rows with valid UniProt IDs."
    )

    # 3. Extract unique UniProt IDs to prevent redundant API calls
    unique_uniprots = filtered_df["UNIPROT_ID"].unique().tolist()
    uniprot_to_string_map = {}

    print(
        f"Fetching STRING IDs for {len(unique_uniprots)} unique UniProt IDs in batches of {batch_size}..."
    )

    session = create_robust_session()
    string_api_url = "https://string-db.org/api/tsv/get_string_ids"

    # 4. Batch query STRING DB REST API
    for i in range(0, len(unique_uniprots), batch_size):
        batch = unique_uniprots[i : i + batch_size]

        payload = {
            "identifiers": "\n".join(batch),
            "echo_query": "1",  # Keeps the query UniProt ID in response for exact matching
            "limit": "1",  # Return top match per UniProt ID
        }

        try:
            response = session.post(
                string_api_url, data=payload, timeout=request_timeout
            )
            response.raise_for_status()

            if response.text.strip():
                # Parse returned TSV string into DataFrame
                res_df = pd.read_csv(io.StringIO(response.text), sep="\t")

                # Map queryItem (UniProt) -> stringId
                if (
                    "queryItem" in res_df.columns
                    and "stringId" in res_df.columns
                ):
                    batch_map = dict(
                        zip(res_df["queryItem"], res_df["stringId"])
                    )
                    uniprot_to_string_map.update(batch_map)

        except requests.exceptions.RequestException as e:
            print(
                f"⚠️ API Warning: Failed batch {i}-{i+len(batch)} after retries: {e}"
            )

        # Pause slightly to avoid triggering STRING API rate limits
        time.sleep(delay_between_batches)

    # 5. Map STRING_ID column back to filtered DataFrame
    filtered_df["STRING_ID"] = filtered_df["UNIPROT_ID"].map(
        uniprot_to_string_map
    )

    # 6. Save output TSV
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    filtered_df.to_csv(output_file, sep="\t", index=False)

    mapped_count = filtered_df["STRING_ID"].notna().sum()
    print(
        f"✅ Mapped {mapped_count}/{len(filtered_df)} rows to STRING IDs. Saved to '{output_file}'."
    )

    return filtered_df


# --- Example Usage ---
if __name__ == "__main__":
    INPUT_FILE = SELECTED_PDB_TO_UNIPROT_MAPPING_FILE
    OUTPUT_FILE = SELECTED_PDB_UNIPROT_STING_FILE 

    mapped_df = map_uniprot_to_string(INPUT_FILE, OUTPUT_FILE)
    print("\nSample Output:")
    print(mapped_df.head(10))