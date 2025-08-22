import pandas as pd
from pathlib import Path

base = Path("../../data/raw/data_for_tree/interpro")
out = Path("../../data/interm/merged")

# --- Collect CSVs (IPR*/**/*.csv + optional scaloptoxins.csv) ---
csv_paths = sorted(base.glob("IPR*/**/*.csv"))
scalop = base / "scaloptoxins.csv"
if scalop.exists():
    csv_paths.append(scalop)

dfs = [pd.read_csv(p, low_memory=False, on_bad_lines="skip") for p in csv_paths]
merged = pd.concat(dfs, ignore_index=True, sort=False)

# Put protein_id first if present
if "protein_id" in merged.columns:
    merged = merged[["protein_id", *[c for c in merged.columns if c != "protein_id"]]]

# --- Normalize helpers ---
def norm_str(s: pd.Series) -> pd.Series:
    return (
        s.astype("string")
         .str.strip()
         .replace({"": pd.NA, "nan": pd.NA, "None": pd.NA})
    )

# Signature column (id preferred over name)
signature_col = "signature_id" if "signature_id" in merged.columns else "signature_name"
merged[signature_col] = (
    merged[signature_col].astype("string")
    .str.replace(r"(?i)^interpro:\s*", "", regex=True)
    .str.strip()
    .replace({"": pd.NA})
)

# Ensure columns exist
for c in ["organism", "scientific_name", "tax_id"]:
    if c not in merged.columns:
        merged[c] = pd.NA

# --- Consolidate organism & scientific_name ---
org_norm = norm_str(merged["organism"])
sci_norm = norm_str(merged["scientific_name"])

# Optional: show differences (keep if you still want the check)
diff_mask = org_norm.notna() & sci_norm.notna() & (org_norm != sci_norm)
if diff_mask.any():
    print("Rows where 'organism' != 'scientific_name':")
    print(merged.loc[diff_mask, ["protein_id", "organism", "scientific_name"]]
                .drop_duplicates()
                .to_string(index=False))

# Prefer organism, else scientific_name
merged["organism"] = merged["organism"].where(org_norm.notna(), merged["scientific_name"])
merged = merged.drop(columns=["scientific_name"], errors="ignore")

# --- Fill missing tax_id by organism (vectorized) ---
tax_map = (
    merged.loc[merged["tax_id"].notna(), ["organism", "tax_id"]]
          .drop_duplicates()
          .set_index("organism")["tax_id"]
)
merged["tax_id"] = merged["tax_id"].fillna(merged["organism"].map(tax_map))

# --- Ordered-unique signatures per protein, spread to 4 cols ---
id_lists = (
    merged.groupby("protein_id", sort=False)[signature_col]
          .apply(lambda s: [x for x in pd.unique(s.dropna())])    # preserves order
          .reset_index(name="id_list")
)

def to_four(xs):
    xs4 = (xs[:4] + [pd.NA]*4)[:4]
    return pd.Series(xs4, index=["interpro_1", "interpro_2", "interpro_3", "interpro_4"])

wide = id_lists.join(id_lists["id_list"].apply(to_four)).drop(columns="id_list")

# --- Metadata: prefer first non-null per column (skipna) ---
# Replace empty strings with NA so GroupBy.first() can skip them
merged = merged.replace("", pd.NA)
meta_cols = [c for c in merged.columns if c not in [signature_col, "protein_id"]]
meta = merged.groupby("protein_id", as_index=False)[meta_cols].first()

# --- Final deduped dataframe ---
final_df = meta.merge(wide, on="protein_id", how="left")

# Drop domain columns if present
final_df = final_df.drop(columns=[c for c in ["domain_start","domain_end","domain_score","interpro_domain_name"] if c in final_df.columns])

# Quick stats
print(f"Rows: {len(merged):,} | Columns: {len(merged.columns)}")
if "protein_id" in merged.columns:
    print(f"Unique protein_id: {merged['protein_id'].nunique():,}")

final_df.to_csv(Path(out,"merged.csv"), index=False)




# --- Download UniProt sequences into FASTA (chunked, ≤100 ORs) ---
import time
import requests

protein_ids = pd.Series(final_df["protein_id"]).dropna().astype(str).unique().tolist()
print(f"Fetching {len(protein_ids):,} sequences from UniProt...")

def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i+n]

fasta_path = Path(out, "sequences.fasta")
missed = []

with open(fasta_path, "w") as fasta_file:
    for i, block in enumerate(chunks(protein_ids, 100), start=1):
        # Build query: accession:ID1 OR accession:ID2 ...
        query = " OR ".join(f"accession:{acc}" for acc in block)
        url = "https://rest.uniprot.org/uniprotkb/stream"
        params = {"format": "fasta", "query": query}

        # simple retry for transient 429/5xx
        for attempt in range(4):
            r = requests.get(url, params=params, timeout=120)
            if r.status_code == 200 and r.text.strip():
                fasta_file.write(r.text)
                break
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2 ** attempt)   # backoff
                continue
            # hard failure (e.g., 400) — record and move on
            missed.extend(block)
            print(f"[Block {i}] Error {r.status_code}: {r.text[:200]}")
            break

print(f"Saved FASTA to {fasta_path}")
if missed:
    print(f"Missed IDs (first 10 of {len(missed)}): {missed[:10]}")

from pathlib import Path
import pandas as pd

fasta_path = Path(out, "sequences.fasta")

# --- read all protein_ids from your DF
requested = pd.Series(final_df["protein_id"]).dropna().astype(str).unique().tolist()
requested_set = set(requested)

# --- parse accessions from fasta headers
def iter_uniprot_accessions_from_fasta(path):
    with open(path, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].strip()
                acc = None
                if "|" in header:               # e.g. sp|P12345|...
                    parts = header.split("|")
                    if len(parts) >= 2:
                        acc = parts[1].strip()
                if not acc:                     # fallback: first token
                    acc = header.split()[0].strip()
                acc = acc.split("-")[0]         # normalize isoforms
                yield acc

found_set = set(iter_uniprot_accessions_from_fasta(fasta_path))

missing = sorted(requested_set - found_set)
extras  = sorted(found_set - requested_set)

print(f"Requested: {len(requested_set):,}")
print(f"Found in FASTA: {len(found_set):,}")
print(f"Missing: {len(missing):,}")
print(f"Extras: {len(extras):,}")

if missing:
    print("\n--- Missing IDs (first 50) ---")
    print("\n".join(missing[:50]))

if extras:
    print("\n--- Extra IDs in FASTA (first 50) ---")
    print("\n".join(extras[:50]))
