import pandas as pd
import requests, time
from pathlib import Path
import taxopy

# Paths
merged_path = Path("../../data/interm/merged/merged.csv")
out_path = Path("../../data/interm/merged/merged_with_lineage.csv")

# Load CSV
df = pd.read_csv(merged_path, low_memory=False)
if "tax_id" not in df.columns:
    df["tax_id"] = pd.NA

# Fetch missing tax_id
to_fetch = df[df["tax_id"].isna()]["protein_id"].dropna().unique().tolist()
print(f"Fetching tax_id for {len(to_fetch):,} proteins…")

def fetch_taxids(protein_ids, batch=100):
    tax_map = {}
    for i in range(0, len(protein_ids), batch):
        block = protein_ids[i : i + batch]
        query = " OR ".join(f"accession:{acc}" for acc in block)
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": query,
            "fields": "accession,organism_id",
            "format": "tsv",
            "size": batch
        }
        for attempt in range(4):
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200:
                for line in r.text.strip().splitlines()[1:]:
                    acc, taxid = line.split("\t")
                    tax_map[acc.split("-")[0]] = taxid
                break
            if r.status_code in (429, 500, 502, 503, 504):
                time.sleep(2**attempt)
                continue
            print(f"Error {r.status_code}: {r.text[:200]}")
            break
    return tax_map

if to_fetch:
    tax_map = fetch_taxids(to_fetch)
    df["tax_id"] = df.apply(
        lambda r: tax_map.get(str(r["protein_id"]).split("-")[0], r["tax_id"]), axis=1
    )

# Initialize taxopy
taxdump_url = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Initialize TaxDb with auto-download enabled
taxdb = taxopy.TaxDb(taxdump_url=taxdump_url, keep_files=True)

def get_lineage(tid):
    try:
        tax = taxopy.Taxon(int(tid), taxdb)
        lineage = tax.rank_name_dictionary
        return pd.Series({
            "kingdom": lineage.get("kingdom"),
            "phylum": lineage.get("phylum"),
            "class": lineage.get("class"),
            "order": lineage.get("order"),
            "family": lineage.get("family"),
            "genus": lineage.get("genus")
        })
    except:
        return pd.Series(dict.fromkeys(["kingdom","phylum","class","order","family","genus"], None))

df = df.join(df["tax_id"].dropna().apply(get_lineage))
df.to_csv(out_path, index=False)
print(f"Saved → {out_path}")
