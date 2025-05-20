import requests
from Bio import SeqIO
import os
from tqdm import tqdm

def check_alphafold(uniprot_id):
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    r = requests.head(url)
    return r.status_code == 200

script_dir = os.path.dirname(os.path.abspath(__file__))
fasta_file = os.path.join(script_dir, "../../data/interm/merged/merged_sanitized.fasta")
fasta_file = os.path.normpath(fasta_file)

records = list(SeqIO.parse(fasta_file, "fasta"))[100:200] # 100 for debugging
missing = []

for record in tqdm(records, desc="Checking AlphaFold entries"):
    uniprot_id = record.id
    if not check_alphafold(uniprot_id):
        missing.append(uniprot_id)

if missing:
    print("\nNo AlphaFold structure for UniProt IDs:")
    print("\n".join(missing))
else:
    print("\nAll entries have AlphaFold structures.")
