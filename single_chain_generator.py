from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select

parser = PDBParser(QUIET=True)
io = PDBIO()

input_folder = Path("pdb_files")
output_folder = Path("pdb_files_pp")
output_folder.mkdir(exist_ok=True)

results = []


class KeepOneChain(Select):
    def accept_chain(self, chain):
        return chain.id == self.chain_id

    def __init__(self, chain_id):
        self.chain_id = chain_id


for pdb_file in input_folder.glob("*.pdb"):
    structure = parser.get_structure(pdb_file.stem, str(pdb_file))
    chains = list(structure[0].get_chains())

    if not chains:
        continue

    original_chains = len(chains)
    first_chain = chains[0]

    io.set_structure(structure)
    io.save(
        str(output_folder / pdb_file.name),
        KeepOneChain(first_chain.id)
    )

    results.append({
        "PDB_ID": pdb_file.stem,
        "original_chains": original_chains,
        "final_chains": 1
    })

pd.DataFrame(results).to_csv(
    output_folder / "chain_summary.csv",
    index=False
)

print(f"Processed {len(results)} PDB files.")
