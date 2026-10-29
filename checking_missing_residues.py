from Bio.PDB import PDBParser
from pathlib import Path
import numpy as np
import pandas as pd

parser = PDBParser(QUIET=True)
THRESHOLD = 4.5

results = []

for pdb_file in Path("pdb_files_pp").glob("*.pdb"):
    structure = parser.get_structure(pdb_file.stem, str(pdb_file))

    for model in structure:
        for chain in model:
            residues = [
                res for res in chain
                if res.id[0] == " " and res.has_id("CA")
            ]

            for prev, curr in zip(residues, residues[1:]):
                distance = np.linalg.norm(
                    curr["CA"].coord - prev["CA"].coord
                )

                if distance > THRESHOLD:
                    results.append({
                        "PDB_ID": pdb_file.stem,
                        "chain": chain.id,
                        "residue_1": prev.id[1],
                        "residue_2": curr.id[1],
                        "CA_distance": round(float(distance), 3)
                    })

pd.DataFrame(results).to_csv("suspicious_pdbs.csv", index=False)

print(f"Found {len(results)} suspicious connections.")
