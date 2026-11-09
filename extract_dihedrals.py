from Bio.PDB import PDBParser, PPBuilder
from Bio.SeqUtils import seq1
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

parser = PDBParser(QUIET=True)
builder = PPBuilder()
MIN_WINDOW = 3
MAX_WINDOW = 17
BACKBONE = ("N", "CA", "C", "O")

# PDB IDs flagged by the previous script
suspicious = set(pd.read_csv("suspicious_pdbs.csv")["PDB_ID"].astype(str))

# Keep only the files that are not suspicious
all_files = list(Path("pdb_files_pp").glob("*.pdb"))
good_files = [
    pdb_file for pdb_file in all_files
    if pdb_file.stem not in suspicious
]

print(f"{len(good_files)} good structures out of {len(all_files)}.")

fragments = []
phis = []
psis = []

for pdb_file in good_files:
    structure = parser.get_structure(pdb_file.stem, str(pdb_file))

    for model in structure:
        for chain in model:
            for peptide in builder.build_peptides(chain):
                # Angles for the whole peptide, so every inner residue has both phi and psi
                angles = peptide.get_phi_psi_list()

                # All angles, for the histogram
                phis += [np.degrees(angle[0]) for angle in angles if angle[0] is not None]
                psis += [np.degrees(angle[1]) for angle in angles if angle[1] is not None]

                for window in range(MIN_WINDOW, MAX_WINDOW + 1):
                    for start in range(len(peptide) - window + 1):
                        fragment = peptide[start:start + window]
                        fragment_angles = angles[start:start + window]

                        # Skip fragments at a chain end (missing phi or psi)
                        if any(phi is None or psi is None for phi, psi in fragment_angles):
                            continue

                        # Skip fragments with incomplete backbone atoms
                        if not all(res.has_id(atom) for res in fragment for atom in BACKBONE):
                            continue

                        fragments.append({
                            "PDB_ID": pdb_file.stem,
                            "model": model.id,
                            "chain": chain.id,
                            "window": window,
                            "start_residue": fragment[0].id[1],
                            "end_residue": fragment[-1].id[1],
                            "sequence": "".join(seq1(res.get_resname()) for res in fragment),
                            "phi": np.degrees([phi for phi, psi in fragment_angles]),
                            "psi": np.degrees([psi for phi, psi in fragment_angles]),
                            "backbone_coords": np.array([res[atom].coord for res in fragment for atom in BACKBONE])
                        })

df = pd.DataFrame(fragments)

# Full data: phi and psi (in degrees) and backbone coordinates of each fragment
df.to_pickle("fragments_dihedrals.pkl")

# Readable summary without the arrays
df.drop(columns=["phi", "psi", "backbone_coords"]).to_csv("fragments_dihedrals.csv", index=False)

print(f"Extracted {len(fragments)} fragments.")

# Histogram of all phi and psi angles
plt.hist(phis, bins=72, label=r"$\phi$", alpha=1)
plt.hist(psis, bins=72, label=r"$\psi$", alpha=0.5)
plt.xlabel("Angle (degrees)")
plt.ylabel("Counts")
plt.legend()
plt.savefig("phi_psi_histogram.png", dpi=150)