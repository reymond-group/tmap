"""
Visualizing the Drugbank chemical space using tmap.
"""

import tmap as tm
import numpy as np
import scipy.stats as ss
from faerun import Faerun
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import AllChem


def main():
    """ The main function """
    mhfp = MHFPEncoder()
    lf = tm.LSHForest(2048, 128)

    fps = []
    smiles = []
    values = []
    with open("drugbank.smi", "r") as f:
        for i, line in enumerate(f):
            smi = line.split("\t")[0]
            mol = AllChem.MolFromSmiles(smi)
            if mol is not None:
                hac = mol.GetNumHeavyAtoms()
                if hac > 2:
                    fps.append(tm.VectorUint(mhfp.encode_mol(mol)))
                    smiles.append(smi)
                    values.append(hac)

    lf.batch_add(fps)
    lf.index()

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf)

    # Rank the molecules between 0.0 and 1.0 on values
    ranked_values = ss.rankdata(np.array(values) / max(values)) / len(values)

    # Create the plot
    faerun = Faerun(view="front", coords=False, title="Drugbank", legend_title="")
    faerun.add_scatter(
        "DRUGBANK",
        {"x": x, "y": y, "c": ranked_values, "labels": smiles},
        point_scale=2.0,
        max_point_size=10,
        shader="smoothCircle",
        has_legend=True,
        legend_title="Heavy Atom Count",
        max_legend_label=str(max(values)),
        min_legend_label=str(min(values)),
    )
    faerun.add_tree(
        "DRUGBANK_tree", {"from": s, "to": t}, point_helper="DRUGBANK", color="#666666"
    )
    faerun.plot("drugbank", template="smiles")


if __name__ == "__main__":
    main()
