"""
Visualizing the Drugbank chemical space using tmap.
"""

import tmap as tm
import numpy as np
import scipy.stats as ss
from faerun import Faerun
from mhfp.encoder import MHFPEncoder
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors


def main():
    """ The main function """
    mhfp = MHFPEncoder()
    lf = tm.LSHForest(2048, 128)

    fps = []
    smiles = []
    values = []
    values_b = []
    values_c = []
    values_d = []
    substruct = AllChem.MolFromSmiles("COC")
    with open("drugbank.smi", "r") as f:
        for i, line in enumerate(f):
            smi = line.split("\t")[0]
            mol = AllChem.MolFromSmiles(smi)
            if mol is not None:
                hac = mol.GetNumHeavyAtoms()
                n_rings = len(AllChem.GetSymmSSSR(mol))
                logp = Descriptors.MolLogP(mol)
                has_subs = mol.HasSubstructMatch(substruct)
                if hac > 2:
                    fps.append(tm.VectorUint(mhfp.encode_mol(mol)))
                    smiles.append(smi)
                    values.append(hac)
                    values_b.append(n_rings)
                    values_c.append(logp)
                    values_d.append((1 if has_subs else 0))

    lf.batch_add(fps)
    lf.index()

    cfg = tm.LayoutConfiguration()
    cfg.k = 100
    cfg.node_size = 1 / 20

    f = tm.VectorUint()
    t = tm.VectorUint()
    w = tm.VectorFloat()

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)

    # Rank the molecules between 0.0 and 1.0 on values
    ranked_values = ss.rankdata(np.array(values) / max(values)) / len(values)
    ranked_values_c = ss.rankdata(np.array(values_c) / max(values_c)) / len(values_c)

    # Create the plot
    faerun = Faerun(view="front", coords=False, title="Drugbank", legend_title="")
    faerun.add_scatter(
        "DRUGBANK",
        {
            "x": x,
            "y": y,
            "c": [ranked_values, values_b, ranked_values_c, values_d],
            "labels": smiles,
        },
        colormap=["viridis", "tab20", "plasma", "Dark2"],
        point_scale=2.0,
        max_point_size=10,
        shader="smoothCircle",
        has_legend=True,
        categorical=[False, True, False, True],
        series_title=[
            "Heavy Atom Count",
            "Number of Rings",
            "Computed logP",
            "Has COC",
        ],
        legend_title=["Drugbank"],
        max_legend_label=[str(max(values)), None, str(round(max(values_c)))],
        min_legend_label=[str(min(values)), None, str(round(min(values_c)))],
    )
    faerun.add_tree(
        "DRUGBANK_tree", {"from": s, "to": t}, point_helper="DRUGBANK", color="#666666"
    )
    faerun.plot("drugbank", template="smiles")


if __name__ == "__main__":
    main()
