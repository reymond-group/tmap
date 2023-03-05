import numpy as np
import pandas as pd
import tmap as tm
import scipy.stats as ss
from faerun import Faerun
from rdkit.Chem import AllChem, Descriptors, Descriptors3D
from mhfp.encoder import MHFPEncoder
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from lipinski import lipinski_pass


def main():
    """The main function"""
    df = pd.read_csv("drugbank.csv").dropna(subset=["SMILES"]).reset_index(drop=True)
    enc = MHFPEncoder()
    lf = tm.LSHForest(2048, 128)

    fps = []
    labels = []
    groups = []
    tpsa = []
    logp = []
    mw = []
    h_acceptors = []
    h_donors = []
    ring_count = []
    is_lipinski = []
    has_coc = []
    has_sa = []
    has_tz = []

    substruct_coc = AllChem.MolFromSmiles("COC")
    substruct_sa = AllChem.MolFromSmiles("NS(=O)=O")
    substruct_tz = AllChem.MolFromSmiles("N1N=NN=C1")

    total = len(df)
    for i, row in df.iterrows():
        if i % 1000 == 0 and i > 0:
            print(f"{round(100 * (i / total))}% done ...")

        smiles = row[6]
        mol = AllChem.MolFromSmiles(smiles)

        if mol and mol.GetNumAtoms() > 5 and smiles.count(".") < 2:
            fps.append(tm.VectorUint(enc.encode_mol(mol, min_radius=0)))
            labels.append(
                f'{smiles}__<a href="https://www.drugbank.ca/drugs/{row[0]}" target="_blank">{row[0]}</a>__{row[1]}'.replace(
                    "'", ""
                )
            )
            groups.append(row[3].split(";")[0])
            tpsa.append(Descriptors.TPSA(mol))
            logp.append(Descriptors.MolLogP(mol))
            mw.append(Descriptors.MolWt(mol))
            h_acceptors.append(Descriptors.NumHAcceptors(mol))
            h_donors.append(Descriptors.NumHDonors(mol))
            ring_count.append(Descriptors.RingCount(mol))
            is_lipinski.append(lipinski_pass(mol))
            has_coc.append(mol.HasSubstructMatch(substruct_coc))
            has_sa.append(mol.HasSubstructMatch(substruct_sa))
            has_tz.append(mol.HasSubstructMatch(substruct_tz))

    # Create the labels and the integer encoded array for the groups,
    # as they're categorical
    labels_groups, groups = Faerun.create_categories(groups)
    tpsa_ranked = ss.rankdata(np.array(tpsa) / max(tpsa)) / len(tpsa)
    logp_ranked = ss.rankdata(np.array(logp) / max(logp)) / len(logp)
    mw_ranked = ss.rankdata(np.array(mw) / max(mw)) / len(mw)
    h_acceptors_ranked = ss.rankdata(np.array(h_acceptors) / max(h_acceptors)) / len(
        h_acceptors
    )
    h_donors_ranked = ss.rankdata(np.array(h_donors) / max(h_donors)) / len(h_donors)
    ring_count_ranked = ss.rankdata(np.array(ring_count) / max(ring_count)) / len(
        ring_count
    )

    lf.batch_add(fps)
    lf.index()
    cfg = tm.LayoutConfiguration()
    cfg.k = 100
    # cfg.sl_extra_scaling_steps = 1
    cfg.sl_repeats = 2
    cfg.mmm_repeats = 2
    cfg.node_size = 2
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)

    # Define a colormap highlighting approved vs non-approved
    custom_cmap = ListedColormap(
        ["#2ecc71", "#9b59b6", "#ecf0f1", "#e74c3c", "#e67e22", "#f1c40f", "#95a5a6"],
        name="custom",
    )

    bin_cmap = ListedColormap(["#e74c3c", "#2ecc71"], name="bin_cmap")

    f = Faerun(
        clear_color="#222222",
        coords=False,
        view="front",
        impress='made with <a href="http://tmap.gdb.tools" target="_blank">tmap</a><br />and <a href="https://github.com/reymond-group/faerun-python" target="_blank">faerun</a><br /><a href="https://gist.github.com/daenuprobst/5cddd0159c0cf4758fb16b4b4acbef89">source</a>',
    )

    f.add_scatter(
        "Drugbank",
        {
            "x": x,
            "y": y,
            "c": [
                groups,
                is_lipinski,
                has_coc,
                has_sa,
                has_tz,
                tpsa_ranked,
                logp_ranked,
                mw_ranked,
                h_acceptors_ranked,
                h_donors_ranked,
                ring_count_ranked,
            ],
            "labels": labels,
        },
        shader="smoothCircle",
        colormap=[
            custom_cmap,
            bin_cmap,
            bin_cmap,
            bin_cmap,
            bin_cmap,
            "viridis",
            "viridis",
            "viridis",
            "viridis",
            "viridis",
            "viridis",
        ],
        point_scale=2.5,
        categorical=[True, True, True, True, True, False, False, False, False, False],
        has_legend=True,
        legend_labels=[
            labels_groups,
            [(0, "No"), (1, "Yes")],
            [(0, "No"), (1, "Yes")],
            [(0, "No"), (1, "Yes")],
            [(0, "No"), (1, "Yes")],
        ],
        selected_labels=["SMILES", "Drugbank ID", "Name"],
        series_title=[
            "Group",
            "Lipinski",
            "Ethers",
            "Sulfonamides",
            "Tetrazoles",
            "TPSA",
            "logP",
            "Mol Weight",
            "H Acceptors",
            "H Donors",
            "Ring Count",
        ],
        max_legend_label=[
            None,
            None,
            None,
            None,
            None,
            str(round(max(tpsa))),
            str(round(max(logp))),
            str(round(max(mw))),
            str(round(max(h_acceptors))),
            str(round(max(h_donors))),
            str(round(max(ring_count))),
        ],
        min_legend_label=[
            None,
            None,
            None,
            None,
            None,
            str(round(min(tpsa))),
            str(round(min(logp))),
            str(round(min(mw))),
            str(round(min(h_acceptors))),
            str(round(min(h_donors))),
            str(round(min(ring_count))),
        ],
        title_index=2,
        legend_title="",
    )

    f.add_tree("drugbanktree", {"from": s, "to": t}, point_helper="Drugbank")

    f.plot("drugbank", template="smiles")


if __name__ == "__main__":
    main()
