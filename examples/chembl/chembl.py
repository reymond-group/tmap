import pickle
from timeit import default_timer as timer

import numpy as np
import pandas as pd
import tmap as tm

from faerun import Faerun


def main():
    """ Main function """

    dims = 512
    lf = tm.LSHForest(dims, 128, store=True)

    # Due to the large data size (> 1GB) the following files are not provided directly
    smiles, target_class, activity, chembl_id = pickle.load(open("chembl.pickle", "rb"))

    labels = []
    for i, s in enumerate(smiles):
        labels.append(
            s
            + "__"
            + chembl_id[i]
            + "__"
            + f'<a target="_blank" href="https://www.ebi.ac.uk/chembl/compound_report_card/{chembl_id[i]}">{chembl_id[i]}</a>'
        )

    lf.restore("chembl.dat")

    target_class_map = dict(
        [(y, x + 1) for x, y in enumerate(sorted(set(target_class)))]
    )

    classes = [
        "enzyme",
        "kinase",
        "protease",
        "cytochrome p450",
        "ion channel",
        "transporter",
        "transcription factor",
        "membrane receptor",
        "epigenetic regulator",
    ]

    i = 0
    for key, value in target_class_map.items():
        if key not in classes:
            target_class_map[key] = 7
        else:
            target_class_map[key] = i
            i += 1
            if i == 7:
                i = 8

    cfg = tm.LayoutConfiguration()
    cfg.node_size = 1 / 70
    cfg.mmm_repeats = 2
    cfg.sl_repeats = 2

    start = timer()
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, cfg)
    end = timer()
    print(end - start)

    activity = np.array(activity)
    activity = np.maximum(0.0, activity)
    activity = np.minimum(100.0, activity)
    activity = 10.0 - activity

    legend_labels = [
        (0, "Cytochrome p450"),
        (1, "Other Enzyme"),
        (2, "Epigenetic Regulator"),
        (3, "Ion Channel"),
        (4, "Kinase"),
        (5, "Membrane Receptor"),
        (6, "Protease"),
        (8, "Transcription Factor"),
        (9, "Transporter"),
        (7, "Other"),
    ]

    vals = [int(target_class_map[x]) for x in target_class]

    faerun = Faerun(view="front", coords=False)
    faerun.add_scatter(
        "chembl",
        {"x": x, "y": y, "c": vals, "labels": labels},
        colormap="tab10",
        point_scale=1.0,
        max_point_size=10,
        has_legend=True,
        categorical=True,
        shader="smoothCircle",
        legend_labels=legend_labels,
        title_index=1,
    )
    faerun.add_tree(
        "chembl_tree", {"from": s, "to": t}, point_helper="chembl", color="#222222"
    )

    faerun.plot("chembl", template="smiles")


if __name__ == "__main__":
    main()
