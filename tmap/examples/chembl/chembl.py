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

    for i in [1]:
        for j in [0]:
            cfg = tm.LayoutConfiguration()
            cfg.node_size = 1 / 70

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

            faerun = Faerun(view="front", coords=False, title="ChEMBL")
            faerun.add_scatter(
                "chembl",
                {"x": x, "y": y, "c": vals, "labels": smiles},
                colormap="tab10",
                point_scale=2.5,
                max_point_size=10,
                has_legend=True,
                categorical=True,
                legend_labels=legend_labels,
            )
            faerun.add_tree(
                "chembl_tree",
                {"from": s, "to": t},
                point_helper="chembl",
                color="#222222",
            )

            faerun.plot("chembl" + str(i) + str(j), template="smiles")


if __name__ == "__main__":
    main()
