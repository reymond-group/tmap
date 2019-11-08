"""
Visualizing RNA sequencing data using tmap.

Data Source:
https://gdc.cancer.gov/about-data/publications/pancanatlas
"""

import numpy as np
import pandas as pd
from faerun import Faerun
import tmap as tm


# Coniguration for the tmap layout
CFG_TMAP = tm.LayoutConfiguration()
CFG_TMAP.k = 50
CFG_TMAP.kc = 50
CFG_TMAP.node_size = 1 / 20


DATA = pd.read_csv("data.csv.xz", index_col=0, sep=",")
LABELS = pd.read_csv("labels.csv", index_col=0, sep=",")

LABELMAP = {"PRAD": 1, "LUAD": 2, "BRCA": 3, "KIRC": 4, "COAD": 5}
LABELS = np.array([int(LABELMAP[v]) for v in LABELS["Class"]], dtype=np.int)


def main():
    """ Main function """

    # Initialize and configure tmap
    dims = 256
    enc = tm.Minhash(len(DATA.columns), 42, dims)
    lf = tm.LSHForest(dims * 2, 32, weighted=True)

    fps = []
    for _, row in DATA.iterrows():
        fps.append(tm.VectorFloat(list(row)))

    lf.batch_add(enc.batch_from_weight_array(fps))
    lf.index()

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, CFG_TMAP)
    lf.clear()

    legend_labels = {(1, "PRAD"), (2, "LUAD"), (3, "BRCA"), (4, "KIRC"), (5, "COAD")}

    # Create the plot
    faerun = Faerun(view="front", coords=False, legend_title="")
    faerun.add_scatter(
        "RNASEQ",
        {"x": x, "y": y, "c": LABELS, "labels": LABELS},
        colormap="tab10",
        point_scale=5.0,
        max_point_size=10,
        shader="smoothCircle",
        has_legend=True,
        categorical=True,
        legend_labels=legend_labels,
        legend_title="Tumor Types",
    )
    faerun.add_tree(
        "RNASEQ_tree", {"from": s, "to": t}, point_helper="RNASEQ", color="#666666"
    )
    faerun.plot("rnaseq")


if __name__ == "__main__":
    main()
