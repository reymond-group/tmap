"""
Visualization of MiniBooNE experimental measurements using tmap.

Data Source:
https://archive.ics.uci.edu/ml/datasets/MiniBooNE+particle+identification
"""

import pandas as pd
from faerun import Faerun
from scipy.spatial.distance import cosine as cosine_distance
from annoy import AnnoyIndex
import tmap as tm


DATA = pd.read_csv("MiniBooNE_PID.txt.xz", sep=",", header=None)[::-1]

LABELS = [0] * 93565
LABELS.extend([1] * 36499)


def main():
    """ Main function """

    # Building a k-nearest neighbor graph using annoy and cosine distance
    annoy = AnnoyIndex(len(DATA.columns), metric="angular")
    annoy_graph = []

    for i, v in enumerate(DATA.values):
        annoy.add_item(i, v)
    annoy.build(10)

    for i in range(len(DATA)):
        for j in annoy.get_nns_by_item(i, 10):
            annoy_graph.append((i, j, cosine_distance(DATA.values[i], DATA.values[j])))

    # Creating the tmap layout
    x, y, s, t, _ = tm.layout_from_edge_list(len(DATA), annoy_graph)

    faerun = Faerun(view="front", coords=False)
    faerun.add_scatter(
        "MINIBOONE",
        {"x": x, "y": y, "c": LABELS, "labels": LABELS},
        shader="smoothCircle",
        colormap="Set1",
        point_scale=2.0,
        max_point_size=20,
        has_legend=True,
        categorical=True,
        legend_labels={(0, "Noise"), (1, "Signal")},
    )
    faerun.add_tree(
        "MINIBOONE_tree",
        {"from": s, "to": t},
        point_helper="MINIBOONE",
        color="#666666",
    )
    faerun.plot("miniboone", template="default")


if __name__ == "__main__":
    main()
