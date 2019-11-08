"""
Create plot and faerun html for the gutenberg data set.
The data has been pre-processed in similar fashion as the NIPS dataset.

Data Source:
https://web.eecs.umich.edu/~lahiri/gutenberg_dataset.html
"""

import gzip

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from faerun import Faerun
import tmap as tm

# Coniguration for the tmap layout
CFG_TMAP = tm.LayoutConfiguration()
CFG_TMAP.k = 50
CFG_TMAP.kc = 50
CFG_TMAP.fme_iterations = 1000
CFG_TMAP.node_size = 1 / 10

# Load and preparte the data
DATA = []

with gzip.open("books_indexed_all.csv.gz", "rt") as f:
    for line in f:
        DATA.append(list(map(int, line.split(","))))

LABELS = []
FAERUN_LABELS = []

with open("books_labels_all.csv", "r") as f:
    for line in f:
        line = line.strip()
        author_title = line.split(",")
        LABELS.append(author_title[0])
        fl = author_title[1] + " (" + author_title[0] + ")"
        FAERUN_LABELS.append(
            fl.replace("'", "")
            + "__"
            + author_title[0].replace("'", "")
            + "__"
            + author_title[1].replace("'", "")
        )

labels_unique = set(LABELS)
labels_map = {}
for i, val in enumerate(labels_unique):
    labels_map[val] = i

# for i, val in enumerate(LABELS):
#     LABELS[i] = labels_map[val]

AUTHORS = [
    "Rudyard Kipling",
    "Herbert George Wells",
    "Charles Darwin",
    "George Bernard Shaw",
    "William Wymark Jacobs",
]

lbl_tmp = []
for label in LABELS:
    if "Rudyard Kipling" in label:
        lbl_tmp.append(1)
    elif "Herbert George Wells" in label:
        lbl_tmp.append(2)
    elif "Charles Darwin" in label:
        lbl_tmp.append(3)
    elif "George Bernard Shaw" in label:
        lbl_tmp.append(4)
    elif "William Wymark Jacobs" in label:
        lbl_tmp.append(5)
    else:
        lbl_tmp.append(0)

LABELS = lbl_tmp


def main():
    """ Main function """

    # Initialize and configure tmap
    dims = 2048
    enc = tm.Minhash(dims)
    lf = tm.LSHForest(dims, 128, store=True)

    fps = []
    # fps_umap = []
    for row in DATA:
        fps.append(tm.VectorUint(list(row)))

    lf.batch_add(enc.batch_from_sparse_binary_array(fps))
    lf.index()

    x_tmap, y_tmap, s, t, _ = tm.layout_from_lsh_forest(lf, CFG_TMAP)
    lf.clear()

    # Prepare custom color map
    tab10 = plt.get_cmap("tab10").colors
    colors_gray = [(0.2, 0.2, 0.2), tab10[0], tab10[1], tab10[2], tab10[3], tab10[4]]
    custom_cm_gray = LinearSegmentedColormap.from_list(
        "custom_cm_gray", colors_gray, N=len(colors_gray)
    )

    legend_labels = [
        (1, "Rudyard Kipling"),
        (2, "Herbert George Wells"),
        (3, "Charles Darwin"),
        (4, "George Bernard Shaw"),
        (5, "William Wymark Jacobs"),
        (0, "Other"),
    ]

    faerun = Faerun(
        clear_color="#111111",
        view="front",
        coords=False,
        alpha_blending=True,
        legend_title="",
    )
    faerun.add_scatter(
        "gutenberg",
        {"x": x_tmap, "y": y_tmap, "c": LABELS, "labels": FAERUN_LABELS},
        colormap=custom_cm_gray,
        point_scale=4.2,
        max_point_size=10,
        has_legend=True,
        categorical=True,
        legend_title="Authors",
        legend_labels=legend_labels,
        shader="smoothCircle",
        selected_labels=["Author", "Title"],
    )
    faerun.add_tree(
        "gutenberg_tree",
        {"from": s, "to": t},
        point_helper="gutenberg",
        color="#222222",
    )
    faerun.plot("gutenberg", template="default")


if __name__ == "__main__":
    main()
