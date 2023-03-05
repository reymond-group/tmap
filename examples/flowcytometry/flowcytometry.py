""" 
An example of tmap visualizing data gathered in a flow cytometry 
experiment. The k-nearest neighbor graph is constructed
using the Annoy library.

Data Source:
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0057002
"""

import numpy as np
import flowkit as fk
import tmap as tm
from faerun import Faerun
from annoy import AnnoyIndex
from scipy.spatial.distance import cosine as cosine_distance


PATHS = [
    "FR-FCM-ZZCF/K562 Cells_No Target Probes_002.fcs",
    "FR-FCM-ZZCF/K562 Cells_BCR-A647_001.fcs",
]
SKIP = 1


def load_data(sample: fk.Sample):
    channels = range(SKIP, len(sample.channels))
    channel_data = []

    for channel in channels:
        channel_data.append(np.array(sample.get_channel_events(channel, source="raw")))

    return np.array(channel_data).T


def load_time(sample: fk.Sample):
    """Assuming the time is channel 0"""
    return np.array(sample.get_channel_events(0, source="raw"))


def main():
    """Main function"""
    data = []
    time = []

    for path in PATHS:
        sample = fk.Sample(path)
        data.append(load_data(sample))
        time.append(load_time(sample))

    sources = []
    for i, e in enumerate(data):
        sources.extend([i] * len(e))

    data = np.concatenate(data, axis=0)
    time = np.concatenate(time, axis=0)

    d = len(data[0])

    # Initialize a new Annoy object and index it using 10 trees
    annoy = AnnoyIndex(d, metric="angular")
    for i, v in enumerate(data):
        annoy.add_item(i, v)
    annoy.build(10)

    # Create the k-nearest neighbor graph (k = 10)
    edge_list = []
    for i in range(len(data)):
        for j in annoy.get_nns_by_item(i, 10):
            edge_list.append((i, j, cosine_distance(data[i], data[j])))

    # Compute the layout from the edge list
    x, y, s, t, _ = tm.layout_from_edge_list(len(data), edge_list)

    legend_labels = [(0, "No Target Probe Negative Control"), (1, "Stained Sample")]

    # Create the plot
    faerun = Faerun(
        view="front",
        coords=False,
        legend_title="RNA Flow Cytometry: evaluation of detection sensitivity in low abundant intracellular RNA ",
    )
    faerun.add_scatter(
        "CYTO",
        {"x": x, "y": y, "c": sources, "labels": sources},
        point_scale=1.0,
        max_point_size=10,
        shader="smoothCircle",
        colormap="Set1",
        has_legend=True,
        categorical=True,
        legend_labels=legend_labels,
        legend_title="Cell Types",
    )
    faerun.add_tree(
        "CYTO_tree", {"from": s, "to": t}, point_helper="CYTO", color="#222222"
    )

    faerun.plot("cyto")


if __name__ == "__main__":
    main()
