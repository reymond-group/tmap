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
    "../flowcytometry/FR-FCM-ZZCF/K562 Cells_No Target Probes_002.fcs",
    "../flowcytometry/FR-FCM-ZZCF/K562 Cells_BCR-A647_001.fcs",
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

    print("Loading data ...")

    for path in PATHS:
        sample = fk.Sample(path)
        data.append(load_data(sample))
        time.append(load_time(sample))

    sources = []
    for i, e in enumerate(data):
        sources.extend([i] * len(e))

    data = np.concatenate(data, axis=0)
    rand_idx = np.random.randint(data.shape[0], size=5000)
    print(rand_idx)
    data = data[rand_idx, :]
    sources = np.array(sources)[rand_idx]

    legend_labels = [(0, "No Target Probe Negative Control"), (1, "Stained Sample")]

    print("Data shape:", data.shape)
    print("Embedding ...")

    te = tm.embed(
        data,
        layout_generator=tm.layout_generators.AnnoyLayoutGenerator(),
        # layout_generator=tm.layout_generators.BuiltinLayoutGenerator(),
        keep_knn=True,
    )

    print("Plotting ...")
    tm.plot(
        te,
        show=True,
        line_kws={"linestyle": "--", "color": "gray"},
        scatter_kws={"s": 5, "c": sources},
    )


if __name__ == "__main__":
    main()
