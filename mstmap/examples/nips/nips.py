"""
An example of tmap visualizing data extracted from
NIPS articles.

Data Source:
https://www.kaggle.com/benhamner/nips-papers#papers.csv
"""

import re

from collections import Counter

import pandas as pd
import tmap as tm

from faerun import Faerun


def main():
    """ The main function """
    df = pd.read_csv("papers.tar.xz")
    df.drop(df.tail(1).index, inplace=True)
    df["title"] = df["title"].apply(lambda t: t.replace("'", '"'))
    enc = tm.Minhash()
    lf = tm.LSHForest()

    ctr = Counter()
    texts = []
    for _, row in df.iterrows():
        text = re.sub(r"[^a-zA-Z-]+", " ", row["paper_text"])
        text = [t.lower() for t in text.split(" ") if len(t) > 2]
        ctr.update(text)
        texts.append(text)

    # Remove the top n words
    n = 6000
    ctr = ctr.most_common()[: -(len(ctr) - n) - 1 : -1]

    # Make it fast using a lookup map
    all_words = {}
    for i, (key, _) in enumerate(ctr):
        all_words[key] = i

    # Create the fingerprints and also check whether the word
    # "deep" is found in the document
    fingerprints = []
    has_word = []
    for text in texts:
        if "deep" in text:
            has_word.append(1)
        else:
            has_word.append(0)

        fingerprint = []
        for t in text:
            if t in all_words:
                fingerprint.append(all_words[t])
        fingerprints.append(tm.VectorUint(fingerprint))

    # Index the article fingerprints
    lf.batch_add(enc.batch_from_sparse_binary_array(fingerprints))
    lf.index()

    # Create the tmap
    config = tm.LayoutConfiguration()
    config.k = 100
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=config)

    faerun = Faerun(
        view="front", coords=False, legend_title="", legend_number_format="{:.0f}"
    )

    # Add a scatter that is bigger than the one above, to add colored
    # circles.
    faerun.add_scatter(
        "NIPS_word",
        {"x": x, "y": y, "c": has_word, "labels": df["title"]},
        colormap="Set1",
        point_scale=7.5,
        max_point_size=25,
        shader="smoothCircle",
        has_legend=True,
        categorical=True,
        legend_title="Contains word<br/>'deep'",
        legend_labels=[(0, "No"), (1, "Yes")],
        interactive=False,
    )

    # Add a scatter that is colored by year on top
    faerun.add_scatter(
        "NIPS",
        {"x": x, "y": y, "c": df["year"], "labels": df["title"]},
        colormap="gray",
        point_scale=5.0,
        max_point_size=20,
        shader="smoothCircle",
        has_legend=True,
        legend_title="Year of<br/>Publication",
    )

    faerun.add_tree(
        "NIPS_tree", {"from": s, "to": t}, point_helper="NIPS", color="#666666"
    )

    faerun.plot("nips_papers")


if __name__ == "__main__":
    main()
