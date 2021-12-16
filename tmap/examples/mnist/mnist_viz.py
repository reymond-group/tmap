"""
Visualization of the MNIST image data using tmap.
"""

import base64
from io import BytesIO

import numpy as np
from faerun import Faerun
from mnist import MNIST
from PIL import Image

import tmap as tm


# Load the data
MN = MNIST("./python-mnist/bin/data")
IMAGES_TRAIN, LABELS_TRAIN = MN.load_training()
IMAGES_TEST, LABELS_TEST = MN.load_testing()

IMAGES = np.concatenate((IMAGES_TRAIN, IMAGES_TEST))
LABELS = np.concatenate((LABELS_TRAIN, LABELS_TEST))
IMAGE_LABELS = []

# Coniguration for the tmap layout
CFG = tm.LayoutConfiguration()
CFG.node_size = 1 / 50


def main():
    """ Main function """

    # Initialize and configure tmap
    dims = 1024
    enc = tm.Minhash(dims)
    lf = tm.LSHForest(dims, 128)

    print("Converting images ...")
    for image in IMAGES:
        img = Image.fromarray(np.uint8(np.split(np.array(image), 28)))
        buffered = BytesIO()
        img.save(buffered, format="JPEG")
        img_str = base64.b64encode(buffered.getvalue())
        IMAGE_LABELS.append(
            "data:image/bmp;base64," + str(img_str).replace("b'", "").replace("'", "")
        )
    tmp = []
    for _, image in enumerate(IMAGES):
        avg = sum(image) / sum([1 if x > 0 else 0 for x in image])
        tmp.append(tm.VectorUchar([1 if x >= avg else 0 for x in image]))
        # tmp.append(tm.VectorUint(image))

    print("Running tmap ...")
    lf.batch_add(enc.batch_from_binary_array(tmp))
    # LF.batch_add(ENC.batch_from_int_weight_array(tmp))
    lf.index()

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, CFG)

    faerun = Faerun(clear_color="#111111", view="front", coords=False)
    faerun.add_scatter(
        "MNIST",
        {"x": x, "y": y, "c": LABELS, "labels": IMAGE_LABELS},
        colormap="tab10",
        shader="smoothCircle",
        point_scale=2.5,
        max_point_size=10,
        has_legend=True,
        categorical=True,
    )
    faerun.add_tree(
        "MNIST_tree", {"from": s, "to": t}, point_helper="MNIST", color="#666666"
    )
    faerun.plot("mnist", template="url_image")


if __name__ == "__main__":
    main()
