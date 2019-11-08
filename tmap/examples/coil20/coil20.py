""" Create plot and faerun html for coil_20"""
import base64
import os
from io import BytesIO

import numpy as np
import tmap as tm
from PIL import Image

from faerun import Faerun


def main():
    """ Main function """

    # Initialize and configure tmap
    dims = 2048
    enc = tm.Minhash(16384, 42, dims)
    lf = tm.LSHForest(dims * 2, 128, weighted=True)

    images = []
    labels = []
    image_labels = []

    for file in os.listdir("coil_20"):
        labels.append(int(file.split("__")[0].replace("obj", "")) - 1)
        images.append(list(Image.open("coil_20/" + file).getdata()))

    for image in images:
        img = Image.fromarray(np.uint8(np.split(np.array(image), 128)))
        buffered = BytesIO()
        img.save(buffered, format="JPEG")
        img_str = base64.b64encode(buffered.getvalue())
        image_labels.append(
            "data:image/bmp;base64," + str(img_str).replace("b'", "").replace("'", "")
        )

    tmp = []
    for _, image in enumerate(images):
        avg = sum(image) / sum([1 if x > 0 else 0 for x in image])
        tmp.append([i / 255 for i in image])

    lf.batch_add(enc.batch_from_weight_array(tmp))
    lf.index()

    x, y, s, t, _ = tm.layout_from_lsh_forest(lf)

    faerun = Faerun(clear_color="#111111", view="front", coords=False)
    faerun.add_scatter(
        "COIL20",
        {"x": x, "y": y, "c": labels, "labels": image_labels},
        colormap="tab20",
        shader="smoothCircle",
        point_scale=2.5,
        max_point_size=10,
        has_legend=True,
        categorical=True,
    )
    faerun.add_tree(
        "COIL20_tree", {"from": s, "to": t}, point_helper="COIL20", color="#666666"
    )
    faerun.plot("coil", template="url_image")


if __name__ == "__main__":
    main()
