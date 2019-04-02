import os
import sys
import pickle
import base64
import PIL.Image

import numpy as np

from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from mstmap import mstmap
from operator import itemgetter
from faerun import Faerun
from matplotlib import pyplot as plt
from io import BytesIO
from fashion_mnist.utils import mnist_reader

dims = 2048
enc = mstmap.Minhash(dims)
# enc = mstmap.Minhash(784, sample_size=dims)
lf = mstmap.LSHForest(dims, 128, store=True)

images, labels = mnist_reader.load_mnist('fashion_mnist/data/fashion', kind='train')
images_2, labels_2 = mnist_reader.load_mnist('fashion_mnist/data/fashion', kind='t10k')

images = np.concatenate((images, images_2))
labels = np.concatenate((labels, labels_2))

# images = images[:2000]
# labels = labels[:2000]

image_count = len(images)

bar = Bar('Generating fingerprint', max=image_count, suffix='%(remaining)d numbers remaining')

tmp = []
for i, image in enumerate(images):
    avg = sum(image) / sum([1 if x > 0 else 0 for x in image])
    tmp.append(mstmap.VectorUchar([1 if x >= avg / 2 else 0 for x in image]))
    # tmp.append(mstmap.VectorFloat([x / 255.0 for x in image]))
    bar.next()

# mhs = enc.batch_from_weight_array(tmp)
lf.batch_add(enc.batch_from_binary_array(tmp))

bar.finish()

lf.index()

lf.store('mnist.dat')

print("Getting knn graph")

config = mstmap.LayoutConfiguration()
config.k = 10
config.kc = 25
config.sl_scaling_x = 1.0
config.sl_scaling_y = 1.0
config.sl_repeats = 5
config.sl_extra_scaling_steps = 2
config.placer = mstmap.Placer.Barycenter
config.merger = mstmap.Merger.LocalBiconnected
config.merger_factor = 2.0
config.merger_adjustment = 0
config.sl_scaling_type = mstmap.ScalingType.RelativeToDrawing
config.node_size = 1
config.mmm_repeats = 10

x, y, s, t = mstmap.layout_from_lsh_forest(lf, config)


image_labels = []

for image in images:
    img = PIL.Image.fromarray(np.uint8(np.split(np.array(image), 28)))
    buffered = BytesIO()
    img.save(buffered, format="JPEG")
    img_str = base64.b64encode(buffered.getvalue())
    image_labels.append('data:image/bmp;base64,' + str(img_str).replace('b\'', '').replace('\'', ''))


legend_labels = [(0, 'T-shirt/top'), (1, 'Trouser'), (2, 'Pullover'), (3, 'Dress'), (4, 'Coat'), 
                 (5, 'Sandal'), (6, 'Shirt'), (7, 'Sneaker'), (8, 'Bag'), (9, 'Ankle boot')]

faerun = Faerun(view='front', coords=False, title='Fashion MNIST')
faerun.add_scatter('fashion_mnist', { 'x': x, 'y': y, 'c': labels, 'labels': image_labels }, colormap='tab10', 
                   point_scale=1.0, categorical=True, has_legend=True, legend_labels=legend_labels)
faerun.add_tree('fashion_mnist_tree', { 'from': s, 'to': t }, point_helper='fashion_mnist', color='#222222')
faerun.plot('fashion_mnist', template='url_image')