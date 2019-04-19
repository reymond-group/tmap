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



out_name = 'mnist'
dims = 2048
enc = mstmap.Minhash(dims)
# enc = mstmap.Minhash(784, sample_size=dims)
lf = mstmap.LSHForest(dims, 128, store=True)

print('Loading MNIST data ...')
mn = MNIST('./python-mnist/data')
images, labels = mn.load_training()
images_2, labels_2 = mn.load_testing()
images.extend(images_2)
labels.extend(labels_2)


# images = images[:2000]
# labels = labels[:2000]

image_count = len(images)

bar = Bar('Generating fingerprint', max=image_count, suffix='%(remaining)d numbers remaining')

tmp = []
for i, image in enumerate(images):
    avg = sum(image) / sum([1 if x > 0 else 0 for x in image])
    tmp.append(mstmap.VectorUchar([1 if x >= avg else 0 for x in image]))
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
config.sl_scaling_min = 1.0
config.sl_scaling_max = 1.0
config.sl_repeats = 5
config.sl_extra_scaling_steps = 2
config.placer = mstmap.Placer.Barycenter
config.merger = mstmap.Merger.LocalBiconnected
config.merger_factor = 2.0
config.merger_adjustment = 0
config.sl_scaling_type = mstmap.ScalingType.RelativeToDrawing
config.node_size = 1
config.mmm_repeats = 10

start = timer()
x, y, s, t, props = mstmap.layout_from_lsh_forest(lf, config)
print(timer() - start)
print(max(props.degrees))
print(sum(props.degrees) / len(props.degrees))


image_labels = []

for image in images:
    img = PIL.Image.fromarray(np.uint8(np.split(np.array(image), 28)))
    buffered = BytesIO()
    img.save(buffered, format="JPEG")
    img_str = base64.b64encode(buffered.getvalue())
    image_labels.append('data:image/bmp;base64,' + str(img_str).replace('b\'', '').replace('\'', ''))

faerun = Faerun(view='front', coords=False, title='MNIST')
faerun.add_scatter('mnist', { 'x': x, 'y': y, 'c': labels, 'labels': image_labels }, colormap='tab10', 
                   point_scale=1.0, categorical=True, has_legend=True)
faerun.add_tree('mnist_tree', { 'from': s, 'to': t }, point_helper='mnist', color='#222222')
faerun.plot('mnist', template='url_image')