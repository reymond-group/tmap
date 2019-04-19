import os
import sys
import pickle
import base64
import numpy as np

from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from tmap import tmap
from operator import itemgetter
from faerun import Faerun
from matplotlib import pyplot as plt
from io import BytesIO
from PIL import Image

dims = 2048
enc = tmap.Minhash(dims)
# enc = tmap.Minhash(784, sample_size=dims)
lf = tmap.LSHForest(dims, 128, store=True)


images = []
labels = []
for file in os.listdir('coil_20'):
    label = int(file.split('__')[0].replace('obj', ''))
    labels.append(label - 1)
    img = Image.open('coil_20/' + file)
    images.append(list(img.getdata()))

# images = images[:2000]
# labels = labels[:2000]

image_count = len(images)

bar = Bar('Generating fingerprint', max=image_count, suffix='%(remaining)d numbers remaining')

tmp = []
for i, image in enumerate(images):
    avg = sum(image) / sum([1 if x > 0 else 0 for x in image])
    tmp.append(tmap.VectorUchar([1 if x >= 122 else 0 for x in image]))
    # tmp.append(tmap.VectorFloat([x / 255.0 for x in image]))
    bar.next()

# mhs = enc.batch_from_weight_array(tmp)
lf.batch_add(enc.batch_from_binary_array(tmp))

bar.finish()

print(len(images))
lf.index()

lf.store('mnist.dat')

print("Getting knn graph")

config = tmap.LayoutConfiguration()
config.k = 10
config.kc = 50
config.sl_scaling_x = 1.0
config.sl_scaling_y = 1.0
config.sl_repeats = 5
config.sl_extra_scaling_steps = 2
config.placer = tmap.Placer.Barycenter
config.merger = tmap.Merger.LocalBiconnected
config.merger_factor = 2.0
config.merger_adjustment = 0
config.sl_scaling_type = tmap.ScalingType.RelativeToDrawing
config.node_size = 50
config.mmm_repeats = 50

x, y, s, t = tmap.layout_from_lsh_forest(lf, config)


image_labels = []

for image in images:
    img = Image.fromarray(np.uint8(np.split(np.array(image), 128)))
    buffered = BytesIO()
    img.save(buffered, format="JPEG")
    img_str = base64.b64encode(buffered.getvalue())
    image_labels.append('data:image/bmp;base64,' + str(img_str).replace('b\'', '').replace('\'', ''))


# legend_labels = [(0, 'T-shirt/top'), (1, 'Trouser'), (2, 'Pullover'), (3, 'Dress'), (4, 'Coat'), 
#                  (5, 'Sandal'), (6, 'Shirt'), (7, 'Sneaker'), (8, 'Bag'), (9, 'Ankle boot')]

faerun = Faerun(view='front', coords=False, title='Fashion MNIST')
faerun.add_scatter('coil20', { 'x': x, 'y': y, 'c': labels, 'labels': image_labels }, colormap='tab20', 
                   point_scale=5.0, categorical=True, has_legend=True) #, legend_labels=legend_labels)
faerun.add_tree('coil20_tree', { 'from': s, 'to': t }, point_helper='coil20', color='#222222')
faerun.plot('coil20', template='url_image')