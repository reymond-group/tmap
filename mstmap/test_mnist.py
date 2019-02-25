import os
import sys
import pickle

import numpy as np
import multiprocessing as mp
import networkx as nx

from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from mstmap import mstmap
from operator import itemgetter
from faerun import Faerun



out_name = 'mnist'
dims = 2048
enc = mstmap.Minhash(dims)
# enc = mstmap.Minhash(784, sample_size=dims)
lf = mstmap.LSHForest(dims, 128, store=True)

print('Loading MNIST data ...')
mn = MNIST('./python-mnist/data')
images, labels = mn.load_training()

# images = images[:20000]
# labels = labels[:20000]

image_count = len(images)

bar = Bar('Generating fingerprint', max=image_count, suffix='%(remaining)d numbers remaining')

tmp = []
for i, image in enumerate(images):
    avg = sum(image) / sum([1 if x > 0 else 0 for x in image])
    tmp.append(mstmap.VectorUchar([1 if x >= avg else 0 for x in image]))
    # tmp.append(mstmap.VectorFloat([x / 255.0 for x in image]))
    bar.next()

mhs = enc.batch_from_binary_array(tmp)
# mhs = enc.batch_from_weight_array(tmp)
lf.batch_add(mhs)

bar.finish()

lf.index()

lf.store('test.dat')

coords = mstmap.layout_from_lsh_forest(lf)

with open('mnist.labels', 'w+') as f:
    for label in labels: f.write(str(label) + '\n')


x = coords[0]
y = coords[1]
z = [0] * len(coords[0])

faerun = Faerun(view='front', shader='legacyCircle', coords=False, point_size=0.25, tree_color='#ff0000')

print("Plotting")
faerun.plot({ 'x': x, 'y': y, 'c': labels, 'smiles': [''] * len(x) }, colormap='tab10')#, tree=edges)