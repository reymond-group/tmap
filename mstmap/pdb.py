import os
import sys
import pickle

import numpy as np
import pandas as pd
import scipy.stats as ss

from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from mstmap import mstmap
from operator import itemgetter
from faerun import Faerun


def convert_the_panda(value):
    return mstmap.VectorUint(list(map(int, value.split(','))))

def convert_the_panda_float(value):
    return list(map(float, value.split(',')))

out_name = 'pdb'
f = 512
enc = mstmap.Minhash(136, 42, f)
lf = mstmap.LSHForest(f, 512, store=True, file_backed=False)
weighted = True

print('Loading pdb data ...')

labels = []
values = [[], [], [], [], [], [], []]
index = 0
chunk_id = 0
for chunk in pd.read_csv('/media/daenu/Even More Data/pdb/full_fp.csv', sep=',', header=None, chunksize=20000):
    print(chunk_id)
    if chunk_id > 9: break
    chunk_id += 1
    fps = []

    for _, record in chunk.iterrows():
        labels.append(record[0])
        for i in range(7):
            values[i].append(record[137 + i])
        
        fps.append(mstmap.VectorFloat([float(i) for i in record[1:137]]))
        index += 1

    start = timer()
    lf.batch_add(enc.batch_from_weight_array(fps))
    end = timer()
    print(end - start)

start = timer()
lf.index()
end = timer()

lf.store('pdb.dat')
with open('pdb.pickle', 'wb+') as handle:
    pickle.dump((labels, values), handle, protocol=pickle.HIGHEST_PROTOCOL)


# lf.restore('pdb.dat')
# labels, values = pickle.load(open('pdb.pickle', 'rb'))

print("Getting knn graph")

config = mstmap.LayoutConfiguration()
config.k = 10
config.kc = 1000
config.sl_scaling_x = 5
config.sl_scaling_y = 10
config.placer = mstmap.Placer.Solar
config.merger = mstmap.Merger.EdgeCover
config.merger_factor = 2.0
config.sl_scaling_type = mstmap.ScalingType.RelativeToAvgLength
# config.fme_iterations = 1000
# config.node_size = 1.0
# config.sl_extra_scaling_steps = 5

start = timer()
x, y, s, t = mstmap.layout_from_lsh_forest(lf, config, True, True, weighted)
lf.clear()
end = timer()
print(end - start)

# Making the first 1028 elements bigger
# sizes = [0.1] * len(x)
# for i in range(1028): sizes[i] = 2.0


for i in range(len(labels)):
    labels[i] = labels[i].lower()
    labels[i] = 'https://cdn.rcsb.org/images/rutgers/' + labels[i][1:3] + '/' + labels[i] + '/' + labels[i] + '.pdb-500.jpg'

for i in range(len(values)):
    print('Writing ' + str(i) + ' ...')
    vals = ss.rankdata(1.0 - np.array(values[i]) / max(values[i])) / len(values[i])

    faerun = Faerun(view='front', coords=False, title=str(i))
    faerun.add_scatter('pdb', { 'x': x, 'y': y, 'c': vals, 'labels': labels }, colormap='rainbow', point_scale=0.75)
    faerun.add_tree('pdb_tree', { 'from': s, 'to': t }, point_helper='pdb')
    faerun.plot('index' + str(i), template='url_image')