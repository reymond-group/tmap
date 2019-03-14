import umap

reducer = umap.UMAP()

import os
import sys
import pickle

import numpy as np
import pandas as pd
import multiprocessing as mp
import networkx as nx
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
    return list(map(int, value.split(',')))

out_name = 'chembl'

print('Loading CHEBML data ...')

lf = mstmap.LSHForest(128, 128, store=True)
enc = MHFPEncoder(128)

# smiles = []
# fps = []
# index = 0
# chunk_id = 0
# for chunk in pd.read_csv('/media/daenu/Even More Data/zinc/zinc.ecfp4', sep=';', header=None, chunksize=1000):
#     print(chunk_id)
#     if chunk_id > 9: break
#     chunk_id += 1
    

#     chunk[2] = chunk[2].apply(convert_the_panda)

#     for _, record in chunk.iterrows():
#         smiles.append(record[0])
#         fps.append(record[2])
#         index += 1

# coords = reducer.fit_transform(fps)


# index = 0
# chunk_id = 0
# fps = []
# smiles = []
# for chunk in pd.read_csv('/media/daenu/Even More Data/zinc/zinc.mhfp6', sep=';', header=None, chunksize=1000):
#     print(chunk_id)
#     if chunk_id > 9: break
#     chunk_id += 1

#     chunk[2] = chunk[2].apply(convert_the_panda)

#     for _, record in chunk.iterrows():
#         smiles.append(record[0])
#         fps.append(mstmap.VectorUint(record[2]))
#         index += 1


index = 0
chunk_id = 0
fps = []
labels = []
values = [[], [], [], [], [], [], []]
for chunk in pd.read_csv('/media/daenu/Even More Data/pdb/full_fp.csv', sep=',', header=None, chunksize=20000):
    print(chunk_id)
    if chunk_id > 9: break
    chunk_id += 1

    for _, record in chunk.iterrows():
        labels.append(record[0])
        for i in range(7):
            values[i].append(record[137 + i])
        
        fps.append([float(i) for i in record[1:137]])
        index += 1

    start = timer()
    end = timer()
    print(end - start)



# lf.batch_add(fps)
# lf.index()
# e_s, e_t = mstmap.mst_from_lsh_forest(lf, 10)
# x, y = mstmap.layout_from_lsh_forest(lf)

coords = reducer.fit_transform(fps)


x = []
y = []
for t in coords:
    x.append(t[0])
    y.append(t[1])

# values = np.array(list(map(len, smiles)))
# values = values / max(values)

for i in range(len(labels)):
    labels[i] = labels[i].lower()
    labels[i] = 'https://cdn.rcsb.org/images/rutgers/' + labels[i][1:3] + '/' + labels[i] + '/' + labels[i] + '.pdb-500.jpg'

faerun = Faerun(view='back', coords=False)
vals = ss.rankdata(1.0 - np.array(values[0]) / max(values[0])) / len(values[0])
faerun.add_scatter('pdb', { 'x': x, 'y': y, 'c': vals, 'labels': labels }, colormap='rainbow', point_scale=0.25)
faerun.plot(template='url_image')

