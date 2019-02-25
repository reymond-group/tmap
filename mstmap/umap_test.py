import umap

reducer = umap.UMAP()

import os
import sys
import pickle

import numpy as np
import pandas as pd
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


def convert_the_panda(value):
    return list(map(int, value.split(',')))

out_name = 'chembl'

print('Loading CHEBML data ...')

smiles = []
fps = []
index = 0
chunk_id = 0
for chunk in pd.read_csv('/media/daenu/Even More Data/zinc/zinc.ecfp4', sep=';', header=None, chunksize=30000):
    print(chunk_id)
    if chunk_id > 9: break
    chunk_id += 1
    

    chunk[2] = chunk[2].apply(convert_the_panda)

    for _, record in chunk.iterrows():
        smiles.append(record[0])
        fps.append(record[2])
        index += 1


coords = reducer.fit_transform(fps)

x = []
y = []
for t in coords:
    x.append(t[0])
    y.append(t[1])

values = np.array(list(map(len, smiles)))
values = values / max(values)

faerun = Faerun(view='front', shader='legacyCircle', coords=False, point_size=0.5, tree_color='#ff0000')
faerun.plot({ 'x': x, 'y': y, 'c': values, 'smiles': smiles }, colormap='rainbow')#, tree=edges)

