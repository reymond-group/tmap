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


def write_edgelist(path, g):
    with open(path, 'w+') as f:
        for edge in g.edges():
            f.write(str(edge[0]) + ' ' +  str(edge[1]) + ' ' + str(g[edge[0]][edge[1]]['weight']) + '\n')


def get_edge_tuple(g):
    edges = []
    for edge in g.edges():
        edges.append((edge[0], edge[1], g[edge[0]][edge[1]]['weight']))
    return edges

def convert_the_panda(value):
    return mstmap.VectorUint(list(map(int, value.split(','))))

out_name = 'chembl'
f = 2048
lf = mstmap.LSHForest(f, 32, store=True)

print('Loading CHEBML data ...')

smiles = []
index = 0
chunk_id = 0
for chunk in pd.read_csv('/media/daenu/Even More Data/chembl_db/chembl.mhfp6', sep=';', header=None, chunksize=10000):
    print(chunk_id)
    if chunk_id > 5: break
    chunk_id += 1
    fps = []

    chunk[2] = chunk[2].apply(convert_the_panda)

    for _, record in chunk.iterrows():
        smiles.append(record[1])
        fps.append(record[2])
        index += 1

    lf.batch_add(fps)

lf.index()

lf.store('test.tmp')

print("Getting knn graph")

config = mstmap.LayoutConfiguration()
config.k = 10
coords = mstmap.layout_from_lsh_forest(lf, config)

x = coords[0]
y = coords[1]
z = [0] * len(coords[0])

faerun = Faerun(view='front', shader='legacyCircle', coords=False, point_size=0.5, tree_color='#ff0000')
faerun.plot({ 'x': x, 'y': y, 'c': [0] * index, 'smiles': smiles }, colormap='tab10')#, tree=edges)