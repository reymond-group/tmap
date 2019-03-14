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
    return mstmap.VectorFloat(list(map(int, value.split(','))))

out_name = 'chembl'
f = 256
lf = mstmap.LSHForest(f, 128, store=True)
enc = mstmap.Minhash(42, 42, f)

print('Loading CHEBML data ...')

smiles = []
values = []
index = 0
chunk_id = 0
for chunk in pd.read_csv('/media/daenu/Even More Data/chembl_db/chembl.mqn', sep=';', header=None, chunksize=200):
    print(chunk_id)
    if chunk_id > 9: break
    chunk_id += 1
    fps = []

    chunk[2] = chunk[2].apply(convert_the_panda)

    for _, record in chunk.iterrows():
        if record[2][11] < 50:
            smiles.append(record[1])
            values.append(record[2][11])
            fps.append(record[2])
            index += 1

    start = timer()
    lf.batch_add(enc.batch_from_weight_array(fps))
    end = timer()
    print(end - start)

start = timer()
lf.index()
end = timer()

print("Getting knn graph")

# coords = mstmap.igraph_layout_from_lsh_forest(lf)
# print(coords)

config = mstmap.LayoutConfiguration()
config.k = 10
config.sl_scaling_x = 5
config.sl_scaling_y = 50
config.merger_factor = 2.0
config.sl_scaling_type = mstmap.ScalingType.RelativeToDesiredLength
# config.fme_iterations = 10
# config.node_size = 1.0
# config.sl_extra_scaling_steps = 5

start = timer()
coords = mstmap.layout_from_lsh_forest(lf, config, True, True)
# coords = lf.get_layout(config)
lf.clear()
end = timer()

# input("Press Enter to continue...")

print(end - start)

x = coords[0]
y = coords[1]


f = mstmap.VectorUint()
t = mstmap.VectorUint()
w = mstmap.VectorFloat()
g = lf.get_knn_graph(f, t, w, 10)

# f, t = mstmap.mst_from_lsh_forest(lf, 10)

# with open('test.edges', 'w+') as f_out:
#     for i in range(len(t)):
#         f_out.write(str(f[i]) + ' ' + str(t[i]) + '\n')


# values = np.array(list(map(len, smiles)))
values = np.array(values) / max(values)

faerun = Faerun(view='front', shader='legacyCircle', coords=False, point_size=0.5, tree_color='#ff0000')
faerun.plot({ 'x': x, 'y': y, 'c': values, 'smiles': smiles }, colormap='rainbow')#, tree=edges)