import os
import sys
import pickle

import numpy as np
import pandas as pd
import scipy.stats as ss

from timeit import default_timer as timer
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from mstmap import mstmap
from operator import itemgetter
from faerun import Faerun
from rdkit.Chem import AllChem, Descriptors


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


f = 128
lf = mstmap.LSHForest(f, 128, store=True)
enc = mstmap.Minhash(f)

print('Loading CHEBML data ...')

smiles = []
index = 0
chunk_id = 0
for chunk in pd.read_csv('/media/daenu/Even More Data/zinc/zinc.mhfp6', sep=';', header=None, chunksize=500000):
    print(chunk_id)
    if chunk_id > 9: break
    chunk_id += 1
    fps = []

    chunk[2] = chunk[2].apply(convert_the_panda)

    for _, record in chunk.iterrows():
        smiles.append(record[0])
        fps.append(record[2])
        index += 1

    start = timer()
    lf.batch_add(fps)
    end = timer()
    print(end - start)

start = timer()
lf.index()
end = timer()

print("Getting knn graph")


config = mstmap.LayoutConfiguration()
config.k = 10
config.kc = 100
config.sl_scaling_x = 5
config.sl_scaling_y = 25
config.placer = mstmap.Placer.Solar
config.merger = mstmap.Merger.EdgeCover
config.merger_factor = 2.0
config.sl_scaling_type = mstmap.ScalingType.RelativeToAvgLength

start = timer()
x, y, s, t = mstmap.layout_from_lsh_forest(lf, config)
lf.store('zinc.dat')
lf.clear()
end = timer()

print(end - start)

vals = []
i = 0
l = len(smiles)
for smile in smiles:
    i += 1
    if l % i == 0: print(i / l)
    mol = AllChem.MolFromSmiles(smile)
    vals.append(Descriptors.MolLogP(mol))

with open('zinc.pickle', 'wb+') as handle:
   pickle.dump((smiles, vals), handle, protocol=pickle.HIGHEST_PROTOCOL)

# smiles, vals = pickle.load(open('zinc.pickle', 'rb'))

vals = ss.rankdata(1.0 - np.array(vals) / max(vals)) / len(vals)

faerun = Faerun(view='front', coords=False, title='ZINC')
faerun.add_scatter('zinc', { 'x': x, 'y': y, 'c': vals, 'labels': smiles }, colormap='rainbow', point_scale=0.75)
faerun.add_tree('zinc_tree', { 'from': s, 'to': t }, point_helper='zinc')

with open('zinc.faerun', 'wb+') as handle:
   pickle.dump(faerun.create_python_data(), handle, protocol=pickle.HIGHEST_PROTOCOL)

faerun.plot('zinc', template='smiles')