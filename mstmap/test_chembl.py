import os
import sys
import pickle

import numpy as np
import pandas as pd
import multiprocessing as mp
import networkx as nx
import scipy.stats as ss

from rdkit.Chem import AllChem, Descriptors
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
f = 512
lf = mstmap.LSHForest(f, 64, store=True)

print('Loading CHEBML data ...')

smiles = []
index = 0
chunk_id = 0
for chunk in pd.read_csv('/media/daenu/Even More Data/chembl_db/chembl.mhfp6', sep=';', header=None, chunksize=200000):
    print(chunk_id)
    if chunk_id > 10: break
    chunk_id += 1
    fps = []

    chunk[2] = chunk[2].apply(convert_the_panda)

    for _, record in chunk.iterrows():
        smiles.append(record[1])
        fps.append(record[2])
        index += 1

    lf.batch_add(fps)

lf.index()

# lf.restore('chembl.dat')
# smiles, vals = pickle.load(open('chembl.pickle', 'rb'))

print("Getting knn graph")

config = mstmap.LayoutConfiguration()
config.k = 10
config.kc = 25
config.sl_scaling_x = 5
config.sl_scaling_y = 25
config.placer = mstmap.Placer.Solar
config.merger = mstmap.Merger.EdgeCover
config.merger_factor = 2.0
config.sl_scaling_type = mstmap.ScalingType.RelativeToAvgLength

x, y, s, t = mstmap.layout_from_lsh_forest(lf, config)

vals = []

for smile in smiles:
    mol = AllChem.MolFromSmiles(smile)
    vals.append(Descriptors.MolLogP(mol))

vals = ss.rankdata(1.0 - np.array(vals) / max(vals)) / len(vals)

lf.store('chembl.dat')
with open('chembl.pickle', 'wb+') as handle:
    pickle.dump((smiles, vals), handle, protocol=pickle.HIGHEST_PROTOCOL)

faerun = Faerun(view='front', coords=False, title='ChEMBL')
faerun.add_scatter('chembl', { 'x': x, 'y': y, 'c': vals, 'labels': smiles }, colormap='rainbow', point_scale=0.75)
faerun.add_tree('chembl_tree', { 'from': s, 'to': t }, point_helper='chembl')
faerun.plot('chembl', template='smiles')