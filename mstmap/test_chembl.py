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
lf = mstmap.LSHForest(f, 128, store=True)

print('Loading CHEBML data ...')

smiles = []
index = 0
chunk_id = 0
chembl_id = []
activity = []
target_class = []

# for chunk in pd.read_csv('/media/daenu/Even More Data/chembl_db/chembl_targets.mhfp6', sep=';', header=None, chunksize=100000):
#     print(chunk_id)
#     # if chunk_id > 10: break
#     chunk_id += 1
#     fps = []

#     chunk[6] = chunk[6].apply(convert_the_panda)

#     for _, record in chunk.iterrows():
#         chembl_id.append(record[0])
#         smiles.append(record[3])
#         activity.append(record[4])
#         fps.append(record[6])
#         tc = record[5].split('  ')[0]
#         if tc == 'enzyme' and len(record[5].split('  ')) > 1: tc = record[5].split('  ')[1]
#         target_class.append(tc)
#         index += 1

#     lf.batch_add(fps)


smiles, target_class, activity, chembl_id = pickle.load(open('chembl.pickle', 'rb'))
lf.restore('chembl.dat')
lf.index()

tmp = pd.DataFrame({ 'target_class': target_class })
print(tmp)
print(tmp['target_class'].value_counts())
# lf.store('chembl.dat')
# with open('chembl.pickle', 'wb+') as handle:
#     pickle.dump((smiles, target_class, activity, chembl_id), handle, protocol=pickle.HIGHEST_PROTOCOL)

# query = 22
# query = 222
query = 2222
distances = lf.get_all_distances(lf.get_hash(query))

target_class_map = dict([(y,x+1) for x, y in enumerate(sorted(set(target_class)))])
classes = ['enzyme', 'kinase', 'protease', 'cytochrome p450', 'ion channel', 'transporter', 'transcription factor', 'membrane receptor', 'epigenetic regulator']
i = 0
for key, value in target_class_map.items():
    if key not in classes:
        target_class_map[key] = 7
    else:
        target_class_map[key] = i
        i += 1
        if i is 7: i = 8

print(target_class_map)

print("Getting knn graph")

# Median - EdgeCover
# Barycenter - Solar

for i in [1]:
    for j in [0]:
        config = mstmap.LayoutConfiguration()
        config.k = 20
        config.kc = 20
        config.sl_scaling_min = 1.0
        config.sl_scaling_max = 1.0
        config.sl_repeats = 1
        config.sl_extra_scaling_steps = 2 # The higher, the sparser the tree
        config.placer = mstmap.Placer.Barycenter
        config.merger = mstmap.Merger.LocalBiconnected
        config.merger_factor = 2.0
        config.merger_adjustment = 0
        config.fme_iterations = 100
        config.fme_precision = 4
        config.sl_scaling_type = mstmap.ScalingType.RelativeToDrawing
        config.node_size = 1 / 65
        config.mmm_repeats = 1


        start = timer()
        x, y, s, t, _ = mstmap.layout_from_lsh_forest(lf, config)
        end = timer()
        print(end - start)

        lf.clear()

        activity = np.array(activity)
        activity = np.maximum(0.0, activity)
        activity = np.minimum(100.0, activity)
        activity = 10.0 - activity

        legend_labels = [(0, 'Cytochrome p450'), (1, 'Other Enzyme'), (2, 'Epigenetic Regulator'), (3, 'Ion Channel'), (4, 'Kinase'), 
                        (5, 'Membrane Receptor'), (6, 'Protease'), (8, 'Transcription Factor'), (9, 'Transporter'), (7, 'Other')]

        vals = [int(target_class_map[x]) for x in target_class]
        # vals = ss.rankdata(1.0 - np.array(distances) / max(distances)) / len(distances)
        # vals = 1 - np.array(distances)

        faerun = Faerun(view='front', coords=False, title='ChEMBL')
        faerun.add_scatter('chembl', { 'x': x, 'y': y, 'c': vals, 'labels': smiles }, colormap='tab10', 
                        point_scale=2.5, max_point_size=10, has_legend=True, categorical=True, legend_labels=legend_labels)
        faerun.add_tree('chembl_tree', { 'from': s, 'to': t }, point_helper='chembl', color='#222222')
        faerun.plot('chembl' + str(i) + str(j), template='smiles')