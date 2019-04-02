import os
import sys
import pickle
import umap

import numpy as np
import pandas as pd
import multiprocessing as mp
import networkx as nx
import scipy.stats as ss

from rdkit.Chem import AllChem, Descriptors, DataStructs
from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from mstmap import mstmap
from operator import itemgetter
from faerun import Faerun

reducer = umap.UMAP(min_dist=0.001)


def convert_the_panda(value):
    return list(map(int, value.split(',')))

print('Loading CHEBML data ...')

smiles = []
index = 0
chunk_id = 0
chembl_id = []
activity = []
target_class = []
fps = []

for chunk in pd.read_csv('/media/daenu/Even More Data/chembl_db/chembl_targets.ecfp4', sep=';', header=None, chunksize=1000):
    chunk_id += 1
    print(chunk_id)
    if chunk_id > 10: break

    chunk[6] = chunk[6].apply(convert_the_panda)

    for _, record in chunk.iterrows():
        chembl_id.append(record[0])
        smiles.append(record[3])
        activity.append(record[4])
        fps.append(record[6])
        tc = record[5].split('  ')[0]
        if tc == 'enzyme' and len(record[5].split('  ')) > 1: tc = record[5].split('  ')[1]
        target_class.append(tc)
        index += 1





tmp = pd.DataFrame({ 'target_class': target_class })

# query = 22
# query = 222
query = 2222

distances = []
for fp in fps:
    distances.append(1.0 - np.bitwise_and(fp, fps[query]).sum() / np.bitwise_or(fp, fps[query]).sum())


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




start = timer()
coords = reducer.fit_transform(fps)
end = timer()
print(end - start)

x = []
y = []
for t in coords:
    x.append(t[0])
    y.append(t[1])


activity = np.array(activity)
activity = np.maximum(0.0, activity)
activity = np.minimum(100.0, activity)
activity = 10.0 - activity

legend_labels = [(0, 'Cytochrome p450'), (1, 'Other Enzyme'), (2, 'Epigenetic Regulator'), (3, 'Ion Channel'), (4, 'Kinase'), 
                (5, 'Membrane Receptor'), (6, 'Protease'), (8, 'Transcription Factor'), (9, 'Transporter'), (7, 'Other')]
# vals = [int(target_class_map[x]) for x in target_class]
# vals = ss.rankdata(1.0 - np.array(distances) / max(distances)) / len(distances)
vals = 1 - np.array(distances)

faerun = Faerun(view='front', coords=False, title='ChEMBL')
faerun.add_scatter('chembl', { 'x': x, 'y': y, 'c': vals, 'labels': smiles }, #colormap='tab10', 
                point_scale=0.5, max_point_size=10, has_legend=True)#, categorical=True, legend_labels=legend_labels)
faerun.plot('chembl_umap', template='smiles')