import os
import sys
import pickle
import math

import umap
import numpy as np
import pandas as pd
import scipy.stats as ss

from timeit import default_timer as timer
from mnist import MNIST
from progress.bar import Bar
from progress.spinner import Spinner
from mhfp.encoder import MHFPEncoder
from tmap import tmap
from operator import itemgetter
from faerun import Faerun
from annoy import AnnoyIndex
from scipy.spatial.distance import cosine as cosine_distance
from sklearn import preprocessing

def convert_the_panda(value):
    return tmap.VectorUint(list(map(int, value.split(','))))

def convert_the_panda_float(value):
    return list(map(float, value.split(',')))

out_name = 'pdb'
f = 512
enc = tmap.Minhash(136, 42, f)
lf = tmap.LSHForest(f, 128, store=True, file_backed=False)
weighted = True

print('Loading pdb data ...')


ANNOY = AnnoyIndex(136) 


labels = []
values = [[], [], [], [], [], [], []]
index = 0
chunk_id = 0
all_fps = []
list_fps = []
for chunk in pd.read_csv('/media/daenu/Even More Data/pdb/full_fp.csv', sep=',', header=None, chunksize=20000):
    print(chunk_id)
    if chunk_id > 9: break
    chunk_id += 1

    for _, record in chunk.iterrows():
        labels.append(record[0])
        for i in range(7):
            values[i].append(record[137 + i])
        
        list_fps.append(tmap.VectorFloat([float(i) for i in record[1:137]]))
        all_fps.append(tmap.VectorFloat([float(i) for i in record[1:137]]))

        ANNOY.add_item(index, record[1:137].values)

        index += 1

ANNOY.build(10)

start = timer()
lf.batch_add(enc.batch_from_weight_array(all_fps))
end = timer()
print(end - start)

# lf.store('pdb.dat')
# lf.restore('pdb.dat')
lf.index()


# ANNOY_GRAPH = []
# for i in range(len(all_fps)):
#     if i % 1000 == 0: print(i / len(all_fps))
#     for j in ANNOY.get_nns_by_item(i, 10):
#         ANNOY_GRAPH.append((i, j, cosine_distance(list_fps[i], list_fps[j])))


# with open('pdb.pickle', 'wb+') as handle:
#     pickle.dump((labels, values), handle, protocol=pickle.HIGHEST_PROTOCOL)

labels, values = pickle.load(open('pdb.pickle', 'rb'))

print("Getting knn graph")

config = tmap.LayoutConfiguration()
config.k = 50
config.kc = 50
config.sl_scaling_min = 1.0
config.sl_scaling_max = 1.0
config.sl_repeats = 1
config.sl_extra_scaling_steps = 2 # The higher, the sparser the tree
config.placer = tmap.Placer.Barycenter
config.merger = tmap.Merger.LocalBiconnected
config.merger_factor = 2.0
config.merger_adjustment = 0
config.fme_iterations = 100
config.fme_precision = 4
config.sl_scaling_type = tmap.ScalingType.RelativeToDrawing
config.node_size = 1 / 60
config.mmm_repeats = 1

values.append(np.array(values[2]) / np.array(values[0]))
print(values[7])

start = timer()
x, y, s, t, _ = tmap.layout_from_lsh_forest(lf, config, True, True, weighted)
# x, y, s, t, _ = tmap.layout_from_edge_list(len(list_fps), ANNOY_GRAPH, config)
lf.clear()
end = timer()
print(end - start)

# Making the first 1028 elements bigger
# sizes = [0.1] * len(x)
# for i in range(1028): sizes[i] = 2.0

# get pos neg ratio

# REDUCER = umap.UMAP(n_neighbors=10)
# start = timer()
# COORDS = REDUCER.fit_transform(all_fps)
# end = timer()
# print(end - start)

# x = []
# y = []
# for coord in COORDS:
#     x.append(coord[0])
#     y.append(coord[1])

# for i in range(len(labels)):
#     labels[i] = labels[i].lower()
    # labels[i] = 'https://cdn.rcsb.org/images/rutgers/' + labels[i][1:3] + '/' + labels[i] + '/' + labels[i] + '.pdb-500.jpg'

for i in range(len(values)):
    print('Writing ' + str(i) + ' ...')
    vals = ss.rankdata(np.array(values[i]) / max(values[i])) / len(values[i])
    l = math.floor(len(vals) / 5)

    sorted_values = sorted(values[i])
    print(min(values[i]), sorted_values[l * 2], sorted_values[l * 3], sorted_values[l * 4], max(values[i]))

    faerun = Faerun(view='front', coords=False, title=str(i))
    faerun.add_scatter('pdb', { 'x': x, 'y': y, 'c': vals, 'labels': labels }, colormap='rainbow', point_scale=2.0, max_point_size=20)
    faerun.add_tree('pdb_tree', { 'from': s, 'to': t }, point_helper='pdb', color='#555555')
    faerun.plot('index_experimental_tree' + str(i), template='url_image')

    with open('pdb_experimental' + str(i) + '.faerun', 'wb+') as handle:
        pickle.dump(faerun.create_python_data(), handle, protocol=pickle.HIGHEST_PROTOCOL)

#     sys.exit()
