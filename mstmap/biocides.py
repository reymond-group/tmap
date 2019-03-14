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
    return mstmap.VectorUint(list(map(int, value.split(','))))

def convert_the_panda_float(value):
    return list(map(float, value.split(',')))

out_name = 'chembl'
f = 512
# enc = mstmap.Minhash(42, 42, f)
lf = mstmap.LSHForest(f, 128, store=True, file_backed=False)
weighted = False

print('Loading CHEBML data ...')

smiles = []
values = [[], [], [], [], [], [], []]
index = 0
chunk_id = 0
# for chunk in pd.read_csv('/media/daenu/Even More Data/biocides/vs_chembl/biocides_chembl.mhfp6', sep=';', header=None, chunksize=2000):
#     print(chunk_id)
#     if chunk_id > 9: break
#     chunk_id += 1
#     fps = []

#     chunk[1] = chunk[1].apply(convert_the_panda)
#     chunk[2] = chunk[2].apply(convert_the_panda_float)

#     for _, record in chunk.iterrows():
#         smiles.append(record[0])
#         for i in range(len(record[2])):
#             values[i].append(record[2][i])
        
#         fps.append(record[1])
#         index += 1

#     start = timer()
#     lf.batch_add(fps)
#     end = timer()
#     print(end - start)

# start = timer()
# lf.index()
# end = timer()

# lf.store('biocides_chembl.dat')
# with open('biocides_chembl.pickle', 'wb+') as handle:
#     pickle.dump((smiles, values), handle, protocol=pickle.HIGHEST_PROTOCOL)


lf.restore('biocides_chembl.dat')
smiles, values = pickle.load(open('biocides_chembl.pickle', 'rb'))

print("Getting knn graph")

config = mstmap.LayoutConfiguration()
config.k = 10
config.kc = 10
config.sl_scaling_x = 5
config.sl_scaling_y = 25
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

with open('biocides_chembl_coords.pickle', 'wb+') as handle:
    pickle.dump((list(x), list(y), list(s), list(t)), handle, protocol=pickle.HIGHEST_PROTOCOL)


species = ['Aspergillus_Brasiliensis', 'Pseudomonas_Aeruginosa', 'Staphylococcus_Aureus', 'Aspergillus_Pseudomonas_Staphylococcus', 'Aspergillus_Pseudomonas', 'Aspergillus_Sraphylococcus', 'Pseudomonas_Staphylococcus']

# Making the first 1028 elements bigger
# sizes = [0.1] * len(x)
# for i in range(1028): sizes[i] = 2.0

for i in range(len(species)):
    print('Writing ' + species[i] + ' ...')
    vals = 1.0 - np.array(values[i]) / max(values[i])
    
    # Split into the two sets

    faerun = Faerun(view='front', coords=False, title=species[i])
    faerun.add_scatter('all', { 'x': x, 'y': y, 'c': vals, 'labels': smiles }, colormap='rainbow', point_scale=0.25, max_point_size=25, shader='circle')
    faerun.add_tree('all_tree', { 'from': s, 'to': t }, point_helper='all', color='#222222')

    faerun.add_scatter('biocides', { 'x': x[:1028], 'y': y[:1028], 'c': vals[:1028], 'labels': smiles[:1028] }, colormap='rainbow', point_scale=5.0, shader='circle')
    faerun.plot(species[i], template='smiles')