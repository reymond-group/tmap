import os
import sys
import random
import pickle
import numpy as np
from timeit import default_timer as timer
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForest, LSHForestHelper
from mstmap import mstmap
from rdkit.Chem import AllChem

enc = MHFPEncoder(512)


fps = []

if not os.path.isfile('fps.dat'):
    with open('drugbank.smi', 'r') as f:
        i = 0
        for line in f:
            smiles = line.split()[0].strip()
            mol = AllChem.MolFromSmiles(smiles)
            if mol:
                fps.append(enc.encode_mol(mol))
            i += 1
            # if i > 2000: break
    pickle.dump(fps, open('fps.dat', 'wb'))
else:
    fps = pickle.load(open('fps.dat', 'rb'))


r = fps[0] # enc.encode('CNCNCNCNC')


lf_classic = LSHForestHelper(512, 128)

start = timer()
for i, e in enumerate(fps):
    lf_classic.add(i, e)
end = timer()
print(end - start)

start = timer()
lf_classic.index()
end = timer()
print(end - start)

start = timer()
for _ in range(100):
    result = lf_classic.query(r, 5, fps)

for i in result:
    print(enc.distance(r, fps[i]))
end = timer()
print(end - start)

print(result)


fps_e = []
for fp in fps:
    fps_e.append(mstmap.VectorUint(fp))

lf = mstmap.LSHForest(512, 128, store=True)
start = timer()
lf.batch_add(fps_e)
end = timer()
print(end - start)


start = timer()
lf.index()
end = timer()
print(end - start)

lf.store('test.bin')
lf.restore('test.bin')

start = timer()
# result_ls = lf.batch_query([mstmap.VectorUint(r)] * 100, 5)
result_ls = lf.query_linear_scan(mstmap.VectorUint(r), 5)
end = timer()
print(end - start)

print(result_ls)
for i in result_ls:
    print(i)

sys.exit()

print("Excluding example")
print(lf.query_linear_scan(mstmap.VectorUint(r), 5))
print(lf.query_linear_scan_exclude(mstmap.VectorUint(r), 5, mstmap.VectorUint([0])))

print("Query by id")
print(lf.query_linear_scan(mstmap.VectorUint(fps[0]), 5))
print(lf.query_linear_scan_by_id(0, 5))

print("kNN graph")
start = timer()
n = lf.get_knn_graph(100)
end = timer()
print(end - start)

print("Get distance by id")
print(lf.get_distance(mstmap.VectorUint(fps[2]), mstmap.VectorUint(fps[6])))
print(lf.get_distance_by_id(2, 6))

# for e in n:
#     print(e[0], e[1], e[2])
# print(len(n))

# print(result_ls)


mh = mstmap.Minhash()
fp_a = mh.from_binary_array(mstmap.VectorUchar([1, 1, 1, 0, 0, 1, 1, 0, 0]))
mh = mstmap.Minhash()
fp_b = mh.from_sparse_binary_array(mstmap.VectorUint([0, 1, 2, 5, 6]))
print(fp_a == fp_b)

mh_classic = MHFPEncoder()

start = timer()
for i in range(10000):
    a = mh.from_string_array(['hello'.encode('utf8'), 'world'.encode('utf8')])
    b = mh.from_string_array(['hello'.encode('utf8'), 'test'.encode('utf8')])
end = timer()
print(end - start)


mh = mstmap.Minhash()
start = timer()
for i in range(10000):
    a = mh.from_string_array(['hello', 'world'])
    b = mh.from_string_array(['hello', 'test'])
end = timer()
print(end - start)
print(mh.get_distance(a, b))

mh = mstmap.Minhash()
start = timer()
a = mh.batch_from_string_array([['hello', 'world']] * 10000)
b = mh.batch_from_string_array([['hello', 'test']] * 10000)
end = timer()
print(end - start)
# print(mh.get_distance(a, b))


print('\nWeighted...')
lf = mstmap.LSHForest(8, store=True)
mh = mstmap.Minhash(8)
a = mh.from_weight_array(mstmap.VectorFloat([0.1, 0.2, 0.6, 0, 0.2, 0.9, 0.6, 0]))
b = mh.from_weight_array(mstmap.VectorFloat([0.9, 0.2, 0.6, 0, 0.2, 0.9, 0.6, 0]))
c = mh.from_weight_array(mstmap.VectorFloat([0.9, 0.5, 0.6, 0, 0.2, 0.9, 0.6, 0]))

dist = mh.get_weighted_distance(a, b)
print(dist)
dist = mh.get_weighted_distance(a, c)
print(dist)
dist = mh.get_weighted_distance(b, c)
print(dist)

# dist = lf.get_weighted_distance(a, b)


lf.add(0, a)
lf.add(1, b)
lf.add(2, c)
lf.index()

print(lf.query_linear_scan(c, 1, 10, True))
