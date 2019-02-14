import os
import sys
import random
import pickle
import numpy as np
from timeit import default_timer as timer
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForest
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
            if i > 2000: break
    pickle.dump(fps, open('fps.dat', 'wb'))
else:
    fps = pickle.load(open('fps.dat', 'rb'))


r = enc.encode('CNCNCNCNC')

# lf_classic = LSHForest(512, 128)
# start = timer()
# for i, e in enumerate(fps):
#     lf_classic.add(i, e)
# end = timer()
# print(end - start)

# start = timer()
# lf_classic.index()
# end = timer()
# print(end - start)

# start = timer()
# for _ in range(100):
#     result = lf_classic.query(r, 5)
# end = timer()
# print(end - start)

# print(result)

lf = mstmap.LSHForest(512, 128)
start = timer()
for i, e in enumerate(fps):
    lf.add(i, mstmap.VectorUint(e))
    # lf.add(i, e)
end = timer()
print(end - start)

start = timer()
lf.index()
end = timer()
print(end - start)

start = timer()
for _ in range(100):
    result = lf.query(mstmap.VectorUint(r), 5)
# result = lf.query(r, 5)
end = timer()

print(end - start)
print(result)
