import numpy as np
from tmap import tmap

np.random.seed(42)

D = 136 #11757

enc = tmap.Minhash(D)
lf = tmap.LSHForest(D, 8, store=True)

data = []
for i in range(100000):
    data.append(tmap.VectorFloat(np.random.uniform(size=D)))

hashes = enc.batch_from_weight_array_experimental(data)
lf.batch_add(hashes)
lf.index()

result = lf.query_linear_scan_by_id(66, 10)
# print(hashes)
