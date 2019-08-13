from timeit import default_timer as timer

import numpy as np
import tmap as tm


def main():
    """ Main function """

    # Use 128 permutations to create the MinHash
    enc = tm.Minhash(128)
    lf = tm.LSHForest(128)

    d = 1000
    n = 10000

    data = []

    # Generating some random data
    start = timer()
    for _ in range(n):
        data.append(tm.VectorUchar(np.random.randint(0, high=2, size=d)))
    print(f"Generating the data took {(timer() - start) * 1000}ms.")

    # Use batch_add to parallelize the insertion of the arrays
    start = timer()
    lf.batch_add(enc.batch_from_binary_array(data))
    print(f"Adding the data took {(timer() - start) * 1000}ms.")

    # Index the added data
    start = timer()
    lf.index()
    print(f"Indexing took {(timer() - start) * 1000}ms.")

    # Construct the k-nearest neighbour graph
    start = timer()
    knng_from = tm.VectorUint()
    knng_to = tm.VectorUint()
    knng_weight = tm.VectorFloat()

    _ = lf.get_knn_graph(knng_from, knng_to, knng_weight, 10)
    print(f"The kNN search took {(timer() - start) * 1000}ms.")


if __name__ == "__main__":
    main()
