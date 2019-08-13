from timeit import default_timer as timer

import numpy as np
import tmap as tm


def main():
    """ Main function """

    # Use 128 permutations to create the MinHash
    enc = tm.Minhash(1024)
    lf = tm.LSHForest(128, file_backed=True)

    # d = 1000
    # n = 1000000
    d = 10000
    n = 1000

    # Generating some random data
    start = timer()
    for _ in range(n):
        # data.append(tm.VectorUint(np.random.randint(0, high=2, size=d)))
        lf.add(
            enc.from_sparse_binary_array(
                tm.VectorUint(np.random.randint(0, high=2, size=d))
            )
        )

    print(f"Generating the data took {(timer() - start) * 1000}ms.")

    # Index the added data
    start = timer()
    lf.index()
    print(f"Indexing took {(timer() - start) * 1000}ms.")

    # Find the 10 nearest neighbors of the first entry
    start = timer()
    knng_from = tm.VectorUint()
    knng_to = tm.VectorUint()
    knng_weight = tm.VectorFloat()

    _ = lf.get_knn_graph(knng_from, knng_to, knng_weight, 10)
    print(f"The kNN search took {(timer() - start) * 1000}ms.")


if __name__ == "__main__":
    main()
