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

    # Use batch_from_binary_array to encode the data
    start = timer()
    data = enc.batch_from_binary_array(data)
    print(f"Encoding the data took {(timer() - start) * 1000}ms.")

    # Use batch_add to parallelize the insertion of the arrays
    start = timer()
    lf.batch_add(data)
    print(f"Adding the data took {(timer() - start) * 1000}ms.")

    # Index the added data
    start = timer()
    lf.index()
    print(f"Indexing took {(timer() - start) * 1000}ms.")

    # Find the 10 nearest neighbors of the first entry
    start = timer()
    _ = lf.query_linear_scan_by_id(0, 10)
    print(f"The kNN search took {(timer() - start) * 1000}ms.")


if __name__ == "__main__":
    main()
