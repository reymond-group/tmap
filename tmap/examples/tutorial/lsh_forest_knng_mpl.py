from timeit import default_timer as timer
from matplotlib import pyplot as plt
import numpy as np
import tmap as tm


def main():
    """ Main function """

    # Use 128 permutations to create the MinHash
    enc = tm.Minhash(128)
    lf = tm.LSHForest(128)

    d = 10000
    n = 1000

    data = []

    # Generating some random data
    start = timer()
    for i in range(n):
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

    # The configuration for the MST plot
    # Distribute the tree more evenly
    cfg = tm.LayoutConfiguration()
    cfg.sl_scaling_min = 1
    cfg.sl_scaling_max = 1
    cfg.node_size = 1 / 50

    # Construct the k-nearest neighbour graph
    start = timer()
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)
    print(f"layout_from_lsh_forest took {(timer() - start) * 1000}ms.")

    # Plot spanning tree layout
    start = timer()
    for i in range(len(s)):
        plt.plot(
            [x[s[i]], x[t[i]]],
            [y[s[i]], y[t[i]]],
            "r-",
            linewidth=1.0,
            alpha=0.5,
            zorder=1,
        )

    plt.scatter(x, y, s=0.1, zorder=2)
    plt.tight_layout()
    plt.savefig("lsh_forest_knng_mpl.png")
    print(f"Plotting using matplotlib took {(timer() - start) * 1000}ms.")


if __name__ == "__main__":
    main()
