import tmap as tm
import numpy as np
from matplotlib import pyplot as plt
from time import sleep


def main():
    n = 25
    edge_list = []

    # Create a random graph
    for i in range(n):
        for j in np.random.randint(0, high=n, size=2):
            edge_list.append([i, j, np.random.rand(1)])

    # Set the initial randomized positioning to True
    # Otherwise, OGDF tends to segfault
    cfg = tm.LayoutConfiguration()
    cfg.fme_randomize = True

    # Compute the layout
    x, y, s, t, gp = tm.layout_from_edge_list(
        n, edge_list, config=cfg, create_mst=False
    )

    # Plot the edges
    for i in range(len(s)):
        plt.plot(
            [x[s[i]], x[t[i]]],
            [y[s[i]], y[t[i]]],
            "k-",
            linewidth=0.5,
            alpha=0.5,
            zorder=1,
        )

    # Plot the vertices
    plt.scatter(x, y, zorder=2)

    plt.savefig("simple_graph.png")


if __name__ == "__main__":
    main()
