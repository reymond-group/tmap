import tmap as tm
import numpy as np
from matplotlib import pyplot as plt


def main():
    """ Main function """

    n = 10
    edge_list = []
    weights = {}

    # Create a random graph
    for i in range(n):
        for j in np.random.randint(0, high=n, size=2):
            # Do not add parallel edges here, to be sure
            # to have the right weight later
            if i in weights and j in weights[i] or j in weights and i in weights[j]:
                continue

            weight = np.random.rand(1)
            edge_list.append([i, j, weight])

            # Store the weights in 2d map for easy access
            if i not in weights:
                weights[i] = {}
            if j not in weights:
                weights[j] = {}

            # Invert weights to make lower ones more visible in the plot
            weights[i][j] = 1.0 - weight
            weights[j][i] = 1.0 - weight

    # Compute the layout
    x, y, s, t, _ = tm.layout_from_edge_list(n, edge_list, create_mst=False)
    x_mst, y_mst, s_mst, t_mst, _ = tm.layout_from_edge_list(
        n, edge_list, create_mst=True
    )

    _, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)

    # Plot graph layout with spanning tree superimposed in red
    for i in range(len(s)):
        ax1.plot(
            [x[s[i]], x[t[i]]],
            [y[s[i]], y[t[i]]],
            "k-",
            linewidth=weights[s[i]][t[i]],
            alpha=0.5,
            zorder=1,
        )

    for i in range(len(s_mst)):
        ax1.plot(
            [x[s_mst[i]], x[t_mst[i]]],
            [y[s_mst[i]], y[t_mst[i]]],
            "r-",
            linewidth=weights[s_mst[i]][t_mst[i]],
            alpha=0.5,
            zorder=2,
        )

    ax1.scatter(x, y, zorder=3)

    # Plot spanning tree layout
    for i in range(len(s_mst)):
        ax2.plot(
            [x_mst[s_mst[i]], x_mst[t_mst[i]]],
            [y_mst[s_mst[i]], y_mst[t_mst[i]]],
            "r-",
            linewidth=weights[s_mst[i]][t_mst[i]],
            alpha=0.5,
            zorder=1,
        )

    ax2.scatter(x_mst, y_mst, zorder=2)

    plt.tight_layout()
    plt.savefig("spanning_tree.png")


if __name__ == "__main__":
    main()
