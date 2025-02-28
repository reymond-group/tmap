Getting Started
---------------

Installation
^^^^^^^^^^^^
TMAP can be installed using the conda package manager that
is distributed with miniconda (or anaconda).

.. code-block:: bash

   conda install tmap

The module is then best imported using a shorter identifier.

.. code-block:: python

   import tmap as tm

Laying out a Simple Graph
^^^^^^^^^^^^^^^^^^^^^^^^^
Even though TMAP is mainly targeted at tasks consisting of
laying out very large data sets, the simplest usage example
is laying out a graph.


.. code-block:: python

    import tmap as tm
    import numpy as np
    from matplotlib import pyplot as plt

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
    x, y, s, t, gp = tm.layout_from_edge_list(n, edge_list, config=cfg,
                                              create_mst=False)

    # Plot the edges
    for i in range(len(s)):
        plt.plot([x[s[i]], x[t[i]]], [y[s[i]], y[t[i]]], 'k-',
                 linewidth=0.5, alpha=0.5, zorder=1)

    # Plot the vertices
    plt.scatter(x, y, zorder=2)

    plt.savefig('simple_graph.png')

.. image:: _static/simple_graph.png
   :alt: A simple graph layed out using TMAP.

When laying out large graphs, it might be useful to discard
some edges in order to create a more interpretable and visually
pleasing layout. This is achieved using the (default) argument
:obj:`create_mst=True`. Following, this is exemplified on a small
graph.

.. code-block:: python

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


    # Set the initial randomized positioning to True
    # Otherwise, OGDF tends to segfault
    cfg = tm.LayoutConfiguration()
    cfg.fme_randomize = True

    # Compute the layout
    x, y, s, t, _ = tm.layout_from_edge_list(n, edge_list, config=cfg,
                                             create_mst=False)
    x_mst, y_mst, s_mst, t_mst, _ = tm.layout_from_edge_list(n, edge_list,
                                                             create_mst=True)

    _, (ax1, ax2) = plt.subplots(ncols=2, sharey=True)

    # Plot graph layout with spanning tree superimposed in red
    for i in range(len(s)):
        ax1.plot([x[s[i]], x[t[i]]], [y[s[i]], y[t[i]]], 'k-',
                 linewidth=weights[s[i]][t[i]], alpha=0.5, zorder=1)

    for i in range(len(s_mst)):
        ax1.plot([x[s_mst[i]], x[t_mst[i]]], [y[s_mst[i]], y[t_mst[i]]], 'r-',
                 linewidth=weights[s_mst[i]][t_mst[i]], alpha=0.5, zorder=2)

    ax1.scatter(x, y, zorder=3)


    # Plot spanning tree layout
    for i in range(len(s_mst)):
        ax2.plot([x_mst[s_mst[i]], x_mst[t_mst[i]]], [y_mst[s_mst[i]], y_mst[t_mst[i]]], 'r-',
                 linewidth=weights[s_mst[i]][t_mst[i]], alpha=0.5, zorder=1)

    ax2.scatter(x_mst, y_mst, zorder=2)

    plt.tight_layout()
    plt.savefig('spanning_tree.png')

.. image:: _static/spanning_tree.png
   :alt: The minimum spanning tree of a simple graph layed out using TMAP.

On a highly connected graph with 1000 vertices, the advantages
of the tree visualizaton method applied by TMAP become obvious.

.. image:: _static/spanning_tree_big.png
   :alt: The minimum spanning tree of a big simple graph layed out using TMAP.

There are a wide array of options to tune the final tree layout
to your linking. See :ref:`layout-doc` for the descriptions of all
available parameters.


MinHash
^^^^^^^
In order to enable the visualization of larger data sets, it is
necessary to speed up the k-nearest neighbor graph generation. While
in general, any approach can be used to create this nearest neighbor
graph (see Laying out a Simple Graph), TMAP provides a built-in LSH Forest
data structure, which enables extremely fast k-nearest neighbor queries.


In order to index data in the LSH forest data structure, it has to be
hashed using a locality sensitive scheme such as MinHash.

TMAP includes the two classes :obj:`MinHash` and :obj:`LSHForest` for
fast k-nearest neighbor search.

The following example shows how to use the :obj:`MinHash` class to estimate
Jaccard distances.

.. code-block:: python

    import tmap as tm

    enc = tm.Minhash()

    mh_a = enc.from_binary_array(tm.VectorUchar([1, 1, 1, 1, 0, 1, 0, 1, 1, 0]))
    mh_b = enc.from_binary_array(tm.VectorUchar([1, 0, 1, 1, 0, 1, 1, 0, 1, 0]))
    mh_c = enc.from_binary_array(tm.VectorUchar([1, 0, 1, 1, 1, 1, 1, 0, 1, 0]))

    dist_a_b = enc.get_distance(mh_a, mh_b)
    dist_b_c = enc.get_distance(mh_b, mh_c)

    print(dist_a_b)
    print(dist_b_c)

>>> 0.390625
>>> 0.140625


An in-depth explanation of MinHash can be found in
`this <https://www.youtube.com/watch?v=96WOGPUgMfw>`_ video by Jeffry D
Ullman.


:obj:`Minhash` also supports encoding strings, indexed binary arrays, and
:obj:`int` and :obj:`float` weighted arrays. See :ref:`minhash-doc` for details.

LSH Forest
^^^^^^^^^^
The hashes generated by :obj:`Minhash` can be indexed using :obj:`LSHForest`
for fast k-nearest neighbor retreival.

.. code-block:: python

    from timeit import default_timer as timer

    import numpy as np
    import tmap as tm

    # Use 128 permutations to create the MinHash
    enc = tm.Minhash(128)
    lf = tm.LSHForest(128)

    d = 1000
    n = 1000000

    data = []

    # Generating some random data
    start = timer()
    for i in range(n):
        data.append(tm.VectorUchar(np.random.randint(0, high=2, size=d)))
    print(f'Generating the data took {(timer() - start) * 1000}ms.')

    # Use batch_add to parallelize the insertion of the arrays
    start = timer()
    lf.batch_add(enc.batch_from_binary_array(data))
    print(f'Adding the data took {(timer() - start) * 1000}ms.')

    # Index the added data
    start = timer()
    lf.index()
    print(f'Indexing took {(timer() - start) * 1000}ms.')

    # Find the 10 nearest neighbors of the first entry
    start = timer()
    result = lf.query_linear_scan_by_id(0, 10)
    print(f'The kNN search took {(timer() - start) * 1000}ms.')

>>> Generating the data took 118498.04133399994ms.
>>> Adding the data took 55051.067827000224ms.
>>> Indexing took 2059.1810410005564ms.
>>> The kNN search took 0.32151699997484684ms.

After indexing the data, the 10 nearest neighbor search on a million
1,000-dimensional vectors took ~0.32ms. In addition, the :obj:`LSHForest`
class also supports the parallelized generation of a k-nearest neighbor graph
using the method :obj:`get_knn_graph()`.

.. code-block:: python

    # ...

    # Construct the k-nearest neighbour graph
    start = timer()
    knng_from = tm.VectorUint()
    knng_to = tm.VectorUint()
    knng_weight = tm.VectorFloat()

    result = lf.get_knn_graph(knng_from, knng_to, knng_weight, 10)
    print(f'The kNN search took {(timer() - start) * 1000}ms.')

>>> The kNN search took 37519.07863999986ms.

Layout
^^^^^^
TMAP ships with the function :obj:`layout_from_lsh_forest()` which
creates a graph / tree layout directly from an :obj:`LSHForest` instance.


The resulting layout can then be plotted using matplotlib / pyplot using
its :obj:`plot()` and :obj:`scatter` methods.

.. code-block:: python

    # ...

    # The configuration for the MST plot
    # Distribute the tree more evenly
    cfg = tm.LayoutConfiguration()
    cfg.sl_scaling_min = 1
    cfg.sl_scaling_max = 1
    cfg.node_size = 1 / 50

    # Construct the k-nearest neighbour graph
    start = timer()
    x, y, s, t, _ = tm.layout_from_lsh_forest(lf, config=cfg)
    print(f'layout_from_lsh_forest took {(timer() - start) * 1000}ms.')

    # Plot spanning tree layout
    start = timer()
    for i in range(len(s)):
        plt.plot([x[s[i]], x[t[i]]], [y[s[i]], y[t[i]]], 'r-',
                 linewidth=1.0, alpha=0.5, zorder=1)

    plt.scatter(x, y, s=0.1, zorder=2)

    plt.savefig('lsh_forest_knng_mpl.png')
    print(f'Plotting using matplotlib took {(timer() - start) * 1000}ms.')

>>> layout_from_lsh_forest took 1218.4765429992694ms.
>>> Plotting using matplotlib took 35739.334431000316ms.

.. image:: _static/lsh_forest_knng_mpl.png
   :alt: Plotting a 10,000 node spanning tree with matplotlib / pyplot.

Using matplotlib / pyplot has tow main disadvantages: It is slow and does
not yield interactive plots. For this reason, we suggest to use
the Python package `Faerun <https://pypi.org/project/faerun/>`_
for large scale data sets. Faerun supports millions of data points in web-based
visualizations.


Together with TMAP, Faerun can easily create visualizations of more than
10 million data points including associated web links and structure drawings
for high dimensional chemical data sets within an hour.

.. image:: _static/fdb.png
   :alt: Plotting millions of data points using Faerun.

.. image:: _static/fdb_zoom.png
   :alt: Faerun plots are fully interactive.
