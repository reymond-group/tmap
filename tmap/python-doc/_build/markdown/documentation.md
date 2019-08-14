# Documentation

## MinHash


#### class tmap.Minhash(self: tmap.tmap.Minhash, d: int=128, seed: int=42, sample_size: int=128)
A generator for MinHash vectors that supports binary, indexed, string and also `int` and `float` weighted vectors as input.

Constructor for the class `Minhash`.


* **Keyword Arguments**

    * **d** (`int`) – The number of permutations used for hashing

    * **seed** (`int`) – The seed used for the random number generator(s)

    * **sample_size** (`int`) – The sample size when generating a weighted MinHash



#### __init__(self: tmap.tmap.Minhash, d: int=128, seed: int=42, sample_size: int=128)
Constructor for the class `Minhash`.


* **Keyword Arguments**

    * **d** (`int`) – The number of permutations used for hashing

    * **seed** (`int`) – The seed used for the random number generator(s)

    * **sample_size** (`int`) – The sample size when generating a weighted MinHash



#### batch_from_binary_array(self: tmap.tmap.Minhash, arg0: List[tmap.tmap.VectorUchar])
Create MinHash vectors from binary arrays (parallelized).


* **Parameters**

    **vec** (`List` of `VectorUchar`) – A list of vectors containing binary values



* **Returns**

    A list of MinHash vectors



* **Return type**

    `List` of `VectorUint`



#### batch_from_int_weight_array(self: tmap.tmap.Minhash, arg0: List[tmap.tmap.VectorUint])
Create MinHash vectors from `int` arrays, where entries are weights rather than indices of ones (parallelized).


* **Parameters**

    **vec** (`List` of `VectorUint`) – A list of vectors containing `int` values



* **Returns**

    A list of MinHash vectors



* **Return type**

    `List` of `VectorUint`



#### batch_from_sparse_binary_array(self: tmap.tmap.Minhash, arg0: List[tmap.tmap.VectorUint])
Create MinHash vectors from sparse binary arrays (parallelized).


* **Parameters**

    **vec** (`List` of `VectorUint`) – A list of vectors containing indices of ones in a binary array



* **Returns**

    A list of MinHash vectors



* **Return type**

    `List` of `VectorUint`



#### batch_from_string_array(self: tmap.tmap.Minhash, arg0: List[List[str]])
Create MinHash vectors from string arrays (parallelized).


* **Parameters**

    **vec** (`List` of `List` of `str`) – A list of list of strings



* **Returns**

    A list of MinHash vectors



* **Return type**

    `List` of `VectorUint`



#### batch_from_weight_array(self: tmap.tmap.Minhash, arg0: List[tmap.tmap.VectorFloat])
Create MinHash vectors from `float` arrays (parallelized).


* **Parameters**

    **vec** (`List` of `VectorFloat`) – A list of vectors containing `float` values



* **Returns**

    A list of MinHash vectors



* **Return type**

    `List` of `VectorUint`



#### from_binary_array(self: tmap.tmap.Minhash, arg0: tmap.tmap.VectorUchar)
Create a MinHash vector from a binary array.


* **Parameters**

    **vec** (`VectorUchar`) – A vector containing binary values



* **Returns**

    A MinHash vector



* **Return type**

    `VectorUint`



#### from_sparse_binary_array(self: tmap.tmap.Minhash, arg0: tmap.tmap.VectorUint)
Create a MinHash vector from a sparse binary array.


* **Parameters**

    **vec** (`VectorUint`) – A vector containing indices of ones in a binary array



* **Returns**

    A MinHash vector



* **Return type**

    `VectorUint`



#### from_string_array(self: tmap.tmap.Minhash, arg0: List[str])
Create a MinHash vector from a string array.


* **Parameters**

    **vec** (`List` of `str`) – A vector containing strings



* **Returns**

    A MinHash vector



* **Return type**

    `VectorUint`



#### from_weight_array(self: tmap.tmap.Minhash, arg0: tmap.tmap.VectorFloat)
Create a MinHash vector from a `float` array.


* **Parameters**

    **vec** (`VectorFloat`) – A vector containing `float` values



* **Returns**

    A MinHash vector



* **Return type**

    `VectorUint`



#### get_distance(self: tmap.tmap.Minhash, arg0: tmap.tmap.VectorUint, arg1: tmap.tmap.VectorUint)
Calculate the Jaccard distance between two MinHash vectors.


* **Parameters**

    * **vec_a** (`VectorUint`) – A MinHash vector

    * **vec_b** (`VectorUint`) – A MinHash vector



* **Returns**

    `float` The Jaccard distance



#### get_weighted_distance(self: tmap.tmap.Minhash, arg0: tmap.tmap.VectorUint, arg1: tmap.tmap.VectorUint)
Calculate the weighted Jaccard distance between two MinHash vectors.


* **Parameters**

    * **vec_a** (`VectorUint`) – A weighted MinHash vector

    * **vec_b** (`VectorUint`) – A weighted MinHash vector



* **Returns**

    `float` The Jaccard distance


## LSH Forest


#### class tmap.LSHForest(self: tmap.tmap.LSHForest, d: int=128, l: int=8, store: bool=True, file_backed: bool=False)
A LSH forest data structure which incorporates optional linear scan to increase the recovery performance. Most query methods are available in parallelized versions named with a `batch_` prefix.

Constructor for the class `LSHForest`.


* **Keyword Arguments**

    * **d** (`int`) – The dimensionality of the MinHashe vectors to be added

    * **l** (`int`) – The number of prefix trees used when indexing data

    * **store** (`bool`) – 

    * **file_backed** (`bool`) Whether to store the data on disk rather than in main memory (experimental) – 



#### __init__(self: tmap.tmap.LSHForest, d: int=128, l: int=8, store: bool=True, file_backed: bool=False)
Constructor for the class `LSHForest`.


* **Keyword Arguments**

    * **d** (`int`) – The dimensionality of the MinHashe vectors to be added

    * **l** (`int`) – The number of prefix trees used when indexing data

    * **store** (`bool`) – 

    * **file_backed** (`bool`) Whether to store the data on disk rather than in main memory (experimental) – 



#### add(self: tmap.tmap.LSHForest, arg0: tmap.tmap.VectorUint)
Add a MinHash vector to the LSH forest.


* **Parameters**

    **vecs** (`VectorUint`) – A MinHash vector that is to be added to the LSH forest



#### batch_add(self: tmap.tmap.LSHForest, arg0: List[tmap.tmap.VectorUint])
Add a list MinHash vectors to the LSH forest (parallelized).


* **Parameters**

    **vecs** (`List` of `VectorUint`) – A list of MinHash vectors that is to be added to the LSH forest



#### batch_query(self: tmap.tmap.LSHForest, arg0: List[tmap.tmap.VectorUint], arg1: int)
Query the LSH forest for k-nearest neighbors (parallelized).


* **Parameters**

    * **vecs** (`List` of `VectorUint`) – The query MinHash vectors

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Returns**

    The results of the queries



* **Return type**

    `List` of `VectorUint`



#### clear(self: tmap.tmap.LSHForest)
Clears all the added data and computed indices from this `LSHForest` instance.


#### get_all_distances(self: tmap.tmap.LSHForest, arg0: tmap.tmap.VectorUint)
Calculate the Jaccard distances of a MinHash vector to all indexed MinHash vectors.


* **Parameters**

    **vec** (`VectorUint`) – The query MinHash vector



* **Returns**

    The Jaccard distances



* **Return type**

    `List` of `float`



#### get_all_nearest_neighbors(self: tmap.tmap.LSHForest, k: int, kc: int=10, weighted: bool=False)
Get the k-nearest neighbors of all indexed MinHash vectors.


* **Parameters**

    **k** (`int`) – The number of nearest neighbors to be retrieved



* **Keyword Arguments**

    * **kc** (`int`) – The factor by which `k` is multiplied for LSH forest retreival

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



* **Returns**

    `VectorUint` The ids of all k-nearest neighbors



#### get_distance(self: tmap.tmap.LSHForest, arg0: tmap.tmap.VectorUint, arg1: tmap.tmap.VectorUint)
Calculate the Jaccard distance between two MinHash vectors.


* **Parameters**

    * **vec_a** (`VectorUint`) – A MinHash vector

    * **vec_b** (`VectorUint`) – A MinHash vector



* **Returns**

    `float` The Jaccard distance



#### get_distance_by_id(self: tmap.tmap.LSHForest, arg0: int, arg1: int)
Calculate the Jaccard distance between two indexed MinHash vectors.


* **Parameters**

    * **a** (`int`) – The id of an indexed MinHash vector

    * **b** (`int`) – The id of an indexed MinHash vector



* **Returns**

    `float` The Jaccard distance



#### get_hash(self: tmap.tmap.LSHForest, arg0: int)
Retrieve the MinHash vector of an indexed entry given its index. The index is defined by order of insertion.


* **Parameters**

    **a** (`int`) – The id of an indexed MinHash vector



* **Returns**

    `VectorUint` The MinHash vector



#### get_knn_graph(self: tmap.tmap.LSHForest, from: tmap.tmap.VectorUint, to: tmap.tmap.VectorUint, weight: tmap.tmap.VectorFloat, k: int, kc: int=10, weighted: bool=False)
Construct the k-nearest neighbor graph of the indexed MinHash vectors. It will be written to out parameters `from`, `to`, and `weight` as an edge list.


* **Parameters**

    * **from** (`VectorUint`) – A vector to which the ids for the from vertices are written

    * **to** (`VectorUint`) – A vector to which the ids for the to vertices are written

    * **weight** (`VectorFloat`) – A vector to which the edge weights are written

    * **k** (`int`) – The number of nearest neighbors to be retrieved during the construction of the k-nearest neighbor graph



* **Keyword Arguments**

    * **kc** (`int`) – The factor by which `k` is multiplied for LSH forest retreival

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



#### get_weighted_distance(self: tmap.tmap.LSHForest, arg0: tmap.tmap.VectorUint, arg1: tmap.tmap.VectorUint)
Calculate the weighted Jaccard distance between two MinHash vectors.


* **Parameters**

    * **vec_a** (`VectorUint`) – A weighted MinHash vector

    * **vec_b** (`VectorUint`) – A weighted MinHash vector



* **Returns**

    `float` The Jaccard distance



#### get_weighted_distance_by_id(self: tmap.tmap.LSHForest, arg0: int, arg1: int)
Calculate the Jaccard distance between two indexed weighted MinHash vectors.


* **Parameters**

    * **a** (`int`) – The id of an indexed weighted MinHash vector

    * **b** (`int`) – The id of an indexed weighted MinHash vector



* **Returns**

    `float` The weighted Jaccard distance



#### index(self: tmap.tmap.LSHForest)
Index the LSH forest. This has to be run after each time new MinHashes were added.


#### is_clean(self: tmap.tmap.LSHForest)
Returns a boolean indicating whether or not the LSH forest has been indexed after the last MinHash vector was added.


* **Returns**

    `True` if `index()` has been run since MinHash vectors have last been added using `add()` or `batch_add()`. `False` otherwise



* **Return type**

    `bool`



#### linear_scan(self: tmap.tmap.LSHForest, vec: tmap.tmap.VectorUint, indices: tmap.tmap.VectorUint, k: int=10, weighted: bool=False)
Query a subset of indexed MinHash vectors using linear scan.


* **Parameters**

    * **vec** (`VectorUint`) – The query MinHash vector

    * **indices** (`VectorUint`) – 



* **Keyword Arguments**

    * **k** (`int`) – The number of nearest neighbors to be retrieved

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



* **Returns**

    The results of the query



* **Return type**

    `List` of `Tuple[float, int]`



#### query(self: tmap.tmap.LSHForest, arg0: tmap.tmap.VectorUint, arg1: int)
Query the LSH forest for k-nearest neighbors.


* **Parameters**

    * **vec** (`VectorUint`) – The query MinHash vector

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Returns**

    The results of the query



* **Return type**

    `VectorUint`



#### query_by_id(self: tmap.tmap.LSHForest, arg0: int, arg1: int)
Query the LSH forest for k-nearest neighbors.


* **Parameters**

    * **id** (`int`) – The id of an indexed MinHash vector

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Returns**

    The results of the query



* **Return type**

    `VectorUint`



#### query_exclude(self: tmap.tmap.LSHForest, arg0: tmap.tmap.VectorUint, arg1: tmap.tmap.VectorUint, arg2: int)
Query the LSH forest for k-nearest neighbors.


* **Parameters**

    * **vec** (`VectorUint`) – The query MinHash vector

    * **exclude** (`List` of `VectorUint`) – 

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Returns**

    The results of the query



* **Return type**

    `VectorUint`



#### query_exclude_by_id(self: tmap.tmap.LSHForest, arg0: int, arg1: tmap.tmap.VectorUint, arg2: int)
Query the LSH forest for k-nearest neighbors.


* **Parameters**

    * **id** (`int`) – The id of an indexed MinHash vector

    * **exclude** (`List` of `VectorUint`) – 

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Returns**

    The results of the query



* **Return type**

    `VectorUint`



#### query_linear_scan(self: tmap.tmap.LSHForest, vec: tmap.tmap.VectorUint, k: int, kc: int=10, weighted: bool=False)
Query k-nearest neighbors with a LSH forest / linear scan combination. `k\`\*:obj:\`kc` nearest neighbors are searched for using LSH forest; from these, the `k` nearest neighbors are retrieved using linear scan.


* **Parameters**

    * **vec** (`VectorUint`) – The query MinHash vector

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Keyword Arguments**

    * **kc** (`int`) – The factor by which `k` is multiplied for LSH forest retreival

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



* **Returns**

    The results of the query



* **Return type**

    `List` of `Tuple[float, int]`



#### query_linear_scan_by_id(self: tmap.tmap.LSHForest, id: int, k: int, kc: int=10, weighted: bool=False)
Query k-nearest neighbors with a LSH forest / linear scan combination. `k\`\*:obj:\`kc` nearest neighbors are searched for using LSH forest; from these, the `k` nearest neighbors are retrieved using linear scan.


* **Parameters**

    * **id** (`int`) – The id of an indexed MinHash vector

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Keyword Arguments**

    * **kc** (`int`) – The factor by which `k` is multiplied for LSH forest retreival

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



* **Returns**

    The results of the query



* **Return type**

    `List` of `Tuple[float, int]`



#### query_linear_scan_exclude(self: tmap.tmap.LSHForest, vec: tmap.tmap.VectorUint, k: int, exclude: tmap.tmap.VectorUint, kc: int=10, weighted: bool=False)
Query k-nearest neighbors with a LSH forest / linear scan combination. `k\`\*:obj:\`kc` nearest neighbors are searched for using LSH forest; from these, the `k` nearest neighbors are retrieved using linear scan.


* **Parameters**

    * **vec** (`VectorUint`) – The query MinHash vector

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Keyword Arguments**

    * **exclude** (`List` of `VectorUint`) – 

    * **kc** (`int`) – The factor by which `k` is multiplied for LSH forest retreival

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



* **Returns**

    The results of the query



* **Return type**

    `List` of `Tuple[float, int]`



#### query_linear_scan_exclude_by_id(self: tmap.tmap.LSHForest, id: int, k: int, exclude: tmap.tmap.VectorUint, kc: int=10, weighted: bool=False)
Query k-nearest neighbors with a LSH forest / linear scan combination. `k\`\*:obj:\`kc` nearest neighbors are searched for using LSH forest; from these, the `k` nearest neighbors are retrieved using linear scan.


* **Parameters**

    * **id** (`int`) – The id of an indexed MinHash vector

    * **k** (`int`) – The number of nearest neighbors to be retrieved



* **Keyword Arguments**

    * **exclude** (`List` of `VectorUint`) – 

    * **kc** (`int`) – The factor by which `k` is multiplied for LSH forest retreival

    * **weighted** (`bool`) – Whether the MinHash vectors in this `LSHForest` instance are weighted



* **Returns**

    The results of the query



* **Return type**

    `List` of `Tuple[float, int]`



#### restore(self: tmap.tmap.LSHForest, arg0: str)
Deserializes a previously serialized (using `store()`) state into this instance of `LSHForest` and recreates the index.


* **Parameters**

    **path** (`str`) – The path to the file which is deserialized



#### size(self: tmap.tmap.LSHForest)
Returns the number of MinHash vectors in this LSHForest instance.


* **Returns**

    The number of MinHash vectors



* **Return type**

    `int`



#### store(self: tmap.tmap.LSHForest, arg0: str)
Serializes the current state of this instance of `LSHForest` to the disk in binary format. The index is not serialized and has to be rebuilt after deserialization.


* **Parameters**

    **path** (`str`) – The path to which to searialize the file


## Layout


#### tmap.layout_from_lsh_forest()
layout_from_lsh_forest(lsh_forest: tmap::LSHForest, config: tmap.tmap.LayoutConfiguration=k: 10
kc: 10
fme_iterations: 1000
fme_randomize: 0
fme_threads: 4
fme_precision: 4
sl_repeats: 1
sl_extra_scaling_steps: 1
sl_scaling_x: 5.000000
sl_scaling_y: 20.000000
sl_scaling_type: RelativeToDrawing
mmm_repeats: 1
placer: Barycenter
merger: LocalBiconnected
merger_factor: 2.000000
merger_adjustment: 0
node_size1.000000, create_mst: bool=True, clear_lsh_forest: bool=False, weighted: bool=False) -> Tuple[tmap.tmap.VectorFloat, tmap.tmap.VectorFloat, tmap.tmap.VectorUint, tmap.tmap.VectorUint, tmap.tmap.GraphProperties]

> Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an `LSHForest` instance.

> Arguments:

>     lsh_forest (`LSHForest`): An `LSHForest` instance

> Keyword Arguments:

>     config (`LayoutConfiguration`, optional): An `LayoutConfiguration` instance
>     create_mst (`bool`): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph
>     clear_lsh_forest (`bool`): Whether to run `clear()` on the `LSHForest` instance after k-nearest negihbor graph and MST creation and before layout
>     weighted (`bool`): Whether the MinHash vectors in the `LSHForest` instance are weighted

> Returns:

>     `Tuple[VectorFloat, VectorFloat, VectorUint, VectorUint, GraphProperties]` The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph


#### tmap.layout_from_edge_list()
layout_from_edge_list(vertex_count: int, edges: List[Tuple[int, int, float]], config: tmap.tmap.LayoutConfiguration=k: 10
kc: 10
fme_iterations: 1000
fme_randomize: 0
fme_threads: 4
fme_precision: 4
sl_repeats: 1
sl_extra_scaling_steps: 1
sl_scaling_x: 5.000000
sl_scaling_y: 20.000000
sl_scaling_type: RelativeToDrawing
mmm_repeats: 1
placer: Barycenter
merger: LocalBiconnected
merger_factor: 2.000000
merger_adjustment: 0
node_size1.000000, create_mst: bool=True) -> Tuple[tmap.tmap.VectorFloat, tmap.tmap.VectorFloat, tmap.tmap.VectorUint, tmap.tmap.VectorUint, tmap.tmap.GraphProperties]

> Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an edge list.

> Arguments:

>     vertex_count (`int`): The number of vertices in the edge list
>     edges (`List` of `Tuple[int, int, float]`): An edge list defining a graph

> Keyword Arguments:

>     config (`LayoutConfiguration`, optional): An `LayoutConfiguration` instance
>     create_mst (`bool`): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph

> Returns:

>     `Tuple[VectorFloat, VectorFloat, VectorUint, VectorUint, GraphProperties]`: The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph


#### tmap.mst_from_lsh_forest(lsh_forest: tmap::LSHForest, k: int, kc: int=10, weighted: bool=False)
Create minimum spanning tree topology from an `LSHForest` instance.


* **Parameters**

    * **lsh_forest** (`LSHForest`) – An `LSHForest` instance

    * **k** (*int*) – The number of nearest neighbors used to create the k-nearest neighbor graph



* **Keyword Arguments**

    * **kc** (*int*) – The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned

    * **weighted** (`bool`) – Whether the MinHash vectors in the `LSHForest` instance are weighted



* **Returns**

    the topology of the minimum spanning tree of the data indexed in the LSH forest



* **Return type**

    `Tuple[VectorUint, VectorUint]`



#### class tmap.ScalingType(self: tmap.tmap.ScalingType, arg0: int)
The scaling types available in OGDF. The class is to be used as an enum.

### Notes

The available values are

`ScalingType.Absolute`: Absolute factor, can be used to scale relative to level size change.

`ScalingType.RelativeToAvgLength`: Scales by a factor relative to the average edge weights.

`ScalingType.RelativeToDesiredLength`: Scales by a factor relative to the disired edge length.

`ScalingType.RelativeToDrawing`: Scales by a factor relative to the drawing.


#### class tmap.Placer(self: tmap.tmap.Placer, arg0: int)
The places available in OGDF. The class is to be used as an enum.

### Notes

The available values are

`Placer.Barycenter`: Places a vertex at the barycenter of its neighbors’ position.

`Placer.Solar`: Uses information of the merging phase of the solar merger. Places a new vertex on the direct line between two suns.

`Placer.Circle`: Places the vertices in a circle around the barycenter and outside of the current drawing

`Placer.Median`: Places a vertex at the median position of the neighbor nodes for each coordinate axis.

`Placer.Random`: Places a vertex at a random position within the smallest circle containing all vertices around the barycenter of the current drawing.

`Placer.Zero`: Places a vertex at the same position as its representative in the previous level.


#### class tmap.Merger(self: tmap.tmap.Merger, arg0: int)
The mergers available in OGDF. The class is to be used as an enum.

### Notes

The available values are

`Merger.EdgeCover`: Based on the matching merger. Computes an edge cover such that each contained edge is incident to at least one unmatched vertex. The cover edges are then used to merge their adjacent vertices.

`Merger.LocalBiconnected`: Based on the edge cover merger. Avoids distortions by checking whether biconnectivity will be lost in the local neighborhood around the potential merging position.

`Merger.Solar`: Vertices are partitioned into solar systems, consisting of sun, planets and moons. The systems are then merged into the sun vertices.

`Merger.IndependentSet`: Uses a maximal independent set filtration. See GRIP for details.


#### class tmap.LayoutConfiguration(self: tmap.tmap.LayoutConfiguration)
A container for configuration options for `layout_from_lsh_forest()` and `layout_from_edge_list()`.


#### int k()
The number of nearest neighbors used to create the k-nearest neighbor graph.


* **Type**

    `int`



#### int kc()
The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned.


* **Type**

    `int`



#### int fme_iterations()
Maximum number of iterations of the fast multipole embedder.


* **Type**

    `int`



#### bool fme_randomize()
Whether or not to randomize the layout at the start.


* **Type**

    `bool`



#### int fme_threads()
The number of threads for the fast multipole embedder.


* **Type**

    `int`



#### int fme_precision()
The number of coefficients of the multipole expansion.


* **Type**

    `int`



#### int sl_repeats()
The number of repeats of the scaling layout algorithm.


* **Type**

    `int`



#### int sl_extra_scaling_steps()
Sets the number of repeats of the scaling.


* **Type**

    `int`



#### double sl_scaling_min()
The minimum scaling factor.


* **Type**

    `float`



#### double sl_scaling_max()
The maximum scaling factor.


* **Type**

    `float`



#### ScalingType sl_scaling_type()
Defines the (relative) scale of the graph.


* **Type**

    `ScalingType`



#### int mmm_repeats()
Number of repeats of the per-level layout algorithm.


* **Type**

    `int`



#### Placer placer()
The  method  by  which  the  initial  positons  of  the  vertices  at  eachlevel are defined.


* **Type**

    `Placer`



#### Merger merger()
The vertex merging strategy applied during the coarsening phaseof the multilevel algorithm.


* **Type**

    `Merger`



#### double merger_factor()
The ratio of the sizes between two levels up to which the mergingis run.  Does not apply to all merging strategies.


* **Type**

    `float`



#### int merger_adjustment()
The  edge  length  adjustment  of  the  merging  algorithm.   Does  notapply to all merging strategies.


* **Type**

    `int`



#### float node_size()
The size of the nodes, which affects the magnitude of their repellingforce. Decreasing  this  value  generally  resolves  overlaps  in  a  verycrowded tree.


* **Type**

    `float`


Constructor for the class `LayoutConfiguration`.


#### __init__(self: tmap.tmap.LayoutConfiguration)
Constructor for the class `LayoutConfiguration`.


#### class tmap.GraphProperties(self: tmap.tmap.GraphProperties)
Contains properties of the minimum spanning tree (or forest) generated by `layout_from_lsh_forest()` and `layout_from_edge_list()`.


#### mst_weight()
The total weight of the minimum spanning tree.


* **Type**

    `float`



#### n_connected_components()
The number of connected components in the minimum spanning forest.


* **Type**

    `int`



#### n_isolated_vertices()

* **Type**

    `int`



#### degrees()
The degrees of all vertices in the minimum spanning tree (or forest).


* **Type**

    `VectorUint`



#### adjacency_list()
The adjaceny lists for all vertices in the minimum spanning tree (or forest).


* **Type**

    `List` of `VectorUint`


Constructor for the class `GraphProperties`.


#### __init__(self: tmap.tmap.GraphProperties)
Constructor for the class `GraphProperties`.
