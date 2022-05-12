# Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`namespace `[`tmap`](#namespacetmap) | 

# namespace `tmap` <a id="namespacetmap"></a>

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`enum `[`Placer`](#layout_8hh_1afdc98947e81dc6f4c30f256e6f42f90b)            | The placers available in OGDF.
`enum `[`Merger`](#layout_8hh_1a8c7bb9956a1a724233182a166cfdc0ff)            | The mergers available in OGDF.
`enum `[`ScalingType`](#layout_8hh_1a50ec215c9e54cf12b9dd0a0056160761)            | The scaling types available in OGDF.
`public std::tuple< std::vector< float >, std::vector< float >, std::vector< uint32_t >, std::vector< uint32_t >, `[`GraphProperties`](#structtmap_1_1GraphProperties)` > `[`LayoutFromLSHForest`](#layout_8hh_1acd9c409403d706202320359541674ba8)`(`[`LSHForest`](#classtmap_1_1LSHForest)` & lsh_forest,`[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration)` config,bool create_mst,bool clear_lsh_forest,bool weighted)`            | Genereates coordinates, edges and properties of a MST (via a kNN graph) from an [LSHForest](#classtmap_1_1LSHForest) instance.
`public std::tuple< std::vector< uint32_t >, std::vector< uint32_t > > `[`MSTFromLSHForest`](#layout_8hh_1a3d7180c2fcc2a1ec72838f18177d6840)`(`[`LSHForest`](#classtmap_1_1LSHForest)` & lsh_forest,uint32_t k,uint32_t kc,bool weighted)`            | Generates an MST (via a kNN graph) from an [LSHForest](#classtmap_1_1LSHForest) instance.
`public std::tuple< std::vector< float >, std::vector< float >, std::vector< uint32_t >, std::vector< uint32_t >, `[`GraphProperties`](#structtmap_1_1GraphProperties)` > `[`LayoutFromEdgeList`](#layout_8hh_1a780993ad8dd7e349b77f55895cc33451)`(uint32_t vertex_count,const std::vector< std::tuple< uint32_t, uint32_t, float >> & edges,`[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration)` config,bool create_mst)`            | Genereates coordinates, edges and properties of a MST from an edge list.
`public std::tuple< std::vector< float >, std::vector< float >, std::vector< uint32_t >, std::vector< uint32_t >, `[`GraphProperties`](#structtmap_1_1GraphProperties)` > `[`LayoutInternal`](#layout_8hh_1a126dbc6ec8355732c528abb2877e60d4)`(ogdf::EdgeWeightedGraph< float > & g,uint32_t vertex_count,`[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration)` config,`[`GraphProperties`](#structtmap_1_1GraphProperties)` & gp)`            | Laying out an OGDF graph.
`class `[`tmap::LSHForest`](#classtmap_1_1LSHForest) | Provides locality sensitive hashing forest functionalities.
`class `[`tmap::Minhash`](#classtmap_1_1Minhash) | An implementation of MinHash and weighted MinHash using SHA1.
`class `[`tmap::Timer`](#classtmap_1_1Timer) | A simple timer class used to check performance during development.
`struct `[`tmap::GraphProperties`](#structtmap_1_1GraphProperties) | The properties of a generated graph. An instance of this struct is returned from the layout functions.
`struct `[`tmap::LayoutConfiguration`](#structtmap_1_1LayoutConfiguration) | A struct containing all the configuration options available for and applied to a layout.
`struct `[`tmap::MapKeyPointer`](#structtmap_1_1MapKeyPointer) | The pointer map used for pointing to the keys from the sorted hash map.
`struct `[`tmap::SimpleHash`](#structtmap_1_1SimpleHash) | Hash struct used for the sparsepp sparse hash map.

## Members

#### `enum `[`Placer`](#layout_8hh_1afdc98947e81dc6f4c30f256e6f42f90b) <a id="layout_8hh_1afdc98947e81dc6f4c30f256e6f42f90b"></a>

 Values                         | Descriptions                                
--------------------------------|---------------------------------------------
Barycenter            | Places a vertex at the barycenter of its neighbors' position.
Solar            | Uses information of the merging phase of the solar merger. Places a new vertex on the direct line between two suns.
Circle            | Places the vertices in a circle around the barycenter and outside of the current drawing.
Median            | Places a vertex at the median position of the neighbor nodes for each coordinate axis.
Random            | Places a vertex at a random position within the smallest circle containing all vertices around the barycenter of the current drawing.
Zero            | Places a vertex at the same position as its representative in the previous level.

The placers available in OGDF.

#### `enum `[`Merger`](#layout_8hh_1a8c7bb9956a1a724233182a166cfdc0ff) <a id="layout_8hh_1a8c7bb9956a1a724233182a166cfdc0ff"></a>

 Values                         | Descriptions                                
--------------------------------|---------------------------------------------
EdgeCover            | Based on the matching merger. Computes an edge cover such that each contained edge is incident to at least one unmatched vertex. The cover edges are then used to merge their adjacent vertices.
LocalBiconnected            | Based on the edge cover merger. Avoids distortions by checking whether biconnectivity will be lost in the local neighborhood around the potential merging position.
Solar            | Vertices are partitioned into solar systems, consisting of sun, planets and moons. The systems are then merged into the sun vertices.
IndependentSet            | Uses a maximal independent set filtration. See GRIP for details.

The mergers available in OGDF.

#### `enum `[`ScalingType`](#layout_8hh_1a50ec215c9e54cf12b9dd0a0056160761) <a id="layout_8hh_1a50ec215c9e54cf12b9dd0a0056160761"></a>

 Values                         | Descriptions                                
--------------------------------|---------------------------------------------
Absolute            | Absolute factor, can be used to scale relative to level size change.
RelativeToAvgLength            | Scales by a factor relative to the average edge weights.
RelativeToDesiredLength            | Scales by a factor relative to the disired edge length.
RelativeToDrawing            | Scales by a factor relative to the drawing.

The scaling types available in OGDF.

#### `public std::tuple< std::vector< float >, std::vector< float >, std::vector< uint32_t >, std::vector< uint32_t >, `[`GraphProperties`](#structtmap_1_1GraphProperties)` > `[`LayoutFromLSHForest`](#layout_8hh_1acd9c409403d706202320359541674ba8)`(`[`LSHForest`](#classtmap_1_1LSHForest)` & lsh_forest,`[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration)` config,bool create_mst,bool clear_lsh_forest,bool weighted)` <a id="layout_8hh_1acd9c409403d706202320359541674ba8"></a>

Genereates coordinates, edges and properties of a MST (via a kNN graph) from an [LSHForest](#classtmap_1_1LSHForest) instance.

#### Parameters
* `lsh_forest` An [LSHForest](#classtmap_1_1LSHForest) instance which is used to construct the kNN graph. 

* `config` A [LayoutConfiguration](#structtmap_1_1LayoutConfiguration) instance. 

* `create_mst` Whether to create an MST before laying out the graph. 

* `clear_lsh_forest` Whether to clear the [LSHForest](#classtmap_1_1LSHForest) after it's use (might save memory). 

* `weighted` Whether the [LSHForest](#classtmap_1_1LSHForest) instance contains weighted MinHash data. 

#### Returns
std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, [GraphProperties](#structtmap_1_1GraphProperties)>

#### `public std::tuple< std::vector< uint32_t >, std::vector< uint32_t > > `[`MSTFromLSHForest`](#layout_8hh_1a3d7180c2fcc2a1ec72838f18177d6840)`(`[`LSHForest`](#classtmap_1_1LSHForest)` & lsh_forest,uint32_t k,uint32_t kc,bool weighted)` <a id="layout_8hh_1a3d7180c2fcc2a1ec72838f18177d6840"></a>

Generates an MST (via a kNN graph) from an [LSHForest](#classtmap_1_1LSHForest) instance.

#### Parameters
* `lsh_forest` An [LSHForest](#classtmap_1_1LSHForest) instance which is used to construct the kNN graph. 

* `k` The number of nearest neighbors used to create the kNN graph. 

* `kc` The factor by which k is multiplied when retrieving nearest neighbors. 

* `weighted` Whether the [LSHForest](#classtmap_1_1LSHForest) instance contains weighted MinHash data. 

#### Returns
std::tuple<std::vector<uint32_t>, std::vector<uint32_t>>

#### `public std::tuple< std::vector< float >, std::vector< float >, std::vector< uint32_t >, std::vector< uint32_t >, `[`GraphProperties`](#structtmap_1_1GraphProperties)` > `[`LayoutFromEdgeList`](#layout_8hh_1a780993ad8dd7e349b77f55895cc33451)`(uint32_t vertex_count,const std::vector< std::tuple< uint32_t, uint32_t, float >> & edges,`[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration)` config,bool create_mst)` <a id="layout_8hh_1a780993ad8dd7e349b77f55895cc33451"></a>

Genereates coordinates, edges and properties of a MST from an edge list.

#### Parameters
* `vertex_count` The number of vertices in the input graph. 

* `edges` An edge list in the form of [(from, to, weight)]. 

* `config` A [LayoutConfiguration](#structtmap_1_1LayoutConfiguration) instance. 

* `create_mst` Whether to create an MST before laying out the graph. 

#### Returns
std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, [GraphProperties](#structtmap_1_1GraphProperties)>

#### `public std::tuple< std::vector< float >, std::vector< float >, std::vector< uint32_t >, std::vector< uint32_t >, `[`GraphProperties`](#structtmap_1_1GraphProperties)` > `[`LayoutInternal`](#layout_8hh_1a126dbc6ec8355732c528abb2877e60d4)`(ogdf::EdgeWeightedGraph< float > & g,uint32_t vertex_count,`[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration)` config,`[`GraphProperties`](#structtmap_1_1GraphProperties)` & gp)` <a id="layout_8hh_1a126dbc6ec8355732c528abb2877e60d4"></a>

Laying out an OGDF graph.

#### Parameters
* `g` An OGDF Graph instance 

* `vertex_count` The number of vertices in the graph. 

* `config` A [LayoutConfiguration](#structtmap_1_1LayoutConfiguration) instance. 

* `gp` An instance of a [GraphProperties](#structtmap_1_1GraphProperties) struct. 

#### Returns
std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, [GraphProperties](#structtmap_1_1GraphProperties)>

# class `tmap::LSHForest` <a id="classtmap_1_1LSHForest"></a>

Provides locality sensitive hashing forest functionalities.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public  `[`LSHForest`](#classtmap_1_1LSHForest_1a153cb1f5090432257a17f2e9dacc32a0)`(unsigned int d,unsigned int l,bool store,bool file_backed)` | Construct a new [LSHForest](#classtmap_1_1LSHForest) object.
`public inline  `[`~LSHForest`](#classtmap_1_1LSHForest_1a3ab5789f702f9dac3f801c7b9d53afd4)`()` | Destroy the [LSHForest](#classtmap_1_1LSHForest) object.
`public void `[`Add`](#classtmap_1_1LSHForest_1a480d0de16bc1e4b1365bf97b9b60223a)`(std::vector< uint32_t > & vec)` | Add a MinHash to this [LSHForest](#classtmap_1_1LSHForest) instance.
`public void `[`BatchAdd`](#classtmap_1_1LSHForest_1ab3f73f59918a37b63662679461828cbb)`(std::vector< std::vector< uint32_t >> & vecs)` | Add Minhashes to this [LSHForest](#classtmap_1_1LSHForest) (parallelized).
`public void `[`Index`](#classtmap_1_1LSHForest_1aba68c9cab8cc3c32e684e08b4f9d0a33)`()` | Create the index (trees).
`public bool `[`IsClean`](#classtmap_1_1LSHForest_1a7785c1a7f17eddd5e943db4b5d6d7cf2)`()` | Check whether the added MinHashes have been indexed.
`public void `[`Store`](#classtmap_1_1LSHForest_1a1731bf94cd09e7ebc4a10dd42145dc51)`(const std::string & path)` | Write / serialize the current LSH forest to the disk.
`public void `[`Restore`](#classtmap_1_1LSHForest_1a869273bd3d4c72c4c296dc42519558c8)`(const std::string & path)` | Read / deserialize a LSH forest instance form the disk. The forest is indexed automatically.
`public std::vector< uint32_t > `[`GetHash`](#classtmap_1_1LSHForest_1a78106dd9e3a9a9ec012e5405445be78c)`(uint32_t id)` | Get the MinHash of an entry at a given index. The index is defined by order of insertion.
`public void `[`GetKNNGraph`](#classtmap_1_1LSHForest_1a11ccbeea4356cce12b579566925d865f)`(std::vector< uint32_t > & from,std::vector< uint32_t > & to,std::vector< float > & weight,unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbor graph of the data stored in this LSH forest instance. It will be written to out parameters as an edge list.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScan`](#classtmap_1_1LSHForest_1a2ba770074cd9c0e6860b30679793c569)`(const std::vector< uint32_t > & vec,unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of a query.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExclude`](#classtmap_1_1LSHForest_1a5afd77e1f9349edcdec32a1d7aa3f38e)`(const std::vector< uint32_t > & vec,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of a query except those defined in the argument exclude.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanById`](#classtmap_1_1LSHForest_1ae4e013129270d53af27091c2f3e4e5d6)`(uint32_t id,unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of an entry.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExcludeById`](#classtmap_1_1LSHForest_1af8cff8cd9cf3b1d30823e9d36938745b)`(uint32_t id,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of an entry except those defined in the argument exclude.
`public std::vector< std::pair< float, uint32_t > > `[`LinearScan`](#classtmap_1_1LSHForest_1a63e22972fbd38e3edc8f321f88a6be8d)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & indices,unsigned int k,bool weighted)` | Get the k-nearest neighbors of a query using linear scan.
`public std::vector< uint32_t > `[`Query`](#classtmap_1_1LSHForest_1a0da6325b50a92db6ff6c49bd62a5e95b)`(const std::vector< uint32_t > & vec,unsigned int k)` | Query the LSH forest for k-nearest neighbors.
`public std::vector< uint32_t > `[`QueryExclude`](#classtmap_1_1LSHForest_1a7aba9b1df0273b71ed2d9233d05ed0db)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & exclude,unsigned int k)` | Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.
`public std::vector< uint32_t > `[`QueryById`](#classtmap_1_1LSHForest_1aa200b72cc60947e5e03fd72ea726a999)`(uint32_t id,unsigned int k)` | Query the LSH forest for k-nearest neighbors.
`public std::vector< uint32_t > `[`QueryExcludeById`](#classtmap_1_1LSHForest_1a0bcdb607c4e08e0e620b4d1d1dd12f86)`(uint32_t id,std::vector< uint32_t > & exclude,unsigned int k)` | Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.
`public std::vector< std::vector< uint32_t > > `[`BatchQuery`](#classtmap_1_1LSHForest_1adec697793677c79683490b776ae8642c)`(const std::vector< std::vector< uint32_t >> & vecs,unsigned int k)` | Query the LSH forest for k-nearest neighbors (parallelized).
`public std::vector< uint32_t > `[`GetAllNearestNeighbors`](#classtmap_1_1LSHForest_1a378f0494bce3354bb0d618558f316c84)`(unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of all LSH forest entries.
`public std::vector< uint32_t > `[`GetData`](#classtmap_1_1LSHForest_1ac4ec080057307f69548e6ca756ce5609)`(uint32_t id)` | Get the MinHash of an entry at a given index. The index is defined by order of insertion. Alias for GetHash.
`public std::vector< float > `[`GetAllDistances`](#classtmap_1_1LSHForest_1a438a46f67fb257ae85c3dd16e8b194df)`(const std::vector< uint32_t > & vec)` | Get the distances of a MinHash to all entries in the LSH forest.
`public float `[`GetDistance`](#classtmap_1_1LSHForest_1ab1c5e002deea04a625ab141f280bab92)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` | Get the distance between two MinHashes.
`public float `[`GetWeightedDistance`](#classtmap_1_1LSHForest_1aa6c035b27040909b3d7a8782ad1c63b8)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` | Get the distance between two weighted MinHashes.
`public float `[`GetDistanceById`](#classtmap_1_1LSHForest_1a8fc81622125b40114951a61cbe90863f)`(uint32_t a,uint32_t b)` | Get the distance between two MinHashes.
`public float `[`GetWeightedDistanceById`](#classtmap_1_1LSHForest_1ab00052289bb6bea152e6024049eebcc5)`(uint32_t a,uint32_t b)` | Get the distance between two weighted MinHashes.
`public void `[`Clear`](#classtmap_1_1LSHForest_1a9ee2595fb0f85d917989234ab4aaee8d)`()` | Remove all entries and the index from the LSH forest.
`public size_t `[`size`](#classtmap_1_1LSHForest_1a8ba5c1f500e915c6717c64ac24744874)`()` | Get the number of entries.

## Members

#### `public  `[`LSHForest`](#classtmap_1_1LSHForest_1a153cb1f5090432257a17f2e9dacc32a0)`(unsigned int d,unsigned int l,bool store,bool file_backed)` <a id="classtmap_1_1LSHForest_1a153cb1f5090432257a17f2e9dacc32a0"></a>

Construct a new [LSHForest](#classtmap_1_1LSHForest) object.

#### Parameters
* `d` The dimensionality of the MinHashes to be added to this [LSHForest](#classtmap_1_1LSHForest). 

* `l` The number of prefix trees used. 

* `store` Whether or not to store the data for later enhanced (using parameter kc) retrievel. 

* `file_backed` Whether to store the data on disk rather than in RAM (experimental).

#### `public inline  `[`~LSHForest`](#classtmap_1_1LSHForest_1a3ab5789f702f9dac3f801c7b9d53afd4)`()` <a id="classtmap_1_1LSHForest_1a3ab5789f702f9dac3f801c7b9d53afd4"></a>

Destroy the [LSHForest](#classtmap_1_1LSHForest) object.

#### `public void `[`Add`](#classtmap_1_1LSHForest_1a480d0de16bc1e4b1365bf97b9b60223a)`(std::vector< uint32_t > & vec)` <a id="classtmap_1_1LSHForest_1a480d0de16bc1e4b1365bf97b9b60223a"></a>

Add a MinHash to this [LSHForest](#classtmap_1_1LSHForest) instance.

#### Parameters
* `vec` A MinHash vector.

#### `public void `[`BatchAdd`](#classtmap_1_1LSHForest_1ab3f73f59918a37b63662679461828cbb)`(std::vector< std::vector< uint32_t >> & vecs)` <a id="classtmap_1_1LSHForest_1ab3f73f59918a37b63662679461828cbb"></a>

Add Minhashes to this [LSHForest](#classtmap_1_1LSHForest) (parallelized).

#### Parameters
* `vecs` A vector containing MinHash vectors.

#### `public void `[`Index`](#classtmap_1_1LSHForest_1aba68c9cab8cc3c32e684e08b4f9d0a33)`()` <a id="classtmap_1_1LSHForest_1aba68c9cab8cc3c32e684e08b4f9d0a33"></a>

Create the index (trees).

#### `public bool `[`IsClean`](#classtmap_1_1LSHForest_1a7785c1a7f17eddd5e943db4b5d6d7cf2)`()` <a id="classtmap_1_1LSHForest_1a7785c1a7f17eddd5e943db4b5d6d7cf2"></a>

Check whether the added MinHashes have been indexed.

#### Returns
true 

#### Returns
false

#### `public void `[`Store`](#classtmap_1_1LSHForest_1a1731bf94cd09e7ebc4a10dd42145dc51)`(const std::string & path)` <a id="classtmap_1_1LSHForest_1a1731bf94cd09e7ebc4a10dd42145dc51"></a>

Write / serialize the current LSH forest to the disk.

#### Parameters
* `path` The location where the LSH forest should be stored on disk.

#### `public void `[`Restore`](#classtmap_1_1LSHForest_1a869273bd3d4c72c4c296dc42519558c8)`(const std::string & path)` <a id="classtmap_1_1LSHForest_1a869273bd3d4c72c4c296dc42519558c8"></a>

Read / deserialize a LSH forest instance form the disk. The forest is indexed automatically.

#### Parameters
* `path` The location from where to load the LSH forest.

#### `public std::vector< uint32_t > `[`GetHash`](#classtmap_1_1LSHForest_1a78106dd9e3a9a9ec012e5405445be78c)`(uint32_t id)` <a id="classtmap_1_1LSHForest_1a78106dd9e3a9a9ec012e5405445be78c"></a>

Get the MinHash of an entry at a given index. The index is defined by order of insertion.

#### Parameters
* `id` The index (order of insertion) of a entry. 

#### Returns
std::vector<uint32_t> The MinHash associated with an index.

#### `public void `[`GetKNNGraph`](#classtmap_1_1LSHForest_1a11ccbeea4356cce12b579566925d865f)`(std::vector< uint32_t > & from,std::vector< uint32_t > & to,std::vector< float > & weight,unsigned int k,unsigned int kc,bool weighted)` <a id="classtmap_1_1LSHForest_1a11ccbeea4356cce12b579566925d865f"></a>

Get the k-nearest neighbor graph of the data stored in this LSH forest instance. It will be written to out parameters as an edge list.

#### Parameters
* `from` A vector to which the from vertices will be written. 

* `to` A vector to which the to vertices will be written. 

* `weight` A vector to which the float weights of the edges will be written. 

* `k` The degree of the nearest neighbor graph. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results are picked to create the k-nearest neighbor graph. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted.

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScan`](#classtmap_1_1LSHForest_1a2ba770074cd9c0e6860b30679793c569)`(const std::vector< uint32_t > & vec,unsigned int k,unsigned int kc,bool weighted)` <a id="classtmap_1_1LSHForest_1a2ba770074cd9c0e6860b30679793c569"></a>

Get the k-nearest neighbors of a query.

#### Parameters
* `vec` The query MinHash. 

* `k` The number of k-nearest neighbors to return. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExclude`](#classtmap_1_1LSHForest_1a5afd77e1f9349edcdec32a1d7aa3f38e)`(const std::vector< uint32_t > & vec,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` <a id="classtmap_1_1LSHForest_1a5afd77e1f9349edcdec32a1d7aa3f38e"></a>

Get the k-nearest neighbors of a query except those defined in the argument exclude.

#### Parameters
* `vec` The query MinHash. 

* `k` The number of k-nearest neighbors to return. 

* `exclude` A list of indices to be excluded from the search 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>>

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanById`](#classtmap_1_1LSHForest_1ae4e013129270d53af27091c2f3e4e5d6)`(uint32_t id,unsigned int k,unsigned int kc,bool weighted)` <a id="classtmap_1_1LSHForest_1ae4e013129270d53af27091c2f3e4e5d6"></a>

Get the k-nearest neighbors of an entry.

#### Parameters
* `id` The id of the query entry. 

* `k` The number of k-nearest neighbors to return. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExcludeById`](#classtmap_1_1LSHForest_1af8cff8cd9cf3b1d30823e9d36938745b)`(uint32_t id,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` <a id="classtmap_1_1LSHForest_1af8cff8cd9cf3b1d30823e9d36938745b"></a>

Get the k-nearest neighbors of an entry except those defined in the argument exclude.

#### Parameters
* `id` The id of the query entry. 

* `k` The number of k-nearest neighbors to return. 

* `exclude` A list of indices to be excluded from the search. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< std::pair< float, uint32_t > > `[`LinearScan`](#classtmap_1_1LSHForest_1a63e22972fbd38e3edc8f321f88a6be8d)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & indices,unsigned int k,bool weighted)` <a id="classtmap_1_1LSHForest_1a63e22972fbd38e3edc8f321f88a6be8d"></a>

Get the k-nearest neighbors of a query using linear scan.

#### Parameters
* `vec` The query MinHash. 

* `indices` A list of indices to on which to run the linear scan. 

* `k` The number of k-nearest neighbors to return. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`Query`](#classtmap_1_1LSHForest_1a0da6325b50a92db6ff6c49bd62a5e95b)`(const std::vector< uint32_t > & vec,unsigned int k)` <a id="classtmap_1_1LSHForest_1a0da6325b50a92db6ff6c49bd62a5e95b"></a>

Query the LSH forest for k-nearest neighbors.

#### Parameters
* `vec` The query MinHash. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`QueryExclude`](#classtmap_1_1LSHForest_1a7aba9b1df0273b71ed2d9233d05ed0db)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & exclude,unsigned int k)` <a id="classtmap_1_1LSHForest_1a7aba9b1df0273b71ed2d9233d05ed0db"></a>

Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.

#### Parameters
* `vec` The query MinHash. 

* `exclude` A list of indices to be excluded from the search. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`QueryById`](#classtmap_1_1LSHForest_1aa200b72cc60947e5e03fd72ea726a999)`(uint32_t id,unsigned int k)` <a id="classtmap_1_1LSHForest_1aa200b72cc60947e5e03fd72ea726a999"></a>

Query the LSH forest for k-nearest neighbors.

#### Parameters
* `id` The id of the query entry. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`QueryExcludeById`](#classtmap_1_1LSHForest_1a0bcdb607c4e08e0e620b4d1d1dd12f86)`(uint32_t id,std::vector< uint32_t > & exclude,unsigned int k)` <a id="classtmap_1_1LSHForest_1a0bcdb607c4e08e0e620b4d1d1dd12f86"></a>

Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.

#### Parameters
* `id` The id of the query entry. 

* `exclude` A list of indices to be excluded from the search. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< std::vector< uint32_t > > `[`BatchQuery`](#classtmap_1_1LSHForest_1adec697793677c79683490b776ae8642c)`(const std::vector< std::vector< uint32_t >> & vecs,unsigned int k)` <a id="classtmap_1_1LSHForest_1adec697793677c79683490b776ae8642c"></a>

Query the LSH forest for k-nearest neighbors (parallelized).

#### Parameters
* `vecs` A vector of MinHashes. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<std::vector<uint32_t>> A vector of the indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`GetAllNearestNeighbors`](#classtmap_1_1LSHForest_1a378f0494bce3354bb0d618558f316c84)`(unsigned int k,unsigned int kc,bool weighted)` <a id="classtmap_1_1LSHForest_1a378f0494bce3354bb0d618558f316c84"></a>

Get the k-nearest neighbors of all LSH forest entries.

#### Parameters
* `k` The number of nearest neighbors to search for. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<uint32_t> The IDs of the nearest neighbors of all LSH forest entries.

#### `public std::vector< uint32_t > `[`GetData`](#classtmap_1_1LSHForest_1ac4ec080057307f69548e6ca756ce5609)`(uint32_t id)` <a id="classtmap_1_1LSHForest_1ac4ec080057307f69548e6ca756ce5609"></a>

Get the MinHash of an entry at a given index. The index is defined by order of insertion. Alias for GetHash.

#### Parameters
* `id` The index (order of insertion) of a entry. 

#### Returns
std::vector<uint32_t> The MinHash associated with an index.

#### `public std::vector< float > `[`GetAllDistances`](#classtmap_1_1LSHForest_1a438a46f67fb257ae85c3dd16e8b194df)`(const std::vector< uint32_t > & vec)` <a id="classtmap_1_1LSHForest_1a438a46f67fb257ae85c3dd16e8b194df"></a>

Get the distances of a MinHash to all entries in the LSH forest.

#### Parameters
* `vec` The query MinHash. 

#### Returns
std::vector<float> The distances form the input MinHash to all the entries in the LSH forest.

#### `public float `[`GetDistance`](#classtmap_1_1LSHForest_1ab1c5e002deea04a625ab141f280bab92)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` <a id="classtmap_1_1LSHForest_1ab1c5e002deea04a625ab141f280bab92"></a>

Get the distance between two MinHashes.

#### Parameters
* `vec_a` A MinHash. 

* `vec_b` A MinHash. 

#### Returns
float

#### `public float `[`GetWeightedDistance`](#classtmap_1_1LSHForest_1aa6c035b27040909b3d7a8782ad1c63b8)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` <a id="classtmap_1_1LSHForest_1aa6c035b27040909b3d7a8782ad1c63b8"></a>

Get the distance between two weighted MinHashes.

#### Parameters
* `vec_a` A weighted MinHash. 

* `vec_b` A weighted MinHash. 

#### Returns
float

#### `public float `[`GetDistanceById`](#classtmap_1_1LSHForest_1a8fc81622125b40114951a61cbe90863f)`(uint32_t a,uint32_t b)` <a id="classtmap_1_1LSHForest_1a8fc81622125b40114951a61cbe90863f"></a>

Get the distance between two MinHashes.

#### Parameters
* `a` The id of an LSH forest entry. 

* `b` The id of an LSH forest entry. 

#### Returns
float

#### `public float `[`GetWeightedDistanceById`](#classtmap_1_1LSHForest_1ab00052289bb6bea152e6024049eebcc5)`(uint32_t a,uint32_t b)` <a id="classtmap_1_1LSHForest_1ab00052289bb6bea152e6024049eebcc5"></a>

Get the distance between two weighted MinHashes.

#### Parameters
* `a` The id of an LSH forest entry. 

* `b` The id of an LSH forest entry. 

#### Returns
float

#### `public void `[`Clear`](#classtmap_1_1LSHForest_1a9ee2595fb0f85d917989234ab4aaee8d)`()` <a id="classtmap_1_1LSHForest_1a9ee2595fb0f85d917989234ab4aaee8d"></a>

Remove all entries and the index from the LSH forest.

#### `public size_t `[`size`](#classtmap_1_1LSHForest_1a8ba5c1f500e915c6717c64ac24744874)`()` <a id="classtmap_1_1LSHForest_1a8ba5c1f500e915c6717c64ac24744874"></a>

Get the number of entries.

#### Returns
size_t

# class `tmap::Minhash` <a id="classtmap_1_1Minhash"></a>

An implementation of MinHash and weighted MinHash using SHA1.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public  `[`Minhash`](#classtmap_1_1Minhash_1ae35f57de5ec2316538384c9d2f588d52)`(unsigned int d,unsigned int seed,unsigned int sample_size)` | Construct a new [Minhash](#classtmap_1_1Minhash) object.
`public std::vector< uint32_t > `[`FromBinaryArray`](#classtmap_1_1Minhash_1a1418049bb8c8f70255c336e58a9b9fec)`(std::vector< uint8_t > & vec)` | Create a MinHash from a binary array.
`public std::vector< std::vector< uint32_t > > `[`BatchFromBinaryArray`](#classtmap_1_1Minhash_1a083c1328a9830ec585613c213b1730cc)`(std::vector< std::vector< uint8_t >> & vecs)` | Create MinHashes from a batch of binary arrays (parallelized).
`public std::vector< uint32_t > `[`FromSparseBinaryArray`](#classtmap_1_1Minhash_1aec48525d1c8006f573b0c534e53d894a)`(std::vector< uint32_t > & vec)` | Create a MinHash from a sparse binary array (values are the indices of 1s).
`public std::vector< std::vector< uint32_t > > `[`BatchFromSparseBinaryArray`](#classtmap_1_1Minhash_1a490cf682e7445393fcf2908d74498ea5)`(std::vector< std::vector< uint32_t >> & vecs)` | Create MinHashes from a vector of sparse binary arrays (values are the indices of 1s) (parallelized).
`public std::vector< uint32_t > `[`FromStringArray`](#classtmap_1_1Minhash_1ab21e92280c7265a8df9477734361b8fc)`(std::vector< std::string > & vec)` | Create a MinHash from an array of strings.
`public std::vector< std::vector< uint32_t > > `[`BatchFromStringArray`](#classtmap_1_1Minhash_1a9382e443b9f622c4564449373051d006)`(std::vector< std::vector< std::string >> & vecs)` | Create MinHashes from a vector of string arrays (parallelized).
`public std::vector< uint32_t > `[`FromWeightArray`](#classtmap_1_1Minhash_1ac77f5302d479a2bc2c23a2304a9cf049)`(std::vector< float > & vec)` | Create a MinHash from an array containing weights.
`public std::vector< std::vector< uint32_t > > `[`BatchFromWeightArray`](#classtmap_1_1Minhash_1a20feda993a498b2d92b8e81ca71f73a9)`(std::vector< std::vector< float >> & vecs)` | Create MinHashes from a vector of weight arrays (parallelized).
`public std::vector< uint8_t > `[`ExpandIntWeightArray`](#classtmap_1_1Minhash_1af515f50f6724d7fb16d39743d7652863)`(std::vector< uint32_t > & vec,std::vector< uint32_t > & max_vec,uint32_t size)` | Expand a integer weight array into a binary array.
`public std::vector< std::vector< uint32_t > > `[`BatchFromIntWeightArray`](#classtmap_1_1Minhash_1ac0112cf3a99b2e882803a07dc0f0a620)`(std::vector< std::vector< uint32_t >> & vecs)` | Create MinHashes from a expanded integer weight array (parallelized).
`public float `[`GetDistance`](#classtmap_1_1Minhash_1a21df254dd86462a1dcbe45285c747e71)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` | Get the distance between two MinHashes.
`public float `[`GetWeightedDistance`](#classtmap_1_1Minhash_1a7a8090c1629a6783fe0e17b227bd59ca)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` | Get the weighted distance between two MinHashes.
`public inline  `[`~Minhash`](#classtmap_1_1Minhash_1ae5e1e056f6a1179651b7fa6c196bf220)`()` | Destroy the [Minhash](#classtmap_1_1Minhash) object.

## Members

#### `public  `[`Minhash`](#classtmap_1_1Minhash_1ae35f57de5ec2316538384c9d2f588d52)`(unsigned int d,unsigned int seed,unsigned int sample_size)` <a id="classtmap_1_1Minhash_1ae35f57de5ec2316538384c9d2f588d52"></a>

Construct a new [Minhash](#classtmap_1_1Minhash) object.

#### Parameters
* `d` The number of permutations used for hashing. 

* `seed` The seed for the random number generator. 

* `sample_size` The sample size when generating a weighted MinHash.

#### `public std::vector< uint32_t > `[`FromBinaryArray`](#classtmap_1_1Minhash_1a1418049bb8c8f70255c336e58a9b9fec)`(std::vector< uint8_t > & vec)` <a id="classtmap_1_1Minhash_1a1418049bb8c8f70255c336e58a9b9fec"></a>

Create a MinHash from a binary array.

#### Parameters
* `vec` A vector containing binary values. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromBinaryArray`](#classtmap_1_1Minhash_1a083c1328a9830ec585613c213b1730cc)`(std::vector< std::vector< uint8_t >> & vecs)` <a id="classtmap_1_1Minhash_1a083c1328a9830ec585613c213b1730cc"></a>

Create MinHashes from a batch of binary arrays (parallelized).

#### Parameters
* `vecs` A vector of a vector containing binary values. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint32_t > `[`FromSparseBinaryArray`](#classtmap_1_1Minhash_1aec48525d1c8006f573b0c534e53d894a)`(std::vector< uint32_t > & vec)` <a id="classtmap_1_1Minhash_1aec48525d1c8006f573b0c534e53d894a"></a>

Create a MinHash from a sparse binary array (values are the indices of 1s).

#### Parameters
* `vec` A vector of indices. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromSparseBinaryArray`](#classtmap_1_1Minhash_1a490cf682e7445393fcf2908d74498ea5)`(std::vector< std::vector< uint32_t >> & vecs)` <a id="classtmap_1_1Minhash_1a490cf682e7445393fcf2908d74498ea5"></a>

Create MinHashes from a vector of sparse binary arrays (values are the indices of 1s) (parallelized).

#### Parameters
* `vecs` A vector of vectors of indices. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint32_t > `[`FromStringArray`](#classtmap_1_1Minhash_1ab21e92280c7265a8df9477734361b8fc)`(std::vector< std::string > & vec)` <a id="classtmap_1_1Minhash_1ab21e92280c7265a8df9477734361b8fc"></a>

Create a MinHash from an array of strings.

#### Parameters
* `vec` A vector of strings. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromStringArray`](#classtmap_1_1Minhash_1a9382e443b9f622c4564449373051d006)`(std::vector< std::vector< std::string >> & vecs)` <a id="classtmap_1_1Minhash_1a9382e443b9f622c4564449373051d006"></a>

Create MinHashes from a vector of string arrays (parallelized).

#### Parameters
* `vecs` A vector of string vectors. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint32_t > `[`FromWeightArray`](#classtmap_1_1Minhash_1ac77f5302d479a2bc2c23a2304a9cf049)`(std::vector< float > & vec)` <a id="classtmap_1_1Minhash_1ac77f5302d479a2bc2c23a2304a9cf049"></a>

Create a MinHash from an array containing weights.

#### Parameters
* `vec` A vector of float weights. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromWeightArray`](#classtmap_1_1Minhash_1a20feda993a498b2d92b8e81ca71f73a9)`(std::vector< std::vector< float >> & vecs)` <a id="classtmap_1_1Minhash_1a20feda993a498b2d92b8e81ca71f73a9"></a>

Create MinHashes from a vector of weight arrays (parallelized).

#### Parameters
* `vecs` A vector of float vector weights. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint8_t > `[`ExpandIntWeightArray`](#classtmap_1_1Minhash_1af515f50f6724d7fb16d39743d7652863)`(std::vector< uint32_t > & vec,std::vector< uint32_t > & max_vec,uint32_t size)` <a id="classtmap_1_1Minhash_1af515f50f6724d7fb16d39743d7652863"></a>

Expand a integer weight array into a binary array.

#### Parameters
* `vec` A vector containing integer weights. 

* `max_vec` The maxima for all columns. 

* `size` The size of the expanded array. 

#### Returns
std::vector<uint8_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromIntWeightArray`](#classtmap_1_1Minhash_1ac0112cf3a99b2e882803a07dc0f0a620)`(std::vector< std::vector< uint32_t >> & vecs)` <a id="classtmap_1_1Minhash_1ac0112cf3a99b2e882803a07dc0f0a620"></a>

Create MinHashes from a expanded integer weight array (parallelized).

#### Parameters
* `vecs` A vector of expanded integer weight vectors. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public float `[`GetDistance`](#classtmap_1_1Minhash_1a21df254dd86462a1dcbe45285c747e71)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` <a id="classtmap_1_1Minhash_1a21df254dd86462a1dcbe45285c747e71"></a>

Get the distance between two MinHashes.

#### Parameters
* `vec_a` A MinHash vector. 

* `vec_b` A MinHash vector. 

#### Returns
float

#### `public float `[`GetWeightedDistance`](#classtmap_1_1Minhash_1a7a8090c1629a6783fe0e17b227bd59ca)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` <a id="classtmap_1_1Minhash_1a7a8090c1629a6783fe0e17b227bd59ca"></a>

Get the weighted distance between two MinHashes.

#### Parameters
* `vec_a` A weighted MinHash vector. 

* `vec_b` A weighted MinHash vector. 

#### Returns
float

#### `public inline  `[`~Minhash`](#classtmap_1_1Minhash_1ae5e1e056f6a1179651b7fa6c196bf220)`()` <a id="classtmap_1_1Minhash_1ae5e1e056f6a1179651b7fa6c196bf220"></a>

Destroy the [Minhash](#classtmap_1_1Minhash) object.

# class `tmap::Timer` <a id="classtmap_1_1Timer"></a>

A simple timer class used to check performance during development.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public inline  `[`Timer`](#classtmap_1_1Timer_1ac6eca9fca91d67ee83361136abd2ef3b)`()` | Construct a new [Timer](#classtmap_1_1Timer) object and start the clock.
`public inline void `[`reset`](#classtmap_1_1Timer_1acf16459c9388f5a3702c934e6135edf3)`()` | Restart the clock.
`public inline double `[`elapsed`](#classtmap_1_1Timer_1a34a0281aaced5d7768c5bc60ebcd0751)`() const` | Return the time elapsed since the timer was started or last reset.

## Members

#### `public inline  `[`Timer`](#classtmap_1_1Timer_1ac6eca9fca91d67ee83361136abd2ef3b)`()` <a id="classtmap_1_1Timer_1ac6eca9fca91d67ee83361136abd2ef3b"></a>

Construct a new [Timer](#classtmap_1_1Timer) object and start the clock.

#### `public inline void `[`reset`](#classtmap_1_1Timer_1acf16459c9388f5a3702c934e6135edf3)`()` <a id="classtmap_1_1Timer_1acf16459c9388f5a3702c934e6135edf3"></a>

Restart the clock.

#### `public inline double `[`elapsed`](#classtmap_1_1Timer_1a34a0281aaced5d7768c5bc60ebcd0751)`() const` <a id="classtmap_1_1Timer_1a34a0281aaced5d7768c5bc60ebcd0751"></a>

Return the time elapsed since the timer was started or last reset.

#### Returns
double

# struct `tmap::GraphProperties` <a id="structtmap_1_1GraphProperties"></a>

The properties of a generated graph. An instance of this struct is returned from the layout functions.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public float `[`mst_weight`](#structtmap_1_1GraphProperties_1a51fa52c1cc9c8d382d507349e99447be) | The total weight of the created spanning tree.
`public uint32_t `[`n_connected_components`](#structtmap_1_1GraphProperties_1a3467571c1e645268e55416504809f7c5) | The number of connected components.
`public uint32_t `[`n_isolated_vertices`](#structtmap_1_1GraphProperties_1a954e8fd087b44b3c568ea07ae2f1efea) | The number of isolated (lone) vertices.
`public std::vector< uint32_t > `[`degrees`](#structtmap_1_1GraphProperties_1af4c85653b3bf56c6dbf19d7a38af40bc) | The degrees of the vertices in the graph.
`public std::vector< std::vector< uint32_t > > `[`adjacency_list`](#structtmap_1_1GraphProperties_1a04c11168f810fdaf8b7bfecf96413ba7) | The adjacency list of the spanning tree.

## Members

#### `public float `[`mst_weight`](#structtmap_1_1GraphProperties_1a51fa52c1cc9c8d382d507349e99447be) <a id="structtmap_1_1GraphProperties_1a51fa52c1cc9c8d382d507349e99447be"></a>

The total weight of the created spanning tree.

#### `public uint32_t `[`n_connected_components`](#structtmap_1_1GraphProperties_1a3467571c1e645268e55416504809f7c5) <a id="structtmap_1_1GraphProperties_1a3467571c1e645268e55416504809f7c5"></a>

The number of connected components.

#### `public uint32_t `[`n_isolated_vertices`](#structtmap_1_1GraphProperties_1a954e8fd087b44b3c568ea07ae2f1efea) <a id="structtmap_1_1GraphProperties_1a954e8fd087b44b3c568ea07ae2f1efea"></a>

The number of isolated (lone) vertices.

#### `public std::vector< uint32_t > `[`degrees`](#structtmap_1_1GraphProperties_1af4c85653b3bf56c6dbf19d7a38af40bc) <a id="structtmap_1_1GraphProperties_1af4c85653b3bf56c6dbf19d7a38af40bc"></a>

The degrees of the vertices in the graph.

#### `public std::vector< std::vector< uint32_t > > `[`adjacency_list`](#structtmap_1_1GraphProperties_1a04c11168f810fdaf8b7bfecf96413ba7) <a id="structtmap_1_1GraphProperties_1a04c11168f810fdaf8b7bfecf96413ba7"></a>

The adjacency list of the spanning tree.

# struct `tmap::LayoutConfiguration` <a id="structtmap_1_1LayoutConfiguration"></a>

A struct containing all the configuration options available for and applied to a layout.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public int `[`k`](#structtmap_1_1LayoutConfiguration_1a80ddc818732d708764fbd83ad7b7d153) | The number of nearest neighbors used to create the k-nearest neighbor graph.
`public int `[`kc`](#structtmap_1_1LayoutConfiguration_1ae63c0a1d5956cbdb837f3aff1978f867) | The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned.
`public int `[`fme_iterations`](#structtmap_1_1LayoutConfiguration_1aa0b26a532aedb8f0fe4490c9c90b0e84) | Maximum number of iterations of the fast multipole embedder.
`public bool `[`fme_randomize`](#structtmap_1_1LayoutConfiguration_1a821bf612fb3063344ca0c6b161424a7b) | Whether or not to randomize the layout at the start.
`public int `[`fme_threads`](#structtmap_1_1LayoutConfiguration_1a7630a7d7513c3f51ea00802a3f67ba92) | The number of threads for the fast multipole embedder.
`public int `[`fme_precision`](#structtmap_1_1LayoutConfiguration_1a9e4d43d8f65c21404cc9912c11a3eba7) | The number of coefficients of the multipole expansion.
`public int `[`sl_repeats`](#structtmap_1_1LayoutConfiguration_1adf81cfbcba521fd87a73dd25eb9c21e7) | The number of repeats of the scaling layout algorithm.
`public int `[`sl_extra_scaling_steps`](#structtmap_1_1LayoutConfiguration_1aefd713cfba563ea8ee9a432c0359c440) | Sets the number of repeats of the scaling.
`public double `[`sl_scaling_min`](#structtmap_1_1LayoutConfiguration_1af2e01075c5fe2a36c6018b40a6919d7a) | The minimum scaling factor.
`public double `[`sl_scaling_max`](#structtmap_1_1LayoutConfiguration_1aa29669e99b7e1df2bcb95d6379da7cc3) | The maximum scaling factor.
`public `[`ScalingType`](#layout_8hh_1a50ec215c9e54cf12b9dd0a0056160761)` `[`sl_scaling_type`](#structtmap_1_1LayoutConfiguration_1a618d286e035eca76e0e464513624beec) | Defines the (relative) scale of the graph.
`public int `[`mmm_repeats`](#structtmap_1_1LayoutConfiguration_1aff2347eb71c98bbc72f16b4de32d4af0) | Number of repeats of the per-level layout algorithm.
`public `[`Placer`](#layout_8hh_1afdc98947e81dc6f4c30f256e6f42f90b)` `[`placer`](#structtmap_1_1LayoutConfiguration_1ae81108ee33f42b2c084b540f902bbb7d) | The method by which the initial positons of the vertices at eachlevel are defined.
`public `[`Merger`](#layout_8hh_1a8c7bb9956a1a724233182a166cfdc0ff)` `[`merger`](#structtmap_1_1LayoutConfiguration_1aeee45308fd8dbda38fbc7b8c7ff9212f) | The vertex merging strategy applied during the coarsening phaseof the multilevel algorithm.
`public double `[`merger_factor`](#structtmap_1_1LayoutConfiguration_1a72fe4f8f738d2d400f70db97c4273a46) | The ratio of the sizes between two levels up to which the mergingis run. Does not apply to all merging strategies.
`public int `[`merger_adjustment`](#structtmap_1_1LayoutConfiguration_1a16109420c8ec0a4c3021345fd943daf6) | The edge length adjustment of the merging algorithm. Does notapply to all merging strategies.
`public float `[`node_size`](#structtmap_1_1LayoutConfiguration_1a9a97e2c0c9edb212190d3afcc3ce2924) | The size of the nodes, which affects the magnitude of their repellingforce. Decreasing this value generally resolves overlaps in a verycrowded tree.
`public inline  `[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration_1a45335a69efe4408b49283554a3bb8875)`()` | Construct a new Layout Configuration object.
`public inline std::string `[`ToString`](#structtmap_1_1LayoutConfiguration_1a498341508ea4806795f44e376af18e11)`() const` | Returns a string describing the set options.

## Members

#### `public int `[`k`](#structtmap_1_1LayoutConfiguration_1a80ddc818732d708764fbd83ad7b7d153) <a id="structtmap_1_1LayoutConfiguration_1a80ddc818732d708764fbd83ad7b7d153"></a>

The number of nearest neighbors used to create the k-nearest neighbor graph.

#### `public int `[`kc`](#structtmap_1_1LayoutConfiguration_1ae63c0a1d5956cbdb837f3aff1978f867) <a id="structtmap_1_1LayoutConfiguration_1ae63c0a1d5956cbdb837f3aff1978f867"></a>

The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned.

#### `public int `[`fme_iterations`](#structtmap_1_1LayoutConfiguration_1aa0b26a532aedb8f0fe4490c9c90b0e84) <a id="structtmap_1_1LayoutConfiguration_1aa0b26a532aedb8f0fe4490c9c90b0e84"></a>

Maximum number of iterations of the fast multipole embedder.

#### `public bool `[`fme_randomize`](#structtmap_1_1LayoutConfiguration_1a821bf612fb3063344ca0c6b161424a7b) <a id="structtmap_1_1LayoutConfiguration_1a821bf612fb3063344ca0c6b161424a7b"></a>

Whether or not to randomize the layout at the start.

#### `public int `[`fme_threads`](#structtmap_1_1LayoutConfiguration_1a7630a7d7513c3f51ea00802a3f67ba92) <a id="structtmap_1_1LayoutConfiguration_1a7630a7d7513c3f51ea00802a3f67ba92"></a>

The number of threads for the fast multipole embedder.

#### `public int `[`fme_precision`](#structtmap_1_1LayoutConfiguration_1a9e4d43d8f65c21404cc9912c11a3eba7) <a id="structtmap_1_1LayoutConfiguration_1a9e4d43d8f65c21404cc9912c11a3eba7"></a>

The number of coefficients of the multipole expansion.

#### `public int `[`sl_repeats`](#structtmap_1_1LayoutConfiguration_1adf81cfbcba521fd87a73dd25eb9c21e7) <a id="structtmap_1_1LayoutConfiguration_1adf81cfbcba521fd87a73dd25eb9c21e7"></a>

The number of repeats of the scaling layout algorithm.

#### `public int `[`sl_extra_scaling_steps`](#structtmap_1_1LayoutConfiguration_1aefd713cfba563ea8ee9a432c0359c440) <a id="structtmap_1_1LayoutConfiguration_1aefd713cfba563ea8ee9a432c0359c440"></a>

Sets the number of repeats of the scaling.

#### `public double `[`sl_scaling_min`](#structtmap_1_1LayoutConfiguration_1af2e01075c5fe2a36c6018b40a6919d7a) <a id="structtmap_1_1LayoutConfiguration_1af2e01075c5fe2a36c6018b40a6919d7a"></a>

The minimum scaling factor.

#### `public double `[`sl_scaling_max`](#structtmap_1_1LayoutConfiguration_1aa29669e99b7e1df2bcb95d6379da7cc3) <a id="structtmap_1_1LayoutConfiguration_1aa29669e99b7e1df2bcb95d6379da7cc3"></a>

The maximum scaling factor.

#### `public `[`ScalingType`](#layout_8hh_1a50ec215c9e54cf12b9dd0a0056160761)` `[`sl_scaling_type`](#structtmap_1_1LayoutConfiguration_1a618d286e035eca76e0e464513624beec) <a id="structtmap_1_1LayoutConfiguration_1a618d286e035eca76e0e464513624beec"></a>

Defines the (relative) scale of the graph.

#### `public int `[`mmm_repeats`](#structtmap_1_1LayoutConfiguration_1aff2347eb71c98bbc72f16b4de32d4af0) <a id="structtmap_1_1LayoutConfiguration_1aff2347eb71c98bbc72f16b4de32d4af0"></a>

Number of repeats of the per-level layout algorithm.

#### `public `[`Placer`](#layout_8hh_1afdc98947e81dc6f4c30f256e6f42f90b)` `[`placer`](#structtmap_1_1LayoutConfiguration_1ae81108ee33f42b2c084b540f902bbb7d) <a id="structtmap_1_1LayoutConfiguration_1ae81108ee33f42b2c084b540f902bbb7d"></a>

The method by which the initial positons of the vertices at eachlevel are defined.

#### `public `[`Merger`](#layout_8hh_1a8c7bb9956a1a724233182a166cfdc0ff)` `[`merger`](#structtmap_1_1LayoutConfiguration_1aeee45308fd8dbda38fbc7b8c7ff9212f) <a id="structtmap_1_1LayoutConfiguration_1aeee45308fd8dbda38fbc7b8c7ff9212f"></a>

The vertex merging strategy applied during the coarsening phaseof the multilevel algorithm.

#### `public double `[`merger_factor`](#structtmap_1_1LayoutConfiguration_1a72fe4f8f738d2d400f70db97c4273a46) <a id="structtmap_1_1LayoutConfiguration_1a72fe4f8f738d2d400f70db97c4273a46"></a>

The ratio of the sizes between two levels up to which the mergingis run. Does not apply to all merging strategies.

#### `public int `[`merger_adjustment`](#structtmap_1_1LayoutConfiguration_1a16109420c8ec0a4c3021345fd943daf6) <a id="structtmap_1_1LayoutConfiguration_1a16109420c8ec0a4c3021345fd943daf6"></a>

The edge length adjustment of the merging algorithm. Does notapply to all merging strategies.

#### `public float `[`node_size`](#structtmap_1_1LayoutConfiguration_1a9a97e2c0c9edb212190d3afcc3ce2924) <a id="structtmap_1_1LayoutConfiguration_1a9a97e2c0c9edb212190d3afcc3ce2924"></a>

The size of the nodes, which affects the magnitude of their repellingforce. Decreasing this value generally resolves overlaps in a verycrowded tree.

#### `public inline  `[`LayoutConfiguration`](#structtmap_1_1LayoutConfiguration_1a45335a69efe4408b49283554a3bb8875)`()` <a id="structtmap_1_1LayoutConfiguration_1a45335a69efe4408b49283554a3bb8875"></a>

Construct a new Layout Configuration object.

#### `public inline std::string `[`ToString`](#structtmap_1_1LayoutConfiguration_1a498341508ea4806795f44e376af18e11)`() const` <a id="structtmap_1_1LayoutConfiguration_1a498341508ea4806795f44e376af18e11"></a>

Returns a string describing the set options.

#### Returns
std::string

# struct `tmap::MapKeyPointer` <a id="structtmap_1_1MapKeyPointer"></a>

The pointer map used for pointing to the keys from the sorted hash map.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public iterator `[`it`](#structtmap_1_1MapKeyPointer_1a7e5f9fcf74a41f693afd2a727a2997d5) | 
`public inline  `[`MapKeyPointer`](#structtmap_1_1MapKeyPointer_1ac96f2c30e0923aa2b24674d1d0a8951f)`(iterator i)` | 
`public inline  `[`MapKeyPointer`](#structtmap_1_1MapKeyPointer_1ad69cc205a77366ea584d3b273798dd0d)`()` | 
`public inline const std::vector< uint8_t > & `[`operator*`](#structtmap_1_1MapKeyPointer_1ade7f29bcf8ac9ec5875f1152d50912e4)`() const` | 
`public inline const std::vector< uint8_t > * `[`operator->`](#structtmap_1_1MapKeyPointer_1ab681e8610b86e353af4feb6ea0ef8a9a)`() const` | 
`typedef `[`iterator`](#structtmap_1_1MapKeyPointer_1a236dc7396d6a4dca2a942196aa536f56) | 

## Members

#### `public iterator `[`it`](#structtmap_1_1MapKeyPointer_1a7e5f9fcf74a41f693afd2a727a2997d5) <a id="structtmap_1_1MapKeyPointer_1a7e5f9fcf74a41f693afd2a727a2997d5"></a>

#### `public inline  `[`MapKeyPointer`](#structtmap_1_1MapKeyPointer_1ac96f2c30e0923aa2b24674d1d0a8951f)`(iterator i)` <a id="structtmap_1_1MapKeyPointer_1ac96f2c30e0923aa2b24674d1d0a8951f"></a>

#### `public inline  `[`MapKeyPointer`](#structtmap_1_1MapKeyPointer_1ad69cc205a77366ea584d3b273798dd0d)`()` <a id="structtmap_1_1MapKeyPointer_1ad69cc205a77366ea584d3b273798dd0d"></a>

#### `public inline const std::vector< uint8_t > & `[`operator*`](#structtmap_1_1MapKeyPointer_1ade7f29bcf8ac9ec5875f1152d50912e4)`() const` <a id="structtmap_1_1MapKeyPointer_1ade7f29bcf8ac9ec5875f1152d50912e4"></a>

#### `public inline const std::vector< uint8_t > * `[`operator->`](#structtmap_1_1MapKeyPointer_1ab681e8610b86e353af4feb6ea0ef8a9a)`() const` <a id="structtmap_1_1MapKeyPointer_1ab681e8610b86e353af4feb6ea0ef8a9a"></a>

#### `typedef `[`iterator`](#structtmap_1_1MapKeyPointer_1a236dc7396d6a4dca2a942196aa536f56) <a id="structtmap_1_1MapKeyPointer_1a236dc7396d6a4dca2a942196aa536f56"></a>

# struct `tmap::SimpleHash` <a id="structtmap_1_1SimpleHash"></a>

Hash struct used for the sparsepp sparse hash map.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public inline size_t `[`operator()`](#structtmap_1_1SimpleHash_1a2246e182ae49e4ff363d6c3344aa97cf)`(std::vector< uint8_t > vec) const` | 

## Members

#### `public inline size_t `[`operator()`](#structtmap_1_1SimpleHash_1a2246e182ae49e4ff363d6c3344aa97cf)`(std::vector< uint8_t > vec) const` <a id="structtmap_1_1SimpleHash_1a2246e182ae49e4ff363d6c3344aa97cf"></a>

Generated by [Moxygen](https://sourcey.com/moxygen)