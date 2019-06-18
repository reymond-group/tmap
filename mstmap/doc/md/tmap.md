# Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`namespace `[`sha1`](#namespacesha1) | 
`namespace `[`tmap`](#namespacetmap) | 
`class `[`LSHForest`](#classLSHForest) | Provides locality sensitive hashing forest functionalities.
`class `[`Minhash`](#classMinhash) | An implementation of MinHash and weighted MinHash using SHA1.
`class `[`Timer`](#classTimer) | A simple timer class used to check performance during development.
`struct `[`GraphProperties`](#structGraphProperties) | The properties of a generated graph. An instance of this struct is returned from the layout functions.
`struct `[`LayoutConfiguration`](#structLayoutConfiguration) | A struct containing all the configuration options available for and applied to a layout.
`struct `[`MapKeyPointer`](#structMapKeyPointer) | The pointer map used for pointing to the keys from the sorted hash map.
`struct `[`SimpleHash`](#structSimpleHash) | Hash struct used for the sparsepp sparse hash map.

# namespace `sha1` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`class `[`sha1::SHA1`](#classsha1_1_1SHA1) | 

# class `sha1::SHA1` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public inline  `[`SHA1`](#classsha1_1_1SHA1_1a440d0c3f17da510418b678f4c9dcbdf2)`()` | 
`public inline virtual  `[`~SHA1`](#classsha1_1_1SHA1_1a7103b57eb1f51a498a92c095b9c21392)`()` | 
`public inline  `[`SHA1`](#classsha1_1_1SHA1_1a1d63d04a7f146a42000e0a2a35314c6f)`(const `[`SHA1`](#classsha1_1_1SHA1)` & s)` | 
`public inline const `[`SHA1`](#classsha1_1_1SHA1)` & `[`operator=`](#classsha1_1_1SHA1_1ae8f2ccbee9a15968482057e0d03a4247)`(const `[`SHA1`](#classsha1_1_1SHA1)` & s)` | 
`public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`reset`](#classsha1_1_1SHA1_1ae2f5530e5be0837e6a222a1f244cbe6d)`()` | 
`public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`processByte`](#classsha1_1_1SHA1_1a8f7edf8ccfea07f988d5f40c4b5c768b)`(uint8_t octet)` | 
`public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`processBlock`](#classsha1_1_1SHA1_1af5e58568ad3bcaf91640c7fd5e448d85)`(const void *const start,const void *const end)` | 
`public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`processBytes`](#classsha1_1_1SHA1_1abf8a680503b98b8f773f98752e34dbca)`(const void *const data,size_t len)` | 
`public inline const uint32_t * `[`getDigest`](#classsha1_1_1SHA1_1a7a097712264451da60aa0f30096dd516)`(digest32_t digest)` | 
`public inline const uint8_t * `[`getDigestBytes`](#classsha1_1_1SHA1_1a4fd38b624cba611b65730a9e747c584e)`(digest8_t digest)` | 
`protected inline void `[`processBlock`](#classsha1_1_1SHA1_1af00112a350d5f1f67a9ac00a0fee6793)`()` | 
`typedef `[`digest32_t`](#classsha1_1_1SHA1_1a15f384f39d235a8912d5042dc920595f) | 
`typedef `[`digest8_t`](#classsha1_1_1SHA1_1a62b6c7838c4cdcb81700a6cc64fde994) | 

## Members

#### `public inline  `[`SHA1`](#classsha1_1_1SHA1_1a440d0c3f17da510418b678f4c9dcbdf2)`()` 

#### `public inline virtual  `[`~SHA1`](#classsha1_1_1SHA1_1a7103b57eb1f51a498a92c095b9c21392)`()` 

#### `public inline  `[`SHA1`](#classsha1_1_1SHA1_1a1d63d04a7f146a42000e0a2a35314c6f)`(const `[`SHA1`](#classsha1_1_1SHA1)` & s)` 

#### `public inline const `[`SHA1`](#classsha1_1_1SHA1)` & `[`operator=`](#classsha1_1_1SHA1_1ae8f2ccbee9a15968482057e0d03a4247)`(const `[`SHA1`](#classsha1_1_1SHA1)` & s)` 

#### `public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`reset`](#classsha1_1_1SHA1_1ae2f5530e5be0837e6a222a1f244cbe6d)`()` 

#### `public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`processByte`](#classsha1_1_1SHA1_1a8f7edf8ccfea07f988d5f40c4b5c768b)`(uint8_t octet)` 

#### `public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`processBlock`](#classsha1_1_1SHA1_1af5e58568ad3bcaf91640c7fd5e448d85)`(const void *const start,const void *const end)` 

#### `public inline `[`SHA1`](#classsha1_1_1SHA1)` & `[`processBytes`](#classsha1_1_1SHA1_1abf8a680503b98b8f773f98752e34dbca)`(const void *const data,size_t len)` 

#### `public inline const uint32_t * `[`getDigest`](#classsha1_1_1SHA1_1a7a097712264451da60aa0f30096dd516)`(digest32_t digest)` 

#### `public inline const uint8_t * `[`getDigestBytes`](#classsha1_1_1SHA1_1a4fd38b624cba611b65730a9e747c584e)`(digest8_t digest)` 

#### `protected inline void `[`processBlock`](#classsha1_1_1SHA1_1af00112a350d5f1f67a9ac00a0fee6793)`()` 

#### `typedef `[`digest32_t`](#classsha1_1_1SHA1_1a15f384f39d235a8912d5042dc920595f) 

#### `typedef `[`digest8_t`](#classsha1_1_1SHA1_1a62b6c7838c4cdcb81700a6cc64fde994) 

# namespace `tmap` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public def `[`get_asset`](#____init_____8py_1a4bc5750c95a4e2269b40c06b1d07da30)`(path)`            | 

## Members

#### `public def `[`get_asset`](#____init_____8py_1a4bc5750c95a4e2269b40c06b1d07da30)`(path)` 

# class `LSHForest` 

Provides locality sensitive hashing forest functionalities.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public  `[`LSHForest`](#classLSHForest_1ae227bb302c173481a129335fa581fa6f)`(unsigned int d,unsigned int l,bool store,bool file_backed)` | Construct a new [LSHForest](#classLSHForest) object.
`public inline  `[`~LSHForest`](#classLSHForest_1ac2e048b53de7268b2b16da47b99a3955)`()` | Destroy the [LSHForest](#classLSHForest) object.
`public void `[`Add`](#classLSHForest_1a000536ff09a94df528a1f72ffa40bac3)`(std::vector< uint32_t > & vec)` | Add a MinHash to this [LSHForest](#classLSHForest) instance.
`public void `[`BatchAdd`](#classLSHForest_1a668dc958bfc4856f26a96a4aa0897e53)`(std::vector< std::vector< uint32_t >> & vecs)` | Add Minhashes to this [LSHForest](#classLSHForest) (parallelized).
`public void `[`Index`](#classLSHForest_1afdd95ab82907622f9a38a386766a53d2)`()` | Create the index (trees).
`public bool `[`IsClean`](#classLSHForest_1af7cbc237713af39831aeab21541dafcf)`()` | Check whether the added MinHashes have been indexed.
`public void `[`Store`](#classLSHForest_1a758c5329128f5ab3723a98b77bbc4634)`(const std::string & path)` | Write / serialize the current LSH forest to the disk.
`public void `[`Restore`](#classLSHForest_1aae4d0534e33a6ac46b56eeb7630916f6)`(const std::string & path)` | Read / deserialize a LSH forest instance form the disk. The forest is indexed automatically.
`public std::vector< uint32_t > `[`GetHash`](#classLSHForest_1aec42a8cb3d1caf12faeb8d0f9ea09529)`(uint32_t id)` | Get the MinHash of an entry at a given index. The index is defined by order of insertion.
`public void `[`GetKNNGraph`](#classLSHForest_1a6c84d67e979b5bfd6123e62e1350dbf2)`(std::vector< uint32_t > & from,std::vector< uint32_t > & to,std::vector< float > & weight,unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbor graph of the data stored in this LSH forest instance. It will be written to out parameters as an edge list.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScan`](#classLSHForest_1aa655ed6c39050b45b73a051abe51e035)`(const std::vector< uint32_t > & vec,unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of a query.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExclude`](#classLSHForest_1a932c426296cbd6da0e7c29cd4212f3ff)`(const std::vector< uint32_t > & vec,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of a query except those defined in the argument exclude.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanById`](#classLSHForest_1afe623496f801357e8f555259620bc174)`(uint32_t id,unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of an entry.
`public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExcludeById`](#classLSHForest_1a5cdb395444dad71ba95e2000bb896454)`(uint32_t id,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of an entry except those defined in the argument exclude.
`public std::vector< std::pair< float, uint32_t > > `[`LinearScan`](#classLSHForest_1a5d5b1675caaa17d9c2a4ba8c95c645a3)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & indices,unsigned int k,bool weighted)` | Get the k-nearest neighbors of a query using linear scan.
`public std::vector< uint32_t > `[`Query`](#classLSHForest_1a6bc39aa54083ede4ab9ba1b0f12c7229)`(const std::vector< uint32_t > & vec,unsigned int k)` | Query the LSH forest for k-nearest neighbors.
`public std::vector< uint32_t > `[`QueryExclude`](#classLSHForest_1ada7ea3fd5c3eb9fc05188a0054de48cf)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & exclude,unsigned int k)` | Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.
`public std::vector< uint32_t > `[`QueryById`](#classLSHForest_1ade573cce99526ba05341dd506673ea8b)`(uint32_t id,unsigned int k)` | Query the LSH forest for k-nearest neighbors.
`public std::vector< uint32_t > `[`QueryExcludeById`](#classLSHForest_1a50da7a1db11f709c54e5e05fd5b08aa8)`(uint32_t id,std::vector< uint32_t > & exclude,unsigned int k)` | Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.
`public std::vector< std::vector< uint32_t > > `[`BatchQuery`](#classLSHForest_1a23e5fd430b95580e09126ba58bde32a4)`(const std::vector< std::vector< uint32_t >> & vecs,unsigned int k)` | Query the LSH forest for k-nearest neighbors (parallelized).
`public std::vector< uint32_t > `[`GetAllNearestNeighbors`](#classLSHForest_1ac741709bb8e322c68a5749c023147bde)`(unsigned int k,unsigned int kc,bool weighted)` | Get the k-nearest neighbors of all LSH forest entries.
`public std::vector< uint32_t > `[`GetData`](#classLSHForest_1aad89848405eebc847c18a75b618da3e1)`(uint32_t id)` | Get the MinHash of an entry at a given index. The index is defined by order of insertion. Alias for GetHash.
`public std::vector< float > `[`GetAllDistances`](#classLSHForest_1abdc3bd3357708bb58e6d957376e92f1c)`(const std::vector< uint32_t > & vec)` | Get the distances of a MinHash to all entries in the LSH forest.
`public float `[`GetDistance`](#classLSHForest_1a08d66568664bdc8e0148c18b18a1b8fa)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` | Get the distance between two MinHashes.
`public float `[`GetWeightedDistance`](#classLSHForest_1acfb878c731daf8da6402e7cebc2b6ef1)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` | Get the distance between two weighted MinHashes.
`public float `[`GetDistanceById`](#classLSHForest_1a49ad1fe0429121a8b572cde4df973d96)`(uint32_t a,uint32_t b)` | Get the distance between two MinHashes.
`public float `[`GetWeightedDistanceById`](#classLSHForest_1a72b5a201bc8c409c0c5742587988ba85)`(uint32_t a,uint32_t b)` | Get the distance between two weighted MinHashes.
`public void `[`Clear`](#classLSHForest_1aec34dc5185166dce9c22a4060ce3914c)`()` | Remove all entries and the index from the LSH forest.
`public size_t `[`size`](#classLSHForest_1af0015fae65afd25bf67875d12dc0d663)`()` | Get the number of entries.

## Members

#### `public  `[`LSHForest`](#classLSHForest_1ae227bb302c173481a129335fa581fa6f)`(unsigned int d,unsigned int l,bool store,bool file_backed)` 

Construct a new [LSHForest](#classLSHForest) object.

#### Parameters
* `d` The dimensionality of the MinHashes to be added to this [LSHForest](#classLSHForest). 

* `l` The number of prefix trees used. 

* `store` Whether or not to store the data for later enhanced (using parameter kc) retrievel. 

* `file_backed` Whether to store the data on disk rather than in RAM (experimental).

#### `public inline  `[`~LSHForest`](#classLSHForest_1ac2e048b53de7268b2b16da47b99a3955)`()` 

Destroy the [LSHForest](#classLSHForest) object.

#### `public void `[`Add`](#classLSHForest_1a000536ff09a94df528a1f72ffa40bac3)`(std::vector< uint32_t > & vec)` 

Add a MinHash to this [LSHForest](#classLSHForest) instance.

#### Parameters
* `vec` A MinHash vector.

#### `public void `[`BatchAdd`](#classLSHForest_1a668dc958bfc4856f26a96a4aa0897e53)`(std::vector< std::vector< uint32_t >> & vecs)` 

Add Minhashes to this [LSHForest](#classLSHForest) (parallelized).

#### Parameters
* `vecs` A vector containing MinHash vectors.

#### `public void `[`Index`](#classLSHForest_1afdd95ab82907622f9a38a386766a53d2)`()` 

Create the index (trees).

#### `public bool `[`IsClean`](#classLSHForest_1af7cbc237713af39831aeab21541dafcf)`()` 

Check whether the added MinHashes have been indexed.

#### Returns
true 

#### Returns
false

#### `public void `[`Store`](#classLSHForest_1a758c5329128f5ab3723a98b77bbc4634)`(const std::string & path)` 

Write / serialize the current LSH forest to the disk.

#### Parameters
* `path` The location where the LSH forest should be stored on disk.

#### `public void `[`Restore`](#classLSHForest_1aae4d0534e33a6ac46b56eeb7630916f6)`(const std::string & path)` 

Read / deserialize a LSH forest instance form the disk. The forest is indexed automatically.

#### Parameters
* `path` The location from where to load the LSH forest.

#### `public std::vector< uint32_t > `[`GetHash`](#classLSHForest_1aec42a8cb3d1caf12faeb8d0f9ea09529)`(uint32_t id)` 

Get the MinHash of an entry at a given index. The index is defined by order of insertion.

#### Parameters
* `id` The index (order of insertion) of a entry. 

#### Returns
std::vector<uint32_t> The MinHash associated with an index.

#### `public void `[`GetKNNGraph`](#classLSHForest_1a6c84d67e979b5bfd6123e62e1350dbf2)`(std::vector< uint32_t > & from,std::vector< uint32_t > & to,std::vector< float > & weight,unsigned int k,unsigned int kc,bool weighted)` 

Get the k-nearest neighbor graph of the data stored in this LSH forest instance. It will be written to out parameters as an edge list.

#### Parameters
* `from` A vector to which the from vertices will be written. 

* `to` A vector to which the to vertices will be written. 

* `weight` A vector to which the float weights of the edges will be written. 

* `k` The degree of the nearest neighbor graph. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results are picked to create the k-nearest neighbor graph. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted.

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScan`](#classLSHForest_1aa655ed6c39050b45b73a051abe51e035)`(const std::vector< uint32_t > & vec,unsigned int k,unsigned int kc,bool weighted)` 

Get the k-nearest neighbors of a query.

#### Parameters
* `vec` The query MinHash. 

* `k` The number of k-nearest neighbors to return. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExclude`](#classLSHForest_1a932c426296cbd6da0e7c29cd4212f3ff)`(const std::vector< uint32_t > & vec,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` 

Get the k-nearest neighbors of a query except those defined in the argument exclude.

#### Parameters
* `vec` The query MinHash. 

* `k` The number of k-nearest neighbors to return. 

* `exclude` A list of indices to be excluded from the search 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>>

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanById`](#classLSHForest_1afe623496f801357e8f555259620bc174)`(uint32_t id,unsigned int k,unsigned int kc,bool weighted)` 

Get the k-nearest neighbors of an entry.

#### Parameters
* `id` The id of the query entry. 

* `k` The number of k-nearest neighbors to return. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< std::pair< float, uint32_t > > `[`QueryLinearScanExcludeById`](#classLSHForest_1a5cdb395444dad71ba95e2000bb896454)`(uint32_t id,unsigned int k,std::vector< uint32_t > & exclude,unsigned int kc,bool weighted)` 

Get the k-nearest neighbors of an entry except those defined in the argument exclude.

#### Parameters
* `id` The id of the query entry. 

* `k` The number of k-nearest neighbors to return. 

* `exclude` A list of indices to be excluded from the search. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< std::pair< float, uint32_t > > `[`LinearScan`](#classLSHForest_1a5d5b1675caaa17d9c2a4ba8c95c645a3)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & indices,unsigned int k,bool weighted)` 

Get the k-nearest neighbors of a query using linear scan.

#### Parameters
* `vec` The query MinHash. 

* `indices` A list of indices to on which to run the linear scan. 

* `k` The number of k-nearest neighbors to return. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<std::pair<float, uint32_t>> The distances and indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`Query`](#classLSHForest_1a6bc39aa54083ede4ab9ba1b0f12c7229)`(const std::vector< uint32_t > & vec,unsigned int k)` 

Query the LSH forest for k-nearest neighbors.

#### Parameters
* `vec` The query MinHash. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`QueryExclude`](#classLSHForest_1ada7ea3fd5c3eb9fc05188a0054de48cf)`(const std::vector< uint32_t > & vec,std::vector< uint32_t > & exclude,unsigned int k)` 

Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.

#### Parameters
* `vec` The query MinHash. 

* `exclude` A list of indices to be excluded from the search. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`QueryById`](#classLSHForest_1ade573cce99526ba05341dd506673ea8b)`(uint32_t id,unsigned int k)` 

Query the LSH forest for k-nearest neighbors.

#### Parameters
* `id` The id of the query entry. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`QueryExcludeById`](#classLSHForest_1a50da7a1db11f709c54e5e05fd5b08aa8)`(uint32_t id,std::vector< uint32_t > & exclude,unsigned int k)` 

Query the LSH forest for k-nearest neighbors. Exclude a list of entries by ID.

#### Parameters
* `id` The id of the query entry. 

* `exclude` A list of indices to be excluded from the search. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<uint32_t> The indices of the k-nearest neighbors.

#### `public std::vector< std::vector< uint32_t > > `[`BatchQuery`](#classLSHForest_1a23e5fd430b95580e09126ba58bde32a4)`(const std::vector< std::vector< uint32_t >> & vecs,unsigned int k)` 

Query the LSH forest for k-nearest neighbors (parallelized).

#### Parameters
* `vecs` A vector of MinHashes. 

* `k` The number of nearest neighbors to search for. 

#### Returns
std::vector<std::vector<uint32_t>> A vector of the indices of the k-nearest neighbors.

#### `public std::vector< uint32_t > `[`GetAllNearestNeighbors`](#classLSHForest_1ac741709bb8e322c68a5749c023147bde)`(unsigned int k,unsigned int kc,bool weighted)` 

Get the k-nearest neighbors of all LSH forest entries.

#### Parameters
* `k` The number of nearest neighbors to search for. 

* `kc` The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned. 

* `weighted` Whether the MinHashes contained within this instance of an LSH forest are weighted. 

#### Returns
std::vector<uint32_t> The IDs of the nearest neighbors of all LSH forest entries.

#### `public std::vector< uint32_t > `[`GetData`](#classLSHForest_1aad89848405eebc847c18a75b618da3e1)`(uint32_t id)` 

Get the MinHash of an entry at a given index. The index is defined by order of insertion. Alias for GetHash.

#### Parameters
* `id` The index (order of insertion) of a entry. 

#### Returns
std::vector<uint32_t> The MinHash associated with an index.

#### `public std::vector< float > `[`GetAllDistances`](#classLSHForest_1abdc3bd3357708bb58e6d957376e92f1c)`(const std::vector< uint32_t > & vec)` 

Get the distances of a MinHash to all entries in the LSH forest.

#### Parameters
* `vec` The query MinHash. 

#### Returns
std::vector<float> The distances form the input MinHash to all the entries in the LSH forest.

#### `public float `[`GetDistance`](#classLSHForest_1a08d66568664bdc8e0148c18b18a1b8fa)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` 

Get the distance between two MinHashes.

#### Parameters
* `vec_a` A MinHash. 

* `vec_b` A MinHash. 

#### Returns
float

#### `public float `[`GetWeightedDistance`](#classLSHForest_1acfb878c731daf8da6402e7cebc2b6ef1)`(const std::vector< uint32_t > & vec_a,const std::vector< uint32_t > & vec_b)` 

Get the distance between two weighted MinHashes.

#### Parameters
* `vec_a` A weighted MinHash. 

* `vec_b` A weighted MinHash. 

#### Returns
float

#### `public float `[`GetDistanceById`](#classLSHForest_1a49ad1fe0429121a8b572cde4df973d96)`(uint32_t a,uint32_t b)` 

Get the distance between two MinHashes.

#### Parameters
* `a` The id of an LSH forest entry. 

* `b` The id of an LSH forest entry. 

#### Returns
float

#### `public float `[`GetWeightedDistanceById`](#classLSHForest_1a72b5a201bc8c409c0c5742587988ba85)`(uint32_t a,uint32_t b)` 

Get the distance between two weighted MinHashes.

#### Parameters
* `a` The id of an LSH forest entry. 

* `b` The id of an LSH forest entry. 

#### Returns
float

#### `public void `[`Clear`](#classLSHForest_1aec34dc5185166dce9c22a4060ce3914c)`()` 

Remove all entries and the index from the LSH forest.

#### `public size_t `[`size`](#classLSHForest_1af0015fae65afd25bf67875d12dc0d663)`()` 

Get the number of entries.

#### Returns
size_t

# class `Minhash` 

An implementation of MinHash and weighted MinHash using SHA1.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public  `[`Minhash`](#classMinhash_1ad07fccee7e95ca368e5fabfb0ee76804)`(unsigned int d,unsigned int seed,unsigned int sample_size)` | Construct a new [Minhash](#classMinhash) object.
`public std::vector< uint32_t > `[`FromBinaryArray`](#classMinhash_1ae47faddc57a5d503257e6cf88dba2e08)`(std::vector< uint8_t > & vec)` | Create a MinHash from a binary array.
`public std::vector< std::vector< uint32_t > > `[`BatchFromBinaryArray`](#classMinhash_1a232f4fd24fcc853934599b666cbfc3be)`(std::vector< std::vector< uint8_t >> & vecs)` | Create MinHashes from a batch of binary arrays (parallelized).
`public std::vector< uint32_t > `[`FromSparseBinaryArray`](#classMinhash_1afe2cf6cc64b2e97ce89db4087febf30f)`(std::vector< uint32_t > & vec)` | Create a MinHash from a sparse binary array (values are the indices of 1s).
`public std::vector< std::vector< uint32_t > > `[`BatchFromSparseBinaryArray`](#classMinhash_1a8f711d80f1cb52f0599927667b123739)`(std::vector< std::vector< uint32_t >> & vecs)` | Create MinHashes from a vector of sparse binary arrays (values are the indices of 1s) (parallelized).
`public std::vector< uint32_t > `[`FromStringArray`](#classMinhash_1a7131b7dbefd40e0d24d7e37601519d62)`(std::vector< std::string > & vec)` | Create a MinHash from an array of strings.
`public std::vector< std::vector< uint32_t > > `[`BatchFromStringArray`](#classMinhash_1adc3ebe293e9999e49a30f9e9b9a8e318)`(std::vector< std::vector< std::string >> & vecs)` | Create MinHashes from a vector of string arrays (parallelized).
`public std::vector< uint32_t > `[`FromWeightArray`](#classMinhash_1a47a107b26e6fface715f5abfbb512484)`(std::vector< float > & vec)` | Create a MinHash from an array containing weights.
`public std::vector< std::vector< uint32_t > > `[`BatchFromWeightArray`](#classMinhash_1ad38f8679778e6291f5e006f76f104312)`(std::vector< std::vector< float >> & vecs)` | Create MinHashes from a vector of weight arrays (parallelized).
`public std::vector< uint8_t > `[`ExpandIntWeightArray`](#classMinhash_1a515153e559f07825616314f52bfd673f)`(std::vector< uint32_t > & vec,std::vector< uint32_t > & max_vec,uint32_t size)` | Expand a integer weight array into a binary array.
`public std::vector< std::vector< uint32_t > > `[`BatchFromIntWeightArray`](#classMinhash_1a6d468d2ee939351ffed9ea3ff8d82643)`(std::vector< std::vector< uint32_t >> & vecs)` | Create MinHashes from a expanded integer weight array (parallelized).
`public float `[`GetDistance`](#classMinhash_1a40dd607c20fa7c8059a3fb3d0c81ba0a)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` | Get the distance between two MinHashes.
`public float `[`GetWeightedDistance`](#classMinhash_1a8b2bd50a845fb4aa513464de10ed3e21)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` | Get the weighted distance between two MinHashes.
`public inline  `[`~Minhash`](#classMinhash_1a15a534f3f3e14c45ee20da1ed5039661)`()` | Destroy the [Minhash](#classMinhash) object.

## Members

#### `public  `[`Minhash`](#classMinhash_1ad07fccee7e95ca368e5fabfb0ee76804)`(unsigned int d,unsigned int seed,unsigned int sample_size)` 

Construct a new [Minhash](#classMinhash) object.

#### Parameters
* `d` The number of permutations used for hashing. 

* `seed` The seed for the random number generator. 

* `sample_size` The sample size when generating a weighted MinHash.

#### `public std::vector< uint32_t > `[`FromBinaryArray`](#classMinhash_1ae47faddc57a5d503257e6cf88dba2e08)`(std::vector< uint8_t > & vec)` 

Create a MinHash from a binary array.

#### Parameters
* `vec` A vector containing binary values. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromBinaryArray`](#classMinhash_1a232f4fd24fcc853934599b666cbfc3be)`(std::vector< std::vector< uint8_t >> & vecs)` 

Create MinHashes from a batch of binary arrays (parallelized).

#### Parameters
* `vecs` A vector of a vector containing binary values. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint32_t > `[`FromSparseBinaryArray`](#classMinhash_1afe2cf6cc64b2e97ce89db4087febf30f)`(std::vector< uint32_t > & vec)` 

Create a MinHash from a sparse binary array (values are the indices of 1s).

#### Parameters
* `vec` A vector of indices. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromSparseBinaryArray`](#classMinhash_1a8f711d80f1cb52f0599927667b123739)`(std::vector< std::vector< uint32_t >> & vecs)` 

Create MinHashes from a vector of sparse binary arrays (values are the indices of 1s) (parallelized).

#### Parameters
* `vecs` A vector of vectors of indices. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint32_t > `[`FromStringArray`](#classMinhash_1a7131b7dbefd40e0d24d7e37601519d62)`(std::vector< std::string > & vec)` 

Create a MinHash from an array of strings.

#### Parameters
* `vec` A vector of strings. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromStringArray`](#classMinhash_1adc3ebe293e9999e49a30f9e9b9a8e318)`(std::vector< std::vector< std::string >> & vecs)` 

Create MinHashes from a vector of string arrays (parallelized).

#### Parameters
* `vecs` A vector of string vectors. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint32_t > `[`FromWeightArray`](#classMinhash_1a47a107b26e6fface715f5abfbb512484)`(std::vector< float > & vec)` 

Create a MinHash from an array containing weights.

#### Parameters
* `vec` A vector of float weights. 

#### Returns
std::vector<uint32_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromWeightArray`](#classMinhash_1ad38f8679778e6291f5e006f76f104312)`(std::vector< std::vector< float >> & vecs)` 

Create MinHashes from a vector of weight arrays (parallelized).

#### Parameters
* `vecs` A vector of float vector weights. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public std::vector< uint8_t > `[`ExpandIntWeightArray`](#classMinhash_1a515153e559f07825616314f52bfd673f)`(std::vector< uint32_t > & vec,std::vector< uint32_t > & max_vec,uint32_t size)` 

Expand a integer weight array into a binary array.

#### Parameters
* `vec` A vector containing integer weights. 

* `max_vec` The maxima for all columns. 

* `size` The size of the expanded array. 

#### Returns
std::vector<uint8_t>

#### `public std::vector< std::vector< uint32_t > > `[`BatchFromIntWeightArray`](#classMinhash_1a6d468d2ee939351ffed9ea3ff8d82643)`(std::vector< std::vector< uint32_t >> & vecs)` 

Create MinHashes from a expanded integer weight array (parallelized).

#### Parameters
* `vecs` A vector of expanded integer weight vectors. 

#### Returns
std::vector<std::vector<uint32_t>>

#### `public float `[`GetDistance`](#classMinhash_1a40dd607c20fa7c8059a3fb3d0c81ba0a)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` 

Get the distance between two MinHashes.

#### Parameters
* `vec_a` A MinHash vector. 

* `vec_b` A MinHash vector. 

#### Returns
float

#### `public float `[`GetWeightedDistance`](#classMinhash_1a8b2bd50a845fb4aa513464de10ed3e21)`(std::vector< uint32_t > & vec_a,std::vector< uint32_t > & vec_b)` 

Get the weighted distance between two MinHashes.

#### Parameters
* `vec_a` A weighted MinHash vector. 

* `vec_b` A weighted MinHash vector. 

#### Returns
float

#### `public inline  `[`~Minhash`](#classMinhash_1a15a534f3f3e14c45ee20da1ed5039661)`()` 

Destroy the [Minhash](#classMinhash) object.

# class `Timer` 

A simple timer class used to check performance during development.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public inline  `[`Timer`](#classTimer_1a5f16e8da27d2a5a5242dead46de05d97)`()` | Construct a new [Timer](#classTimer) object and start the clock.
`public inline void `[`reset`](#classTimer_1a9020542d73357a4eef512eefaf57524b)`()` | Restart the clock.
`public inline double `[`elapsed`](#classTimer_1a6a89a613c2af9b0d1e5f7e4ba9e46c54)`() const` | Return the time elapsed since the timer was started or last reset.

## Members

#### `public inline  `[`Timer`](#classTimer_1a5f16e8da27d2a5a5242dead46de05d97)`()` 

Construct a new [Timer](#classTimer) object and start the clock.

#### `public inline void `[`reset`](#classTimer_1a9020542d73357a4eef512eefaf57524b)`()` 

Restart the clock.

#### `public inline double `[`elapsed`](#classTimer_1a6a89a613c2af9b0d1e5f7e4ba9e46c54)`() const` 

Return the time elapsed since the timer was started or last reset.

#### Returns
double

# struct `GraphProperties` 

The properties of a generated graph. An instance of this struct is returned from the layout functions.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public float `[`mst_weight`](#structGraphProperties_1a1924c22c71ea4f60ed7e31db768d5a9b) | 
`public uint32_t `[`n_connected_components`](#structGraphProperties_1a268ebbd3f3ef13fb56e44adfbd2583ca) | 
`public uint32_t `[`n_isolated_vertices`](#structGraphProperties_1ad3eec03efec4480e80c5510b896660ae) | 
`public std::vector< uint32_t > `[`degrees`](#structGraphProperties_1ae1515b9c7e47a9ed3034bb219871d1c1) | 
`public std::vector< std::vector< uint32_t > > `[`adjacency_list`](#structGraphProperties_1ab567d40199d8d3e6a4a96e8c145df635) | 

## Members

#### `public float `[`mst_weight`](#structGraphProperties_1a1924c22c71ea4f60ed7e31db768d5a9b) 

#### `public uint32_t `[`n_connected_components`](#structGraphProperties_1a268ebbd3f3ef13fb56e44adfbd2583ca) 

#### `public uint32_t `[`n_isolated_vertices`](#structGraphProperties_1ad3eec03efec4480e80c5510b896660ae) 

#### `public std::vector< uint32_t > `[`degrees`](#structGraphProperties_1ae1515b9c7e47a9ed3034bb219871d1c1) 

#### `public std::vector< std::vector< uint32_t > > `[`adjacency_list`](#structGraphProperties_1ab567d40199d8d3e6a4a96e8c145df635) 

# struct `LayoutConfiguration` 

A struct containing all the configuration options available for and applied to a layout.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public int `[`k`](#structLayoutConfiguration_1aeaf7404656862b423dfc63035f5e8bee) | 
`public int `[`kc`](#structLayoutConfiguration_1ac3e7def530d8c31537b8f9de6cd89fdc) | 
`public int `[`fme_iterations`](#structLayoutConfiguration_1afe033986b179ba0e419166b6fcb1dd35) | 
`public bool `[`fme_randomize`](#structLayoutConfiguration_1a44383f49d302581f3e8300fb4677663e) | 
`public int `[`fme_threads`](#structLayoutConfiguration_1ab356194e6ee8b8ad276fd196b540ac02) | 
`public int `[`fme_precision`](#structLayoutConfiguration_1a086e941fb5a53c068de377f740e330d7) | 
`public int `[`sl_repeats`](#structLayoutConfiguration_1ac48f682786c633d2e43c7ae83fa400e0) | 
`public int `[`sl_extra_scaling_steps`](#structLayoutConfiguration_1a5f21efd9d16cb86fde8c9f273ae8a05f) | 
`public double `[`sl_scaling_min`](#structLayoutConfiguration_1a267835cca2b8e0d694b7709009d7eaa5) | 
`public double `[`sl_scaling_max`](#structLayoutConfiguration_1a96f6f38497727a3844acab6be2ae021a) | 
`public `[`ScalingType`](#layout_8hh_1ae327227c361ab0e868a1f25017cb3ae2)` `[`sl_scaling_type`](#structLayoutConfiguration_1a7975f6fe7f0ca315af9c631a325b98af) | 
`public int `[`mmm_repeats`](#structLayoutConfiguration_1a73bfc92692894cafdfa063b0a45bc065) | 
`public `[`Placer`](#layout_8hh_1a93e50260439be3f5fe75b271c0ce2c96)` `[`placer`](#structLayoutConfiguration_1a139a9d88f1bcce6769b440f0f49130f0) | 
`public `[`Merger`](#layout_8hh_1a87e3986b1a6733e81a1c0b4bbd6aba18)` `[`merger`](#structLayoutConfiguration_1a70222497c34b2ffa597cd364d0a1d318) | 
`public double `[`merger_factor`](#structLayoutConfiguration_1a98f6187e2dc15b0f06bcbfef5562beae) | 
`public int `[`merger_adjustment`](#structLayoutConfiguration_1a382a084c8d4785151b9328221c4ba132) | 
`public float `[`node_size`](#structLayoutConfiguration_1a54a32d5173963abca63aae5bfa9d68e1) | 
`public inline  `[`LayoutConfiguration`](#structLayoutConfiguration_1a76742074edbb0cf0fad8d8c2d2f32be4)`()` | Construct a new Layout Configuration object.
`public inline std::string `[`ToString`](#structLayoutConfiguration_1a8be8ea09a3143cf9ba54a5069f0934d1)`() const` | Returns a string describing the set options.

## Members

#### `public int `[`k`](#structLayoutConfiguration_1aeaf7404656862b423dfc63035f5e8bee) 

#### `public int `[`kc`](#structLayoutConfiguration_1ac3e7def530d8c31537b8f9de6cd89fdc) 

#### `public int `[`fme_iterations`](#structLayoutConfiguration_1afe033986b179ba0e419166b6fcb1dd35) 

#### `public bool `[`fme_randomize`](#structLayoutConfiguration_1a44383f49d302581f3e8300fb4677663e) 

#### `public int `[`fme_threads`](#structLayoutConfiguration_1ab356194e6ee8b8ad276fd196b540ac02) 

#### `public int `[`fme_precision`](#structLayoutConfiguration_1a086e941fb5a53c068de377f740e330d7) 

#### `public int `[`sl_repeats`](#structLayoutConfiguration_1ac48f682786c633d2e43c7ae83fa400e0) 

#### `public int `[`sl_extra_scaling_steps`](#structLayoutConfiguration_1a5f21efd9d16cb86fde8c9f273ae8a05f) 

#### `public double `[`sl_scaling_min`](#structLayoutConfiguration_1a267835cca2b8e0d694b7709009d7eaa5) 

#### `public double `[`sl_scaling_max`](#structLayoutConfiguration_1a96f6f38497727a3844acab6be2ae021a) 

#### `public `[`ScalingType`](#layout_8hh_1ae327227c361ab0e868a1f25017cb3ae2)` `[`sl_scaling_type`](#structLayoutConfiguration_1a7975f6fe7f0ca315af9c631a325b98af) 

#### `public int `[`mmm_repeats`](#structLayoutConfiguration_1a73bfc92692894cafdfa063b0a45bc065) 

#### `public `[`Placer`](#layout_8hh_1a93e50260439be3f5fe75b271c0ce2c96)` `[`placer`](#structLayoutConfiguration_1a139a9d88f1bcce6769b440f0f49130f0) 

#### `public `[`Merger`](#layout_8hh_1a87e3986b1a6733e81a1c0b4bbd6aba18)` `[`merger`](#structLayoutConfiguration_1a70222497c34b2ffa597cd364d0a1d318) 

#### `public double `[`merger_factor`](#structLayoutConfiguration_1a98f6187e2dc15b0f06bcbfef5562beae) 

#### `public int `[`merger_adjustment`](#structLayoutConfiguration_1a382a084c8d4785151b9328221c4ba132) 

#### `public float `[`node_size`](#structLayoutConfiguration_1a54a32d5173963abca63aae5bfa9d68e1) 

#### `public inline  `[`LayoutConfiguration`](#structLayoutConfiguration_1a76742074edbb0cf0fad8d8c2d2f32be4)`()` 

Construct a new Layout Configuration object.

#### `public inline std::string `[`ToString`](#structLayoutConfiguration_1a8be8ea09a3143cf9ba54a5069f0934d1)`() const` 

Returns a string describing the set options.

#### Returns
std::string

# struct `MapKeyPointer` 

The pointer map used for pointing to the keys from the sorted hash map.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public iterator `[`it`](#structMapKeyPointer_1ab507d3d701fc46728b0ee1643e510edb) | 
`public inline  `[`MapKeyPointer`](#structMapKeyPointer_1a96ac4dd95d58f8e79be1e5c9c247b349)`(iterator i)` | 
`public inline  `[`MapKeyPointer`](#structMapKeyPointer_1a3a3f7e5e4b49c2d2ca1a3d40e43accf4)`()` | 
`public inline const std::vector< uint8_t > & `[`operator*`](#structMapKeyPointer_1a481cd84d48f25f2ff4bd102103a226d2)`() const` | 
`public inline const std::vector< uint8_t > * `[`operator->`](#structMapKeyPointer_1adf3870c8bdd7c9953fb07bb8708a7e73)`() const` | 
`typedef `[`iterator`](#structMapKeyPointer_1ae2724283aeda91caa9e63c5653be4a65) | 

## Members

#### `public iterator `[`it`](#structMapKeyPointer_1ab507d3d701fc46728b0ee1643e510edb) 

#### `public inline  `[`MapKeyPointer`](#structMapKeyPointer_1a96ac4dd95d58f8e79be1e5c9c247b349)`(iterator i)` 

#### `public inline  `[`MapKeyPointer`](#structMapKeyPointer_1a3a3f7e5e4b49c2d2ca1a3d40e43accf4)`()` 

#### `public inline const std::vector< uint8_t > & `[`operator*`](#structMapKeyPointer_1a481cd84d48f25f2ff4bd102103a226d2)`() const` 

#### `public inline const std::vector< uint8_t > * `[`operator->`](#structMapKeyPointer_1adf3870c8bdd7c9953fb07bb8708a7e73)`() const` 

#### `typedef `[`iterator`](#structMapKeyPointer_1ae2724283aeda91caa9e63c5653be4a65) 

# struct `SimpleHash` 

Hash struct used for the sparsepp sparse hash map.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public inline size_t `[`operator()`](#structSimpleHash_1abb89901eb1591020804460a3dccbde58)`(std::vector< uint8_t > vec) const` | 

## Members

#### `public inline size_t `[`operator()`](#structSimpleHash_1abb89901eb1591020804460a3dccbde58)`(std::vector< uint8_t > vec) const` 

Generated by [Moxygen](https://sourcey.com/moxygen)