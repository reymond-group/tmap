/**
 * @file lshforest.hh
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief An LSH forest algorithm implementation.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#ifndef LSHFOREST_H
#define LSHFOREST_H

#include <algorithm>
#include <chrono>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <set>
#include <stdexcept>
#include <stdint.h>
#include <tuple>
#include <vector>
#include <limits>

#include "cereal/archives/binary.hpp"
#include "cereal/types/map.hpp"
#include "cereal/types/tuple.hpp"
#include "cereal/types/vector.hpp"

#include "layout.hh"
#include "sparsepp/spp.h"

namespace tmap {
/**
 * @brief Hash struct used for the sparsepp sparse hash map.
 *
 */
struct SimpleHash
{
  size_t operator()(std::vector<uint8_t> vec) const
  {
    std::size_t seed = vec.size();
    for (auto& i : vec) {
      seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    return seed;
  }
};

/**
 * @brief The pointer map used for pointing to the keys from the sorted hash
 * map.
 *
 */
struct MapKeyPointer
{
  typedef spp::sparse_hash_map<std::vector<uint8_t>,
                               std::vector<uint32_t>>::iterator iterator;
  MapKeyPointer(iterator i)
    : it(i)
  {}
  MapKeyPointer() {}
  const std::vector<uint8_t>& operator*() const { return it->first; }
  const std::vector<uint8_t>* operator->() const { return &it->first; }
  iterator it;
};

/**
 * @brief A simple timer class used to check performance during development.
 *
 */
class Timer
{
public:
  /**
   * @brief Construct a new Timer object and start the clock.
   *
   */
  Timer()
    : beg_(clock_::now())
  {}

  /**
   * @brief Restart the clock.
   *
   */
  void reset() { beg_ = clock_::now(); }

  /**
   * @brief Return the time elapsed since the timer was started or last reset.
   *
   * @return double
   */
  double elapsed() const
  {
    return std::chrono::duration_cast<second_>(clock_::now() - beg_).count() *
           1000;
  }

private:
  typedef std::chrono::high_resolution_clock clock_;
  typedef std::chrono::duration<double, std::ratio<1>> second_;
  std::chrono::time_point<clock_> beg_;
};

/**
 * @brief Provides locality sensitive hashing forest functionalities.
 *
 */
class LSHForest
{
public:
  /**
   * @brief Construct a new LSHForest object.
   *
   * @param d The dimensionality of the MinHashes to be added to this LSHForest.
   * @param l The number of prefix trees used.
   * @param store Whether or not to store the data for later enhanced (using
   * parameter kc) retrievel.
   * @param file_backed Whether to store the data on disk rather than in RAM
   * (experimental).
   */
  LSHForest(unsigned int d = 128,
            unsigned int l = 8,
            bool store = true,
            bool file_backed = false,
            bool weighted = false);

  /**
   * @brief Destroy the LSHForest object.
   *
   */
  ~LSHForest() {}

  /**
   * @brief Add a MinHash to this LSHForest instance.
   *
   * @param vec A MinHash vector.
   */
  void Add(std::vector<uint32_t>& vec);

  /**
   * @brief Add Minhashes to this LSHForest (parallelized).
   *
   * @param vecs A vector containing MinHash vectors.
   */
  void BatchAdd(std::vector<std::vector<uint32_t>>& vecs);

  /**
   * @brief Add Minhashes with labels to this LSHForest (parallelized).
   *
   * @param vecs A vector containing MinHash vectors.
   * @param labels A vector containing labels.
   */
  void Fit(std::vector<std::vector<uint32_t>>& vecs, std::vector<uint32_t>& labels);

  /**
   * @brief Predict labels of Minhashes using the kNN algorithm (parallelized).
   *
   * @param vecs A vector containing MinHash vectors.
   * @param k The degree of the kNN algorithm.
   * @param kc The scalar by which k is multiplied before querying the LSH
   * @param weighted Whether distances are used as weights by the knn algorithm
   * 
   * @return std::vector<uint32_t> The predicted labels.
   */
  std::vector<uint32_t> Predict(std::vector<std::vector<uint32_t>>& vecs,
                                unsigned int k = 10,
                                unsigned int kc = 10,
                                bool weighted = false);

  /**
   * @brief Create the index (trees).
   *
   */
  void Index();

  /**
   * @brief Check whether the added MinHashes have been indexed.
   *
   * @return true
   * @return false
   */
  bool IsClean();

  /**
   * @brief Write / serialize the current LSH forest to the disk.
   *
   * @param path The location where the LSH forest should be stored on disk.
   */
  void Store(const std::string& path);

  /**
   * @brief Read / deserialize a LSH forest instance form the disk. The forest
   * is indexed automatically.
   *
   * @param path The location from where to load the LSH forest.
   */
  void Restore(const std::string& path);

  /**
   * @brief Get the MinHash of an entry at a given index. The index is defined
   * by order of insertion.
   *
   * @param id The index (order of insertion) of a entry.
   * @return std::vector<uint32_t> The MinHash associated with an index.
   */
  std::vector<uint32_t> GetHash(uint32_t id);

  /**
   * @brief Get the k-nearest neighbor graph of the data stored in this LSH
   * forest instance. It will be written to out parameters as an edge list.
   *
   * @param[out] from A vector to which the from vertices will be written.
   * @param[out] to A vector to which the to vertices will be written.
   * @param[out] weight A vector to which the float weights of the edges will be
   * written.
   * @param k The degree of the nearest neighbor graph.
   * @param kc The scalar by which k is multiplied before querying the LSH
   * forest. The results are then ordered decreasing based on linear-scan
   * distances and the top k results are picked to create the k-nearest neighbor
   * graph.
   */
  void GetKNNGraph(std::vector<uint32_t>& from,
                   std::vector<uint32_t>& to,
                   std::vector<float>& weight,
                   unsigned int k,
                   unsigned int kc = 10);

  /**
   * @brief Get the k-nearest neighbors of a query.
   *
   * @param vec The query MinHash.
   * @param k The number of k-nearest neighbors to return.
   * @param kc The scalar by which k is multiplied before querying the LSH
   * forest. The results are then ordered decreasing based on linear-scan
   * distances and the top k results returned.
   * @return std::vector<std::pair<float, uint32_t>> The distances and indices
   * of the k-nearest neighbors.
   */
  std::vector<std::pair<float, uint32_t>> QueryLinearScan(
    const std::vector<uint32_t>& vec,
    unsigned int k,
    unsigned int kc = 10);

  /**
   * @brief Get the k-nearest neighbors of a query except those defined in the
   * argument exclude.
   *
   * @param vec The query MinHash.
   * @param k The number of k-nearest neighbors to return.
   * @param exclude A list of indices to be excluded from the search
   * @param kc The scalar by which k is multiplied before querying the LSH
   * forest. The results are then ordered decreasing based on linear-scan
   * distances and the top k results returned.
   * @return std::vector<std::pair<float, uint32_t>>
   */
  std::vector<std::pair<float, uint32_t>> QueryLinearScanExclude(
    const std::vector<uint32_t>& vec,
    unsigned int k,
    std::vector<uint32_t>& exclude,
    unsigned int kc = 10);

  /**
   * @brief Get the k-nearest neighbors of an entry.
   *
   * @param id The id of the query entry.
   * @param k The number of k-nearest neighbors to return.
   * @param kc The scalar by which k is multiplied before querying the LSH
   * forest. The results are then ordered decreasing based on linear-scan
   * distances and the top k results returned.
   * @return std::vector<std::pair<float, uint32_t>> The distances and indices
   * of the k-nearest neighbors.
   */
  std::vector<std::pair<float, uint32_t>> QueryLinearScanById(
    uint32_t id,
    unsigned int k,
    unsigned int kc = 10);

  /**
   * @brief Get the k-nearest neighbors of an entry except those defined in the
   * argument exclude.
   *
   * @param id The id of the query entry.
   * @param k The number of k-nearest neighbors to return.
   * @param exclude A list of indices to be excluded from the search.
   * @param kc The scalar by which k is multiplied before querying the LSH
   * forest. The results are then ordered decreasing based on linear-scan
   * distances and the top k results returned.
   * @return std::vector<std::pair<float, uint32_t>> The distances and indices
   * of the k-nearest neighbors.
   */
  std::vector<std::pair<float, uint32_t>> QueryLinearScanExcludeById(
    uint32_t id,
    unsigned int k,
    std::vector<uint32_t>& exclude,
    unsigned int kc = 10);

  /**
   * @brief Get the k-nearest neighbors of a query using linear scan.
   *
   * @param vec The query MinHash.
   * @param indices A list of indices to on which to run the linear scan.
   * @param k The number of k-nearest neighbors to return.
   * @return std::vector<std::pair<float, uint32_t>> The distances and indices
   * of the k-nearest neighbors.
   */
  std::vector<std::pair<float, uint32_t>> LinearScan(
    const std::vector<uint32_t>& vec,
    std::vector<uint32_t>& indices,
    unsigned int k = 10);

  /**
   * @brief Query the LSH forest for k-nearest neighbors.
   *
   * @param vec The query MinHash.
   * @param k The number of nearest neighbors to search for.
   * @return std::vector<uint32_t> The indices of the k-nearest neighbors.
   */
  std::vector<uint32_t> Query(const std::vector<uint32_t>& vec, unsigned int k);

  /**
   * @brief Query the LSH forest for k-nearest neighbors. Exclude a list of
   * entries by ID.
   *
   * @param vec The query MinHash.
   * @param exclude A list of indices to be excluded from the search.
   * @param k The number of nearest neighbors to search for.
   * @return std::vector<uint32_t> The indices of the k-nearest neighbors.
   */
  std::vector<uint32_t> QueryExclude(const std::vector<uint32_t>& vec,
                                     std::vector<uint32_t>& exclude,
                                     unsigned int k);

  /**
   * @brief Query the LSH forest for k-nearest neighbors.
   *
   * @param id The id of the query entry.
   * @param k The number of nearest neighbors to search for.
   * @return std::vector<uint32_t> The indices of the k-nearest neighbors.
   */
  std::vector<uint32_t> QueryById(uint32_t id, unsigned int k);

  /**
   * @brief Query the LSH forest for k-nearest neighbors. Exclude a list of
   * entries by ID.
   *
   * @param id The id of the query entry.
   * @param exclude A list of indices to be excluded from the search.
   * @param k The number of nearest neighbors to search for.
   * @return std::vector<uint32_t> The indices of the k-nearest neighbors.
   */
  std::vector<uint32_t> QueryExcludeById(uint32_t id,
                                         std::vector<uint32_t>& exclude,
                                         unsigned int k);

  /**
   * @brief Query the LSH forest for k-nearest neighbors (parallelized).
   *
   * @param vecs A vector of MinHashes.
   * @param k The number of nearest neighbors to search for.
   * @return std::vector<std::vector<uint32_t>> A vector of the indices of the
   * k-nearest neighbors.
   */
  std::vector<std::vector<uint32_t>> BatchQuery(
    const std::vector<std::vector<uint32_t>>& vecs,
    unsigned int k);

  /**
   * @brief Get the k-nearest neighbors of all LSH forest entries.
   *
   * @param k The number of nearest neighbors to search for.
   * @param kc The scalar by which k is multiplied before querying the LSH
   * forest. The results are then ordered decreasing based on linear-scan
   * distances and the top k results returned.
   * @return std::vector<uint32_t> The IDs of the nearest neighbors of all LSH
   * forest entries.
   */
  std::vector<uint32_t> GetAllNearestNeighbors(unsigned int k,
                                               unsigned int kc = 10);

  /**
   * @brief Get the MinHash of an entry at a given index. The index is defined
   * by order of insertion. Alias for GetHash.
   *
   * @param id The index (order of insertion) of a entry.
   * @return std::vector<uint32_t> The MinHash associated with an index.
   */
  std::vector<uint32_t> GetData(uint32_t id);

  /**
   * @brief Get the distances of a MinHash to all entries in the LSH forest.
   *
   * @param vec The query MinHash.
   * @return std::vector<float> The distances form the input MinHash to all the
   * entries in the LSH forest.
   */
  std::vector<float> GetAllDistances(const std::vector<uint32_t>& vec);

  /**
   * @brief Get the distance between two MinHashes.
   *
   * @param vec_a A MinHash.
   * @param vec_b A MinHash.
   * @return float
   */
  float GetDistance(const std::vector<uint32_t>& vec_a,
                    const std::vector<uint32_t>& vec_b);

  /**
   * @brief Get the distance between two weighted MinHashes.
   *
   * @param vec_a A weighted MinHash.
   * @param vec_b A weighted MinHash.
   * @return float
   */
  float GetWeightedDistance(const std::vector<uint32_t>& vec_a,
                            const std::vector<uint32_t>& vec_b);

  /**
   * @brief Get the distance between two MinHashes.
   *
   * @param a The id of an LSH forest entry.
   * @param b The id of an LSH forest entry.
   * @return float
   */
  float GetDistanceById(uint32_t a, uint32_t b);

  /**
   * @brief Get the distance between two weighted MinHashes.
   *
   * @param a The id of an LSH forest entry.
   * @param b The id of an LSH forest entry.
   * @return float
   */
  float GetWeightedDistanceById(uint32_t a, uint32_t b);

  // std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>,
  // std::vector<uint32_t>, GraphProperties> GetLayout(LayoutConfiguration
  // config = LayoutConfiguration(), bool create_mst = true, bool mem_dump =
  // true);

  /**
   * @brief Remove all entries and the index from the LSH forest.
   *
   */
  void Clear();

  /**
   * @brief Get the number of entries.
   *
   * @return size_t
   */
  size_t size();

private:
  unsigned int d_, l_, k_;
  unsigned long size_;
  bool clean_, store_, file_backed_, weighted_;
  std::vector<spp::sparse_hash_map<std::vector<uint8_t>,
                                   std::vector<uint32_t>,
                                   SimpleHash>>
    hashtables_;
  std::vector<std::tuple<uint32_t, uint32_t>> hashranges_;
  std::vector<std::vector<MapKeyPointer>> sorted_hashtable_pointers_;
  std::vector<std::vector<uint32_t>> data_;
  std::vector<uint32_t> labels_;

  /**
   * @brief Byteswap a 32-bit unsigned integer.
   *
   * @param i A 32-bit unsigned integer.
   * @return uint32_t
   */
  uint32_t Swap(uint32_t i);

  /**
   * @brief Byteswap all 32-bit unsigned integers in a vector.
   *
   * @param vec A vector of 32-bit unsigned integers.
   * @return std::vector<uint32_t>
   */
  std::vector<uint32_t> Swap(std::vector<uint32_t> vec);

  /**
   * @brief Bytesap all 32-bit unsigned integers in a vector of vectors.
   *
   * @param vecs A vector of vectors of 32-bit unsigned integers.
   * @return std::vector<std::vector<uint32_t>>
   */
  std::vector<std::vector<uint32_t>> Swap(
    std::vector<std::vector<uint32_t>> vecs);

  /**
   * @brief Get the keys from a spp sparse hash map.
   *
   * @param hashtable A spp sparse hash map.
   * @return std::vector<std::vector<uint8_t>>
   */
  std::vector<std::vector<uint8_t>> GetKeysFromHashtable(
    spp::sparse_hash_map<std::vector<uint8_t>,
                         std::vector<uint32_t>,
                         SimpleHash> hashtable);

  /**
   * @brief Hash a vector of 32-bit unsigned integers.
   *
   * @param vec A vector of 32-bit unsigned integers.
   * @return std::vector<uint8_t>
   */
  std::vector<uint8_t> Hash(std::vector<uint32_t> vec);

  /**
   * @brief Hash a vector of vectors of 32-bit unsigned integers.
   *
   * @param vecs A vector of vectors of 32-bit unsigned integers.
   * @return std::vector<uint8_t>
   */
  std::vector<uint8_t> Hash(std::vector<std::vector<uint32_t>> vecs);

  /**
   * @brief Helper method to run a binary search.
   *
   * @param n The maximum.
   * @param fn The function that is evaluated in the search.
   * @return unsigned int
   */
  unsigned int BinarySearch(unsigned int n,
                            const std::function<bool(unsigned int)>& fn);

  /**
   * @brief The internal LSH forest query.
   *
   * @param vec The query MinHash.
   * @param r The tree depth / search depth.
   * @param results[out] A vector to which the results of the search are written
   * to.
   * @param k The number of nearest neighbors to search for.
   */
  void QueryInternal(const std::vector<uint32_t>& vec,
                     unsigned int r,
                     std::set<uint32_t>& results,
                     unsigned int k);

  /**
   * @brief The internal LSH forest query with exclusions.
   *
   * @param vec The query MinHash.
   * @param r The tree depth / search depth.
   * @param results
   * @param k A vector to which the results of the search are written to.
   * @param exclude A vector of entry IDs to be excluded from the search.
   */
  void QueryInternalExclude(const std::vector<uint32_t>& vec,
                            unsigned int r,
                            std::set<uint32_t>& results,
                            unsigned int k,
                            std::vector<uint32_t>& exclude);
};
}; // namespace tmap
#endif