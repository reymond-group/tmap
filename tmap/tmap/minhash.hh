/**
 * @file minhash.hh
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief A MinHash algorithm implementation.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#ifndef MINHASH_H
#define MINHASH_H

#include "fnv.hh"
#include <algorithm>
#include <iostream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <stdint.h>
#include <math.h>
#include <vector>

#include "omp.h"

namespace tmap {
/**
 * @brief An implementation of MinHash and weighted MinHash using SHA1.
 *
 */
class Minhash
{
public:
  /**
   * @brief Construct a new Minhash object.
   *
   * @param d The number of permutations used for hashing.
   * @param seed The seed for the random number generator.
   * @param sample_size The sample size when generating a weighted MinHash.
   */
  Minhash(unsigned int d = 128,
          unsigned int seed = 42,
          unsigned int sample_size = 128);

  /**
   * @brief Create a MinHash from a binary array.
   *
   * @param vec A vector containing binary values.
   * @return std::vector<uint32_t>
   */
  std::vector<uint32_t> FromBinaryArray(std::vector<uint8_t>& vec);

  /**
   * @brief Create MinHashes from a batch of binary arrays (parallelized).
   *
   * @param vecs A vector of a vector containing binary values.
   * @return std::vector<std::vector<uint32_t>>
   */
  std::vector<std::vector<uint32_t>> BatchFromBinaryArray(
    std::vector<std::vector<uint8_t>>& vecs);

  /**
   * @brief Create a MinHash from a sparse binary array (values are the indices
   * of 1s).
   *
   * @param vec A vector of indices.
   * @return std::vector<uint32_t>
   */
  std::vector<uint32_t> FromSparseBinaryArray(std::vector<uint32_t>& vec);

  /**
   * @brief Create MinHashes from a vector of sparse binary arrays (values are
   * the indices of 1s) (parallelized).
   *
   * @param vecs A vector of vectors of indices.
   * @return std::vector<std::vector<uint32_t>>
   */
  std::vector<std::vector<uint32_t>> BatchFromSparseBinaryArray(
    std::vector<std::vector<uint32_t>>& vecs);

  /**
   * @brief Create a MinHash from an array of strings.
   *
   * @param vec A vector of strings.
   * @return std::vector<uint32_t>
   */
  std::vector<uint32_t> FromStringArray(std::vector<std::string>& vec);

  /**
   * @brief Create MinHashes from a vector of string arrays (parallelized).
   *
   * @param vecs A vector of string vectors.
   * @return std::vector<std::vector<uint32_t>>
   */
  std::vector<std::vector<uint32_t>> BatchFromStringArray(
    std::vector<std::vector<std::string>>& vecs);

  /**
   * @brief Create a MinHash from an array containing weights.
   *
   * @param vec A vector of float weights.
   * @param method The method which to use to calculate the weighted
   * fingerprint. Options are "I2CWS" and "ICWS".
   * @return std::vector<uint32_t>
   */
  std::vector<uint32_t> FromWeightArray(std::vector<float>& vec,
                                        const std::string& method = "ICWS");

  /**
   * @brief Create MinHashes from a vector of weight arrays (parallelized).
   *
   * @param vecs A vector of float vector weights.
   * @return std::vector<std::vector<uint32_t>>
   */
  std::vector<std::vector<uint32_t>> BatchFromWeightArray(
    std::vector<std::vector<float>>& vecs,
    const std::string& method = "ICWS");

  /**
   * @brief Expand a integer weight array into a binary array.
   *
   * @param vec A vector containing integer weights.
   * @param max_vec The maxima for all columns.
   * @param min_vec The minima for all columns.
   * @param size The size of the expanded array.
   * @return std::vector<uint8_t>
   */
  std::vector<uint8_t> ExpandIntWeightArray(std::vector<uint32_t>& vec,
                                            std::vector<uint32_t>& min_vec,
                                            std::vector<uint32_t>& max_vec,
                                            uint32_t size);

  /**
   * @brief Create MinHashes from a expanded integer weight array
   * (parallelized).
   *
   * @param vecs A vector of expanded integer weight vectors.
   * @return std::vector<std::vector<uint32_t>>
   */
  std::vector<std::vector<uint32_t>> BatchFromIntWeightArray(
    std::vector<std::vector<uint32_t>>& vecs, uint8_t divide_by = 0);

  /**
   * @brief Get the distance between two MinHashes.
   *
   * @param vec_a A MinHash vector.
   * @param vec_b A MinHash vector.
   * @return float
   */
  float GetDistance(std::vector<uint32_t>& vec_a, std::vector<uint32_t>& vec_b);

  /**
   * @brief Get the weighted distance between two MinHashes.
   *
   * @param vec_a A weighted MinHash vector.
   * @param vec_b A weighted MinHash vector.
   * @return float
   */
  float GetWeightedDistance(std::vector<uint32_t>& vec_a,
                            std::vector<uint32_t>& vec_b);

  /**
   * @brief Destroy the Minhash object.
   *
   */
  ~Minhash() {}

private:
  unsigned int d_, sample_size_;
  uint64_t prime_ = 2305843009213693951UL;
  uint32_t max_hash_ = 4294967295;
  std::vector<uint32_t> perms_a_;
  std::vector<uint32_t> perms_b_;
  std::vector<std::vector<float>> rs_;
  std::vector<std::vector<float>> rs_2_;
  std::vector<std::vector<float>> ln_cs_;
  std::vector<std::vector<float>> cs_;
  std::vector<std::vector<float>> betas_;
  std::vector<std::vector<float>> betas_2_;
  std::hash<std::string> hasher;

  uint64_t fast_mod_long(const uint64_t input, const uint64_t ceil)
  {
    return input >= ceil ? input % ceil : input;
  }
};
}; // namespace tmap
#endif