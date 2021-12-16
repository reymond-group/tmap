/**
 * @file minhash.cc
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief A MinHash algorithm implementation.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#include "minhash.hh"

tmap::Minhash::Minhash(unsigned int d,
                       unsigned int seed,
                       unsigned int sample_size)
  : d_(d)
  , sample_size_(sample_size)
  , perms_a_(d, 0)
  , perms_b_(d, 0)
  , rs_(sample_size)
  , rs_2_(sample_size)
  , ln_cs_(sample_size)
  , cs_(sample_size)
  , betas_(sample_size)
  , betas_2_(sample_size)
{
  

  // Permutations for the standard minhash
  std::mt19937 rand;
  rand.seed(seed);

  std::uniform_int_distribution<std::mt19937::result_type> dist_a(1, max_hash_);
  std::uniform_int_distribution<std::mt19937::result_type> dist_b(0, max_hash_);

  std::vector<uint32_t> perms_a_tmp(d, 0);
  std::vector<uint32_t> perms_b_tmp(d, 0);

  for (unsigned int i = 0; i < d; i++) {
    uint32_t a = dist_a(rand);
    uint32_t b = dist_b(rand);

    while (std::find(perms_a_tmp.begin(), perms_a_tmp.end(), a) !=
           perms_a_tmp.end()) {
      a = dist_a(rand);
    }

    while (std::find(perms_b_tmp.begin(), perms_b_tmp.end(), b) !=
           perms_b_tmp.end()) {
      b = dist_a(rand);
    }

    perms_a_tmp[i] = a;
    perms_b_tmp[i] = b;

    perms_a_[i] = a;
    perms_b_[i] = b;
  }

  std::random_device rd;
  std::mt19937 gen(rd());

  // This is for the weighted minhash
  std::mt19937 rand_gamma;
  rand_gamma.seed(seed);
  std::mt19937 rand_gamma_2;
  rand_gamma_2.seed(seed * 2);
  std::mt19937 rand_gamma_3;
  rand_gamma_3.seed(seed * 4);
  std::gamma_distribution<float> gamma_dist(2.0, 1.0);
  std::gamma_distribution<float> gamma_dist_2(2.0, 1.0);
  std::gamma_distribution<float> gamma_dist_3(2.0, 1.0);

  std::uniform_real_distribution<float> dist_beta(0.0f, 1.0f);
  std::uniform_real_distribution<float> dist_beta_2(0.0f, 1.0f);

  for (unsigned int i = 0; i < sample_size_; i++) {
    rs_[i] = std::vector<float>(d_, 0.0f);
    rs_2_[i] = std::vector<float>(d_, 0.0f);
    ln_cs_[i] = std::vector<float>(d_, 0.0f);
    cs_[i] = std::vector<float>(d_, 0.0f);
    betas_[i] = std::vector<float>(d_, 0.0f);
    betas_2_[i] = std::vector<float>(d_, 0.0f);

    for (unsigned int j = 0; j < d_; j++) {
      rs_[i][j] = gamma_dist(rand_gamma);
      rs_2_[i][j] = gamma_dist_2(rand_gamma_2);
      ln_cs_[i][j] = std::log(gamma_dist_3(rand_gamma_3));
      cs_[i][j] = gamma_dist_3(rand_gamma_3);
      betas_[i][j] = dist_beta(rand);
      betas_2_[i][j] = dist_beta(rand);
    }
  }
}

std::vector<uint32_t>
tmap::Minhash::FromBinaryArray(std::vector<uint8_t>& vec)
{
  std::vector<uint32_t> mh(d_, max_hash_);
  std::vector<uint32_t> tmp(d_);

  for (uint32_t i = 0; i < vec.size(); i++) {
    if (vec[i] == 0)
      continue;

    for (size_t j = 0; j < d_; j++)
      tmp[j] =
        (fast_mod_long((perms_a_[j] * i + perms_b_[j]), prime_)) & max_hash_;

    for (size_t j = 0; j < d_; j++)
      mh[j] = std::min(tmp[j], mh[j]);
  }

  return std::vector<uint32_t>(std::begin(mh), std::end(mh));
}

std::vector<std::vector<uint32_t>>
tmap::Minhash::BatchFromBinaryArray(std::vector<std::vector<uint8_t>>& vecs)
{
  std::vector<std::vector<uint32_t>> results(vecs.size());

#pragma omp parallel for
  for (int i = 0; i < vecs.size(); i++)
    results[i] = FromBinaryArray(vecs[i]);

  return results;
}

std::vector<uint32_t>
tmap::Minhash::FromSparseBinaryArray(std::vector<uint32_t>& vec)
{
  std::vector<uint32_t> mh(d_, max_hash_);
  std::vector<uint32_t> tmp(d_);

  for (uint32_t i = 0; i < vec.size(); i++) {
    for (size_t j = 0; j < d_; j++)
      tmp[j] = (fast_mod_long((perms_a_[j] * vec[i] + perms_b_[j]), prime_)) &
               max_hash_;

    for (size_t j = 0; j < d_; j++)
      mh[j] = std::min(tmp[j], mh[j]);
  }

  return std::vector<uint32_t>(std::begin(mh), std::end(mh));
}

std::vector<std::vector<uint32_t>>
tmap::Minhash::BatchFromSparseBinaryArray(
  std::vector<std::vector<uint32_t>>& vecs)
{
  std::vector<std::vector<uint32_t>> results(vecs.size());

#pragma omp parallel for
  for (int i = 0; i < vecs.size(); i++)
    results[i] = FromSparseBinaryArray(vecs[i]);

  return results;
}

std::vector<uint32_t>
tmap::Minhash::FromStringArray(std::vector<std::string>& vec)
{
  std::vector<uint32_t> mh(d_, max_hash_);
  std::vector<uint32_t> tmp(d_);

  for (uint32_t i = 0; i < vec.size(); i++) {
    for (size_t j = 0; j < d_; j++)
      tmp[j] = (fast_mod_long((perms_a_[j] * FNV::fnv1a(vec[i]) + perms_b_[j]),
                              prime_)) &
               max_hash_;

    for (size_t j = 0; j < d_; j++)
      mh[j] = std::min(tmp[j], mh[j]);
  }

  return std::vector<uint32_t>(std::begin(mh), std::end(mh));
}

std::vector<std::vector<uint32_t>>
tmap::Minhash::BatchFromStringArray(std::vector<std::vector<std::string>>& vecs)
{
  std::vector<std::vector<uint32_t>> results(vecs.size());

#pragma omp parallel for
  for (int i = 0; i < vecs.size(); i++)
    results[i] = FromStringArray(vecs[i]);

  return results;
}

std::vector<uint8_t>
tmap::Minhash::ExpandIntWeightArray(std::vector<uint32_t>& vec,
                                    std::vector<uint32_t>& min_vec,
                                    std::vector<uint32_t>& max_vec,
                                    uint32_t size)
{
  std::vector<uint8_t> vec_expanded(size, 0);
  size_t index = 0;

  for (size_t i = 0; i < vec.size(); i++) {
    // LOL (I might leave this in as a monument to why it's a bad idea to code
    // at 2am on a Saturday)
    for (size_t j = min_vec[i]; j < max_vec[i] - vec[i]; j++)
      index++;

    for (size_t j = 0; j < vec[i]; j++)
      vec_expanded[index++] = 1;
  }

  return vec_expanded;
}

std::vector<std::vector<uint32_t>>
tmap::Minhash::BatchFromIntWeightArray(std::vector<std::vector<uint32_t>>& vecs, 
                                       uint8_t divide_by)
{
  size_t size = vecs[0].size();

  // Get the maxima to expand the array by
  std::vector<uint32_t> max_vec(size, 0);
  std::vector<uint32_t> min_vec(size, UINT32_MAX);
  
  for (size_t i = 0; i < vecs.size(); i++)
    for (size_t j = 0; j < size; j++) {
      max_vec[j] = std::max(vecs[i][j], max_vec[j]);
      min_vec[j] = std::min(vecs[i][j], min_vec[j]);
    }

  uint32_t extended_size = std::accumulate(max_vec.begin(), max_vec.end(), 0);

  std::vector<std::vector<uint8_t>> tmp_results(vecs.size());
  std::vector<std::vector<uint32_t>> results(vecs.size());

  bool is_divide_by_pow_two = std::ceil(std::log2(divide_by)) == floor(std::log2(divide_by));
  if (is_divide_by_pow_two)
    divide_by = std::log2(divide_by); // will be -inf when divide_by == 0

#pragma omp parallel for
  for (int i = 0; i < vecs.size(); i++) {
    if (divide_by > 0) {
      if (is_divide_by_pow_two) {
        for(size_t j = 0; j < vecs[i].size(); j++)
          vecs[i][j] = vecs[i][j] >> divide_by;
      } else {
        for(size_t j = 0; j < vecs[i].size(); j++)
          vecs[i][j] /= divide_by;
      }
    }
    tmp_results[i] = ExpandIntWeightArray(vecs[i], min_vec, max_vec, extended_size);
  }

#pragma omp parallel for
  for (int i = 0; i < vecs.size(); i++)
    results[i] = FromBinaryArray(tmp_results[i]);

  return results;
}

std::vector<uint32_t>
tmap::Minhash::FromWeightArray(std::vector<float>& vec,
                               const std::string& method)
{
  std::vector<uint32_t> mh(sample_size_ * 2);

  if (method == "ICWS") {
    // ICWS
    for (unsigned int s = 0; s < sample_size_; s++) {
      std::vector<float> vec_tmp(vec);
      std::vector<float> ln_a(d_, std::numeric_limits<float>::max());
      std::vector<float> rs = rs_[s];
      std::vector<float> betas = betas_[s];
      std::vector<float> ln_cs = ln_cs_[s];

      for (size_t i = 0; i < vec_tmp.size(); i++) {
        if (vec_tmp[i] <= 0) continue;
        vec_tmp[i] = std::floor((vec_tmp[i] / rs[i]) + betas[i]);
        float ln_y = (vec_tmp[i] - betas[i]) * rs[i];
        ln_a[i] = ln_cs[i] - ln_y - rs[i];
      }

      uint32_t k_star = std::min_element(ln_a.begin(), ln_a.end()) - ln_a.begin();

      mh[2 * s] = k_star;
      mh[2 * s + 1] = (uint32_t)vec_tmp[k_star];
    }
  } else {
    // I2CWS
    for (unsigned int s = 0; s < sample_size_; s++) {
      std::vector<float> vec_tmp(d_, 0);
      std::vector<float> a(d_, std::numeric_limits<float>::max());
      std::vector<float> rs = rs_[s];
      std::vector<float> rs_2 = rs_2_[s];
      std::vector<float> betas = betas_[s];
      std::vector<float> betas_2 = betas_2_[s];
      std::vector<float> cs = cs_[s];

      for (size_t i = 0; i < vec_tmp.size(); i++) {
        if (vec[i] <= 0) continue;
        vec_tmp[i] = std::floor((log(vec[i]) / rs_2[i]) + betas_2[i]);
        float z = std::exp(rs_2[i] * (vec_tmp[i] - betas_2[i] + 1));
        a[i] = cs[i] / z;
      }

      uint32_t k_star = std::min_element(a.begin(), a.end()) - a.begin();
      uint32_t t_k = std::floor((log(vec[k_star]) / rs[k_star]) + betas[k_star]);

      mh[2 * s] = k_star;
      mh[2 * s + 1] = std::exp(rs[k_star] * (t_k - betas[k_star]));
    }
  }

  return mh;
}

std::vector<std::vector<uint32_t>>
tmap::Minhash::BatchFromWeightArray(std::vector<std::vector<float>>& vecs,
                                    const std::string& method)
{
  std::vector<std::vector<uint32_t>> results(vecs.size());

#pragma omp parallel for
  for (int i = 0; i < vecs.size(); i++)
    results[i] = FromWeightArray(vecs[i], method);

  return results;
}

float
tmap::Minhash::GetDistance(std::vector<uint32_t>& vec_a,
                           std::vector<uint32_t>& vec_b)
{
  float intersect = 0;
  size_t length = vec_a.size();

  for (unsigned int i = 0; i < length; i++)
    if (vec_a[i] == vec_b[i])
      intersect++;

  return 1.0f - intersect / length;
}

float
tmap::Minhash::GetWeightedDistance(std::vector<uint32_t>& vec_a,
                                   std::vector<uint32_t>& vec_b)
{
  float intersect = 0.0f;
  size_t length = vec_a.size();

  for (unsigned int i = 0; i < length; i += 2)
    if (vec_a[i] == vec_b[i] && vec_a[i + 1] == vec_b[i + 1])
      intersect += 1.0f;

  return 1.0f - ((2.0f * intersect) / (float)length);
}
