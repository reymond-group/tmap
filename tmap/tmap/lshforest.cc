/**
 * @file lshforest.cc
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief An LSH forest algorithm implementation.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#include "lshforest.hh"

tmap::LSHForest::LSHForest(unsigned int d,
                           unsigned int l,
                           bool store,
                           bool file_backed,
                           bool weighted)
  : d_(d)
  , l_(l)
  , size_(0)
  , store_(store)
  , file_backed_(file_backed)
  , weighted_(weighted)
  , hashtables_(l)
  , hashranges_(l)
  , sorted_hashtable_pointers_(l)
{
  if (l_ > d_)
    throw std::invalid_argument("l must be equal to or greater than d");

  k_ = (unsigned int)(d_ / l_);

  for (unsigned int i = 0; i < l_; i++) {
    hashtables_[i] = spp::sparse_hash_map<std::vector<uint8_t>,
                                          std::vector<uint32_t>,
                                          SimpleHash>();
    hashranges_[i] = std::make_tuple(i * k_, (i + 1) * k_);
  }
}

void
tmap::LSHForest::Add(std::vector<uint32_t>& vec)
{
  if (store_) {
    if (file_backed_) {
      std::ofstream fout("data.dat",
                         std::ios::ate | std::ios::app | std::ios::binary);
      fout.write((char*)&vec[0], vec.size() * sizeof(uint32_t));
      fout.close();
    } else
      data_.emplace_back(vec);
  }

  for (size_t i = 0; i < hashtables_.size(); i++) {
    std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                vec.begin() + std::get<1>(hashranges_[i]));

    hashtables_[i][Hash(Swap(range))].emplace_back(size_);
  }

  size_++;
  clean_ = false;
}

void
tmap::LSHForest::BatchAdd(std::vector<std::vector<uint32_t>>& vecs)
{
  size_t length = vecs.size();
  size_t data_length = size_;

  std::vector<uint32_t> keys(length);
  for (size_t i = 0; i < length; i++)
    keys[i] = data_length + i;

  if (store_) {
    if (file_backed_) {
      std::ofstream fout("data.dat",
                         std::ios::ate | std::ios::app | std::ios::binary);
      for (size_t i = 0; i < length; i++)
        fout.write((char*)&vecs[i][0], vecs[i].size() * sizeof(uint32_t));
      fout.close();
    } else {
      for (size_t i = 0; i < length; i++)
        data_.emplace_back(vecs[i]);
    }
  }

  int i, j;
#pragma omp parallel for private(j)
  for (i = 0; i < hashtables_.size(); i++) {
    for (j = 0; j < keys.size(); j++) {
      std::vector<uint32_t> range(vecs[j].begin() + std::get<0>(hashranges_[i]),
                                  vecs[j].begin() +
                                    std::get<1>(hashranges_[i]));

      hashtables_[i][Hash(Swap(range))].emplace_back(keys[j]);
    }
  }

  if (file_backed_) {
    // Free the memory
    std::vector<std::vector<uint32_t>>().swap(vecs);
  }

  size_ += length;
  clean_ = false;
}

void
tmap::LSHForest::Fit(std::vector<std::vector<uint32_t>>& vecs,
                     std::vector<uint32_t>& labels)
{
  if (labels_.size() != data_.size())
    throw std::runtime_error("LSHForest contains unlabelled entries.");
  
  if (vecs.size() != labels.size())
    throw std::runtime_error("The input sizes of vectors and labels has to match.");

  BatchAdd(vecs);
  for (size_t i = 0; i < labels.size(); i++)
    labels_.emplace_back(labels[i]);
  clean_ = false;
}

std::vector<uint32_t> 
tmap::LSHForest::Predict(std::vector<std::vector<uint32_t>>& vecs,
                         unsigned int k,
                         unsigned int kc,
                         bool weighted)
{
  std::vector<uint32_t> pred_labels(vecs.size());

  if (!weighted) {
#pragma omp parallel for
    for (size_t i = 0; i < vecs.size(); i++) {
      auto nn = QueryLinearScan(vecs[i], k, kc);

      std::sort(nn.begin(), nn.end(), [this](auto &left, auto &right) {
        return labels_[left.second] < labels_[right.second];
      });

      uint32_t max_element = labels_[nn[0].second];
      uint32_t max_count = 1;
      uint32_t count = 1;


      for (size_t j = 1; j < nn.size(); j++) {
        if (labels_[nn[j].second] == labels_[nn[j-1].second]) {
          count++;
          if (count > max_count) {
            max_count = count;
            max_element = labels_[nn[j].second];
          }
        } else {
          count = 1;
        }
      }

      pred_labels[i] = max_element;
    }
  } else {
#pragma omp parallel for
    for (size_t i = 0; i < vecs.size(); i++) {
      auto nn = QueryLinearScan(vecs[i], k, kc);

      std::sort(nn.begin(), nn.end(), [this](auto &left, auto &right) {
        return labels_[left.second] < labels_[right.second];
      });

      uint32_t max_element = labels_[nn[0].second];
      double max_count = 1.0 / (nn[0].first * nn[0].first);
      double count = max_count;


      for (size_t j = 1; j < nn.size(); j++) {
        if (labels_[nn[j].second] == labels_[nn[j-1].second]) {
          count += 1.0 / (nn[j].first * nn[j].first);
          if (count > max_count) {
            max_count = count;
            max_element = labels_[nn[j].second];
          }
        } else {
          count = 1.0 / (nn[j].first * nn[j].first);
        }
      }

      pred_labels[i] = max_element;
    }
  }

  return pred_labels;
}

void
tmap::LSHForest::Index()
{
#pragma omp parallel for
  for (int i = 0; i < hashtables_.size(); i++) {
    size_t j = 0;
    sorted_hashtable_pointers_[i].resize(hashtables_[i].size());
    for (auto it = hashtables_[i].begin(); it != hashtables_[i].end(); it++)
      sorted_hashtable_pointers_[i][j++] = MapKeyPointer(it);

    std::sort(sorted_hashtable_pointers_[i].begin(),
              sorted_hashtable_pointers_[i].end(),
              [](MapKeyPointer a, MapKeyPointer b) { return *a < *b; });
  }

  clean_ = true;
}

std::vector<uint32_t>
tmap::LSHForest::GetData(uint32_t id)
{
  if (file_backed_) {
    std::ifstream fin("data.dat", std::ios::in | std::ios::binary);

    size_t pos = id * d_ * sizeof(uint32_t);
    fin.seekg(pos);
    std::vector<uint32_t> result(d_);
    fin.read((char*)&result[0], result.size() * sizeof(uint32_t));
    fin.close();
    return result;
  } else {
    return data_[id];
  }
}

bool
tmap::LSHForest::IsClean()
{
  return clean_;
}

void
tmap::LSHForest::Store(const std::string& path)
{
  std::ofstream file(path, std::ios::binary);
  cereal::BinaryOutputArchive output(file);
  output(hashtables_, hashranges_, data_, labels_, store_, l_, d_, k_, clean_, size_);
  file.close();
}

void
tmap::LSHForest::Restore(const std::string& path)
{
  Clear();
  std::ifstream file(path, std::ios::binary);
  cereal::BinaryInputArchive input(file);
  input(hashtables_, hashranges_, data_, labels_, store_, l_, d_, k_, clean_, size_);
  file.close();

  sorted_hashtable_pointers_ = std::vector<std::vector<MapKeyPointer>>(l_);

  Index();
}

std::vector<uint32_t>
tmap::LSHForest::GetHash(uint32_t id)
{
  return GetData(id);
}

std::vector<std::pair<float, uint32_t>>
tmap::LSHForest::QueryLinearScan(const std::vector<uint32_t>& vec,
                                 unsigned int k,
                                 unsigned int kc)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  auto tmp = Query(vec, k * kc);
  return LinearScan(vec, tmp, k);
}

std::vector<std::pair<float, uint32_t>>
tmap::LSHForest::QueryLinearScanExclude(const std::vector<uint32_t>& vec,
                                        unsigned int k,
                                        std::vector<uint32_t>& exclude,
                                        unsigned int kc)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  auto tmp = QueryExclude(vec, exclude, k * kc);
  return LinearScan(vec, tmp, k);
}

std::vector<std::pair<float, uint32_t>>
tmap::LSHForest::QueryLinearScanById(uint32_t id,
                                     unsigned int k,
                                     unsigned int kc)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  return QueryLinearScan(GetData(id), k, kc);
}

std::vector<std::pair<float, uint32_t>>
tmap::LSHForest::QueryLinearScanExcludeById(uint32_t id,
                                            unsigned int k,
                                            std::vector<uint32_t>& exclude,
                                            unsigned int kc)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  return QueryLinearScanExclude(GetData(id), k, exclude, kc);
}

std::vector<std::pair<float, uint32_t>>
tmap::LSHForest::LinearScan(const std::vector<uint32_t>& vec,
                            std::vector<uint32_t>& indices,
                            unsigned int k)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  if (k == 0 || k > indices.size())
    k = indices.size();

  std::vector<std::pair<float, uint32_t>> result(indices.size());

  for (size_t i = 0; i < indices.size(); i++) {
    auto data = GetData(indices[i]);
    if (weighted_)
      result[i] =
        std::pair<float, uint32_t>(GetWeightedDistance(vec, data), indices[i]);
    else
      result[i] =
        std::pair<float, uint32_t>(GetDistance(vec, data), indices[i]);
  }

  std::sort(result.begin(), result.end());
  result.erase(result.begin() + k, result.end());

  return result;
}

// Does not always return k items. Is this expected?
std::vector<uint32_t>
tmap::LSHForest::Query(const std::vector<uint32_t>& vec, unsigned int k)
{
  std::set<uint32_t> results;

  for (int r = k_; r > 0; r--) {
    QueryInternal(vec, r, results, k);

    if (results.size() >= k)
      return std::vector<uint32_t>(results.begin(), results.end());
  }

  return std::vector<uint32_t>(results.begin(), results.end());
}

std::vector<uint32_t>
tmap::LSHForest::QueryExclude(const std::vector<uint32_t>& vec,
                              std::vector<uint32_t>& exclude,
                              unsigned int k)
{
  std::set<uint32_t> results;

  for (int r = k_; r > 0; r--) {
    QueryInternalExclude(vec, r, results, k, exclude);

    if (results.size() >= k)
      return std::vector<uint32_t>(results.begin(), results.end());
  }

  return std::vector<uint32_t>(results.begin(), results.end());
}

std::vector<uint32_t>
tmap::LSHForest::QueryById(uint32_t id, unsigned int k)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  return Query(GetData(id), k);
}

std::vector<uint32_t>
tmap::LSHForest::QueryExcludeById(uint32_t id,
                                  std::vector<uint32_t>& exclude,
                                  unsigned int k)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  return QueryExclude(GetData(id), exclude, k);
}

std::vector<std::vector<uint32_t>>
tmap::LSHForest::BatchQuery(const std::vector<std::vector<uint32_t>>& vecs,
                            unsigned int k)
{
  std::vector<std::vector<uint32_t>> results(vecs.size());

  for (unsigned int i = 0; i < vecs.size(); i++) {
    results[i] = Query(vecs[i], k);
  }

  return results;
}

std::vector<uint32_t>
tmap::LSHForest::GetAllNearestNeighbors(unsigned int k,
                                        unsigned int kc)
{
  std::vector<uint32_t> results(size_, 0);

#pragma omp parallel for
  for (int i = 0; i < size_; i++) {
    auto result = QueryLinearScanById(i, k, kc);
    if (result.size() > 1)
      results[i] = result[1].second;
  }

  return results;
}

void
tmap::LSHForest::GetKNNGraph(std::vector<uint32_t>& from,
                             std::vector<uint32_t>& to,
                             std::vector<float>& weight,
                             unsigned int k,
                             unsigned int kc)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");
  if (!clean_)
    throw std::runtime_error("LSHForest has not been (re)indexed");

  from.resize(size_ * k);
  to.resize(size_ * k);
  weight.resize(size_ * k);

  // Set weights that won't be set to negative one.
  for (size_t i = 0; i < weight.size(); i++)
    weight[i] = -1.0;

#pragma omp parallel for
  for (int i = 0; i < size_; i++) {
    auto result = QueryLinearScan(GetData(i), k, kc);

    for (size_t j = 0; j < result.size(); j++) {
      if (result[j].second == i)
        continue;

      from[k * i + j] = i;
      to[k * i + j] = result[j].second;
      weight[k * i + j] = result[j].first;
    }
  }
}

void
tmap::LSHForest::QueryInternal(const std::vector<uint32_t>& vec,
                               unsigned int r,
                               std::set<uint32_t>& results,
                               unsigned int k)
{
  std::vector<std::vector<uint8_t>> prefixes;

  for (size_t i = 0; i < hashranges_.size(); i++) {
    std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                vec.begin() + std::get<0>(hashranges_[i]) + r);

    prefixes.emplace_back(Hash(Swap(range)));
  }

  std::size_t prefix_size = prefixes[0].size();

  for (size_t i = 0; i < hashtables_.size(); i++) {
    auto& hashtable = hashtables_[i];
    auto& sorted_hashtable = sorted_hashtable_pointers_[i];
    auto& prefix = prefixes[i];

    unsigned int j = BinarySearch(sorted_hashtable.size(), [&](unsigned int x) {
      auto& sh = *sorted_hashtable[x];
      std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

      return range >= prefix;
    });

    if (j < sorted_hashtable.size()) {
      auto& sh = *sorted_hashtable[j];
      std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

      if (range != prefix)
        continue;
    }

    for (unsigned int l = j; l < sorted_hashtable.size(); l++) {
      auto& sh = *sorted_hashtable[l];
      std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

      if (range != prefix)
        break;

      auto& keys = hashtable[*sorted_hashtable[l]];

      for (size_t m = 0; m < keys.size(); m++) {
        results.insert(keys[m]);

        if (results.size() >= k)
          return;
      }
    }
  }
}

void
tmap::LSHForest::QueryInternalExclude(const std::vector<uint32_t>& vec,
                                      unsigned int r,
                                      std::set<uint32_t>& results,
                                      unsigned int k,
                                      std::vector<uint32_t>& exclude)
{
  std::vector<std::vector<uint8_t>> prefixes;

  for (size_t i = 0; i < hashranges_.size(); i++) {
    std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                vec.begin() + std::get<0>(hashranges_[i]) + r);

    prefixes.emplace_back(Hash(Swap(range)));
  }

  std::size_t prefix_size = prefixes[0].size();

  for (size_t i = 0; i < hashtables_.size(); i++) {
    auto& hashtable = hashtables_[i];
    auto& sorted_hashtable = sorted_hashtable_pointers_[i];
    auto& prefix = prefixes[i];

    unsigned int j = BinarySearch(sorted_hashtable.size(), [&](unsigned int x) {
      auto& sh = *sorted_hashtable[x];
      std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

      return range >= prefix;
    });

    if (j < sorted_hashtable.size()) {
      auto& sh = *sorted_hashtable[j];
      std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

      if (range != prefix)
        continue;
    }

    for (unsigned int l = j; l < sorted_hashtable.size(); l++) {
      auto& sh = *sorted_hashtable[l];
      std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

      if (range != prefix)
        break;

      auto& keys = hashtable[*sorted_hashtable[l]];

      for (size_t m = 0; m < keys.size(); m++) {
        if (std::find(exclude.begin(), exclude.end(), keys[m]) == exclude.end())
          results.insert(keys[m]);

        if (results.size() >= k)
          return;
      }
    }
  }
}

std::vector<uint8_t>
tmap::LSHForest::Hash(std::vector<uint32_t> vec)
{
  auto length = sizeof(vec[0]) * vec.size();

  std::vector<uint8_t> s(length);
  std::memcpy(s.data(), vec.data(), length);

  return s;
}

std::vector<uint8_t>
tmap::LSHForest::Hash(std::vector<std::vector<uint32_t>> vecs)
{
  // Linearize input matrix
  size_t stride = vecs[0].size();
  std::vector<uint32_t> lin(vecs.size() * stride);

  for (size_t i = 0; i < vecs.size(); i++)
    for (size_t j = 0; j < stride; j++)
      lin[i * stride + j] = vecs[i][j];

  return Hash(lin);
}

uint32_t
tmap::LSHForest::Swap(uint32_t i)
{
  return ((i >> 24) & 0xff) | ((i << 8) & 0xff0000) | ((i >> 8) & 0xff00) |
         ((i << 24) & 0xff000000);
}

std::vector<uint32_t>
tmap::LSHForest::Swap(std::vector<uint32_t> vec)
{
  std::vector<uint32_t> vec_out(vec.size());

  for (size_t i = 0; i < vec.size(); i++)
    vec_out[i] = Swap(vec[i]);

  return vec_out;
}

std::vector<std::vector<uint32_t>>
tmap::LSHForest::Swap(std::vector<std::vector<uint32_t>> vecs)
{
  std::vector<std::vector<uint32_t>> vecs_out(vecs.size());

  for (size_t i = 0; i < vecs.size(); i++) {
    vecs_out[i] = std::vector<uint32_t>(vecs[i].size());

    for (size_t j = 0; j < vecs[i].size(); j++)
      vecs_out[i][j] = Swap(vecs[i][j]);
  }

  return vecs_out;
}

unsigned int
tmap::LSHForest::BinarySearch(unsigned int n,
                              const std::function<bool(unsigned int)>& fn)
{
  unsigned int i = 0;
  unsigned int j = n;

  while (i < j) {
    unsigned int h = (unsigned int)(i + (j - i) / 2);

    if (!fn(h))
      i = h + 1;
    else
      j = h;
  }

  return i;
}

std::vector<float>
tmap::LSHForest::GetAllDistances(const std::vector<uint32_t>& vec)
{
  std::vector<float> dists(size_);

#pragma omp parallel for
  for (int i = 0; i < size_; i++) {
    dists[i] = GetDistance(vec, GetData(i));
  }

  return dists;
}

float
tmap::LSHForest::GetDistance(const std::vector<uint32_t>& vec_a,
                             const std::vector<uint32_t>& vec_b)
{
  float intersect = 0;

  for (unsigned int i = 0; i < d_; i++)
    if (vec_a[i] == vec_b[i])
      intersect++;

  return 1.0f - intersect / d_;
}

float
tmap::LSHForest::GetWeightedDistance(const std::vector<uint32_t>& vec_a,
                                     const std::vector<uint32_t>& vec_b)
{
  float intersect = 0.0f;
  for (unsigned int i = 0; i < d_ * 2; i += 2)
    if (vec_a[i] == vec_b[i] && vec_a[i + 1] == vec_b[i + 1])
      intersect++;

  return 1.0f - 2.0f * intersect / (float)d_;
}

float
tmap::LSHForest::GetDistanceById(uint32_t a, uint32_t b)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");

  return GetDistance(GetData(a), GetData(b));
}

float
tmap::LSHForest::GetWeightedDistanceById(uint32_t a, uint32_t b)
{
  if (!store_)
    throw std::runtime_error("LSHForest was not instantiated with store=true");

  return GetWeightedDistance(GetData(a), GetData(b));
}

size_t
tmap::LSHForest::size()
{
  return size_;
}

void
tmap::LSHForest::Clear()
{
  hashtables_ = std::vector<spp::sparse_hash_map<std::vector<uint8_t>,
                                                 std::vector<uint32_t>,
                                                 SimpleHash>>();
  hashranges_ = std::vector<std::tuple<uint32_t, uint32_t>>();
  data_ = std::vector<std::vector<uint32_t>>();
  labels_ = std::vector<uint32_t>();
  sorted_hashtable_pointers_ = std::vector<std::vector<MapKeyPointer>>();

  hashtables_.clear();
  hashtables_.shrink_to_fit();

  hashranges_.clear();
  hashranges_.shrink_to_fit();

  data_.clear();
  data_.shrink_to_fit();

  labels_.clear();
  labels_.shrink_to_fit();

  sorted_hashtable_pointers_.clear();
  sorted_hashtable_pointers_.shrink_to_fit();

  std::vector<spp::sparse_hash_map<std::vector<uint8_t>,
                                   std::vector<uint32_t>,
                                   SimpleHash>>()
    .swap(hashtables_);
  std::vector<std::tuple<uint32_t, uint32_t>>().swap(hashranges_);
  std::vector<std::vector<uint32_t>>().swap(data_);
  std::vector<uint32_t>().swap(labels_);
  std::vector<std::vector<MapKeyPointer>>().swap(sorted_hashtable_pointers_);
}

std::vector<std::vector<uint8_t>>
tmap::LSHForest::GetKeysFromHashtable(
  spp::sparse_hash_map<std::vector<uint8_t>, std::vector<uint32_t>, SimpleHash>
    hashtable)
{
  std::vector<std::vector<uint8_t>> keys;

  for (auto pair : hashtable)
    keys.emplace_back(pair.first);

  return keys;
}
