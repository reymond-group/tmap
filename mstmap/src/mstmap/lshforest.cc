#include "lshforest.hh"

LSHForest::LSHForest(unsigned int d, unsigned int l, bool store) : d_(d), l_(l), store_(store), hashtables_(l),
                                                                   hashranges_(l), sorted_hashtables_(l)
{
    if (l_ > d_)
        throw std::invalid_argument("l must be equal to or greater than d");

    k_ = (unsigned int)(d_ / l_);

    for (unsigned int i = 0; i < l_; i++)
    {
        hashtables_[i] = std::map<std::vector<uint8_t>, std::vector<uint32_t>>();
        hashranges_[i] = std::make_tuple(i * k_, (i + 1) * k_);
        sorted_hashtables_[i] = std::vector<std::vector<uint8_t>>();
    }
}

void LSHForest::Add(std::vector<uint32_t> &vec)
{
    if (store_)
        data_.emplace_back(vec);

    for (size_t i = 0; i < hashtables_.size(); i++)
    {
        std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                    vec.begin() + std::get<1>(hashranges_[i]));
        
        hashtables_[i][Hash(Swap(range))].emplace_back(data_.size() - 1);
    }

    clean_ = false;
}

void LSHForest::BatchAdd(std::vector<std::vector<uint32_t>> &vecs)
{
    size_t length = vecs.size();
    size_t data_length = data_.size();
    std::vector<uint32_t> keys(length);

    if (store_) 
    {
        for (size_t i = 0; i < length; i++)
        {
            data_.emplace_back(vecs[i]);
            keys[i] = data_length + i;
        }
    }


    size_t i, j;
    #pragma omp parallel for private(j)
    for (i = 0; i < hashtables_.size(); i++)
    {
        for (j = 0; j < keys.size(); j++)
        {
            std::vector<uint32_t> range(vecs[j].begin() + std::get<0>(hashranges_[i]),
                                        vecs[j].begin() + std::get<1>(hashranges_[i]));
        
            hashtables_[i][Hash(Swap(range))].emplace_back(keys[j]);
        }
    }

    clean_ = false;
}

void LSHForest::Index()
{
    #pragma omp parallel for
    for (size_t i = 0; i < hashtables_.size(); i++)
    {
        sorted_hashtables_[i] = GetKeysFromHashtable(hashtables_[i]);
        std::sort(sorted_hashtables_[i].begin(), sorted_hashtables_[i].end());
    }

    clean_ = true;
}

bool LSHForest::IsClean()
{
    return clean_;
}

void LSHForest::Store(const std::string &path)
{
    std::ofstream file(path, std::ios::binary);
    cereal::BinaryOutputArchive output(file);
    output(hashtables_, hashranges_, sorted_hashtables_, data_, store_);
}

void LSHForest::Restore(const std::string &path)
{
    std::ifstream file(path, std::ios::binary);
    cereal::BinaryInputArchive input(file);
    input(hashtables_, hashranges_, sorted_hashtables_, data_, store_);
}

std::vector<uint32_t> LSHForest::GetHash(uint32_t id)
{
    return data_[id];
}

std::vector<std::pair<float, uint32_t>> LSHForest::QueryLinearScan(std::vector<uint32_t> &vec, unsigned int k, unsigned int kc, bool weighted)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");
    
    auto tmp = Query(vec, k * kc);
    return LinearScan(vec, tmp, k, weighted);
}

std::vector<std::pair<float, uint32_t>> LSHForest::QueryLinearScanExclude(std::vector<uint32_t> &vec, unsigned int k, std::vector<uint32_t> &exclude, unsigned int kc, bool weighted)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");
    
    auto tmp = QueryExclude(vec, exclude, k * kc);
    return LinearScan(vec, tmp, k, weighted);
}

std::vector<std::pair<float, uint32_t>> LSHForest::QueryLinearScanById(uint32_t id, unsigned int k, unsigned int kc, bool weighted)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");
    
    return QueryLinearScan(data_[id], k, kc, weighted);
}

std::vector<std::pair<float, uint32_t>> LSHForest::QueryLinearScanExcludeById(uint32_t id, unsigned int k, std::vector<uint32_t> &exclude, unsigned int kc, bool weighted)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");
    
    return QueryLinearScanExclude(data_[id], k, exclude, kc, weighted);
}

std::vector<std::pair<float, uint32_t>> LSHForest::LinearScan(std::vector<uint32_t> &vec, std::vector<uint32_t> indices, unsigned int k, bool weighted)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");

    if (k == 0 || k > indices.size())
        k = indices.size();

    std::vector<std::pair<float, uint32_t>> result(indices.size());

    for (size_t i = 0; i < indices.size(); i++)
    {
        if (weighted)
            result[i] = std::pair<float, uint32_t>(GetWeightedDistance(vec, data_[indices[i]]), indices[i]);
        else
            result[i] = std::pair<float, uint32_t>(GetDistance(vec, data_[indices[i]]), indices[i]);
    }

    std::sort(result.begin(), result.end());
    result.erase(result.begin() + k, result.end());

    return result;
}

// TODO: DOES NOT ALWAYS RETURN K!!!
std::vector<uint32_t> LSHForest::Query(std::vector<uint32_t> &vec, unsigned int k)
{
    std::set<uint32_t> results;

    for (int r = k_; r > 0; r--)
    {
        QueryInternal(vec, r, results, k);

        if (results.size() >= k)
            return std::vector<uint32_t>(results.begin(), results.end());

        r--;
    }

    return std::vector<uint32_t>(results.begin(), results.end());
}

std::vector<uint32_t> LSHForest::QueryExclude(std::vector<uint32_t> &vec, std::vector<uint32_t> &exclude, unsigned int k)
{
    std::set<uint32_t> results;

    for (int r = k_; r > 0; r--)
    {
        QueryInternalExclude(vec, r, results, k, exclude);

        if (results.size() >= k)
            return std::vector<uint32_t>(results.begin(), results.end());

        r--;
    }

    return std::vector<uint32_t>(results.begin(), results.end());
}

std::vector<uint32_t> LSHForest::QueryById(uint32_t id, unsigned int k)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");

    return Query(data_[id], k);
}

std::vector<uint32_t> LSHForest::QueryExcludeById(uint32_t id, std::vector<uint32_t> &exclude, unsigned int k)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");

    return QueryExclude(data_[id], exclude, k);
}

std::vector<std::vector<uint32_t>> LSHForest::BatchQuery(std::vector<std::vector<uint32_t>> &vecs, unsigned int k)
{
    std::vector<std::vector<uint32_t>> results(vecs.size());

    for (unsigned int i = 0; i < vecs.size(); i++)
    {
        results[i] = Query(vecs[i], k);
    }

    return results;
}

std::vector<std::tuple<uint32_t, uint32_t, float>> LSHForest::GetKNNGraph(unsigned int k, unsigned int kc, bool weighted)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");

    std::vector<std::tuple<uint32_t, uint32_t, float>> knn_graph(data_.size() * k);

    size_t j = 0;

    #pragma omp parallel for private(j)
    for(size_t i = 0; i < data_.size(); i++)
    {
        auto result = QueryLinearScan(data_[i], k, kc, weighted);

        for (j = 0; j < result.size(); j++)
        {
            knn_graph[k * i + j] = std::tuple<uint32_t, uint32_t, float>(i, result[j].second, result[j].first);
        } 
    }
    
    return knn_graph;
}

void LSHForest::QueryInternal(std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k)
{
    std::vector<std::vector<uint8_t>> prefixes;

    for (size_t i = 0; i < hashranges_.size(); i++)
    {
        std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                    vec.begin() + std::get<0>(hashranges_[i]) + r);

        prefixes.emplace_back(Hash(Swap(range)));
    }

    std::size_t prefix_size = prefixes[0].size();
    
    for (size_t i = 0; i < sorted_hashtables_.size(); i++)
    {
        auto &hashtable = hashtables_[i];
        auto &sorted_hashtable = sorted_hashtables_[i];
        auto &prefix = prefixes[i];

        unsigned int j = BinarySearch(sorted_hashtable.size(), [&](unsigned int x) {
            auto &sh = sorted_hashtable[x];
            std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

            return range >= prefix;
        });

        
        if (j < sorted_hashtable.size())
        {
            auto &sh = sorted_hashtable[j];
            std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);
            
            if (range != prefix) 
                continue;
        }


        for (unsigned int l = j; l < sorted_hashtable.size(); l++)
        {
            auto &sh = sorted_hashtable[l];
            std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);
            
            if (range != prefix) 
                break;

            auto &keys = hashtable[sorted_hashtable[l]];
            
            for (size_t m = 0; m < keys.size(); m++)
            {
                results.insert(keys[m]);

                if (results.size() >= k)
                    return;
            }
        }
    }
}

void LSHForest::QueryInternalExclude(std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k, std::vector<uint32_t> &exclude)
{
    std::vector<std::vector<uint8_t>> prefixes;

    for (size_t i = 0; i < hashranges_.size(); i++)
    {
        std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                    vec.begin() + std::get<0>(hashranges_[i]) + r);

        prefixes.emplace_back(Hash(Swap(range)));
    }

    std::size_t prefix_size = prefixes[0].size();
    
    for (size_t i = 0; i < sorted_hashtables_.size(); i++)
    {
        auto &hashtable = hashtables_[i];
        auto &sorted_hashtable = sorted_hashtables_[i];
        auto &prefix = prefixes[i];

        unsigned int j = BinarySearch(sorted_hashtable.size(), [&](unsigned int x) {
            auto &sh = sorted_hashtable[x];
            std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);

            return range >= prefix;
        });

        
        if (j < sorted_hashtable.size())
        {
            auto &sh = sorted_hashtable[j];
            std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);
            
            if (range != prefix) 
                continue;
        }


        for (unsigned int l = j; l < sorted_hashtable.size(); l++)
        {
            auto &sh = sorted_hashtable[l];
            std::vector<uint8_t> range(sh.begin(), sh.begin() + prefix_size);
            
            if (range != prefix) 
                break;

            auto &keys = hashtable[sorted_hashtable[l]];
            
            for (size_t m = 0; m < keys.size(); m++)
            {
                if (std::find(exclude.begin(), exclude.end(), keys[m]) == exclude.end())
                    results.insert(keys[m]);

                if (results.size() >= k)
                    return;
            }
        }
    }
}

std::vector<uint8_t> LSHForest::Hash(std::vector<uint32_t> vec)
{
    auto length = sizeof(vec[0]) * vec.size();

    std::vector<uint8_t> s(length);
    std::memcpy(s.data(), vec.data(), length);
    
    return s;
}

std::vector<uint8_t> LSHForest::Hash(std::vector<std::vector<uint32_t>> vecs)
{
    // Linearize input matrix
    size_t stride = vecs[0].size();
    std::vector<uint32_t> lin(vecs.size() * stride);

    for (size_t i = 0; i < vecs.size(); i++)
        for (size_t j = 0; j < stride; j++)
            lin[i * stride + j] = vecs[i][j];

    return Hash(lin);
}

uint32_t LSHForest::Swap(uint32_t i)
{
    return ((i >> 24) & 0xff) | ((i << 8) & 0xff0000) | ((i >> 8) & 0xff00) | ((i << 24) & 0xff000000);
}

std::vector<uint32_t> LSHForest::Swap(std::vector<uint32_t> vec)
{
    std::vector<uint32_t> vec_out(vec.size());

    for (size_t i = 0; i < vec.size(); i++)
        vec_out[i] = Swap(vec[i]);

    return vec_out;
}

std::vector<std::vector<uint32_t>> LSHForest::Swap(std::vector<std::vector<uint32_t>> vecs)
{
    std::vector<std::vector<uint32_t>> vecs_out(vecs.size());

    for (size_t i = 0; i < vecs.size(); i++)
    {
        vecs_out[i] = std::vector<uint32_t>(vecs[i].size());

        for (size_t j = 0; j < vecs[i].size(); j++)
            vecs_out[i][j] = Swap(vecs[i][j]);
    }

    return vecs_out;
}

std::vector<std::vector<uint8_t>> LSHForest::GetKeysFromHashtable(std::map<std::vector<uint8_t>,
                                                                           std::vector<uint32_t>>
                                                                           hashtable)
{
    std::vector<std::vector<uint8_t>> keys;

    for (auto pair : hashtable)
        keys.emplace_back(pair.first);

    return keys;
}

unsigned int LSHForest::BinarySearch(unsigned int n, const std::function<bool(unsigned int)> &fn)
{
    unsigned int i = 0;
    unsigned int j = n;

    while (i < j)
    {
        unsigned int h = (unsigned int)(i + (j - i) / 2);
        
        if (!fn(h))
            i = h + 1;
        else
            j = h;
    }

    return i;
}

float LSHForest::GetDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b)
{
    unsigned int intersect = 0;
    size_t length = vec_a.size();

    for (unsigned int i = 0; i < length; i++)
    {
        if (vec_a[i] == vec_b[i])
            intersect++;
    }

    return 1.0f - intersect / (float)length;
}

float LSHForest::GetWeightedDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b)
{
    unsigned int intersect = 0;
    size_t length = vec_a.size();

    for (unsigned int i = 0; i < length; i += 2)
    {
        if (vec_a[i] == vec_b[i] && vec_a[i + 1] == vec_b[i + 1])
            intersect++;
    }

    return 1.0f -  2.0f * intersect / (float)length;
}

float LSHForest::GetDistanceById(uint32_t a, uint32_t b)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");

    return GetDistance(data_[a], data_[b]);
}

float LSHForest::GetWeightedDistanceById(uint32_t a, uint32_t b)
{
    if (!store_)
        throw std::runtime_error("LSHForest was not instantiated with store=true");

    return GetWeightedDistance(data_[a], data_[b]);
}

size_t LSHForest::size()
{
    return data_.size();
}