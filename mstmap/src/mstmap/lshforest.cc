#include "lshforest.hh"

LSHForest::LSHForest(unsigned int d, unsigned int l) : d_(d), l_(l), hashtables_(l), hashranges_(l),
                                                       sorted_hashtables_(l)
{
    k_ = (unsigned int)(d_ / l_);

    for (unsigned int i = 0; i < l_; i++)
    {
        hashtables_[i] = std::map<std::vector<uint8_t>, std::vector<uint32_t>>();
        hashranges_[i] = std::make_tuple(i * k_, (i + 1) * k_);
        sorted_hashtables_[i] = std::vector<std::vector<uint8_t>>();
    }
}

void LSHForest::Add(uint32_t key, std::vector<uint32_t> &vec)
{
    std::vector<uint32_t> v(vec);
    data_.insert(std::pair<uint32_t, std::vector<uint32_t>>(key, v)); 

    for (size_t i = 0; i < hashranges_.size(); i++)
    {
        std::vector<uint32_t> range(vec.begin() + std::get<0>(hashranges_[i]),
                                    vec.begin() + std::get<1>(hashranges_[i]));
        keys_[key].emplace_back(Hash(Swap(range)));
    }

    for (size_t i = 0; i < hashtables_.size(); i++)
    {
        hashtables_[i][keys_[key][i]].emplace_back(key);
    }

    clean_ = false;
}

void LSHForest::Index()
{
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

std::set<uint32_t> LSHForest::Query(std::vector<uint32_t> &vec, unsigned int k)
{
    std::set<uint32_t> results;
    unsigned int r = k_;
    
    while (r > 0)
    {
        QueryInternal(vec, r, results, k);

        if (results.size() >= k)
            return results;

        r--;
    }
    return results;
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

        unsigned int j = BinarySearch(sorted_hashtable.size(), [&sorted_hashtable, &prefix, &prefix_size](unsigned int x) {
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
            bool skip = false;

            for (unsigned int m = 0; m < prefix_size; m++)
            {
                if (sorted_hashtable[l][m] != prefix[m])
                {
                    skip = true;
                    break;
                }
            }

            if (skip)
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

std::vector<uint8_t> LSHForest::Hash(std::vector<uint32_t> vec)
{
    auto length = sizeof(vec[0]) * vec.size();

    std::vector<uint8_t> s(length);
    std::memcpy(s.data(), vec.data(), length);
    
    return s;
}

uint32_t LSHForest::Swap(uint32_t i)
{
    return ((i >> 24) & 0xff) | ((i << 8) & 0xff0000) | ((i >> 8) & 0xff00) | ((i << 24) & 0xff000000);
}

std::vector<uint32_t> LSHForest::Swap(std::vector<uint32_t> vec)
{
    std::vector<uint32_t> vec_out(vec.size());

    for (size_t i = 0; i < vec.size(); i++)
    {
        vec_out[i] = Swap(vec[i]);
    }

    return vec_out;
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