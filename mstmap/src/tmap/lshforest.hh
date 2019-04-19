#ifndef LSHFOREST_H
#define LSHFOREST_H

#include <stdexcept>
#include <vector>
#include <map>
#include <tuple>
#include <set>
#include <stdint.h>
#include <cstring>
#include <algorithm>
#include <functional>
#include <iostream>
#include <fstream>
#include <chrono>

#include "cereal/archives/binary.hpp"
#include "cereal/types/vector.hpp"
#include "cereal/types/map.hpp"
#include "cereal/types/tuple.hpp"


#include "sparsepp/spp.h"
#include "layout.hh"

struct MyHash {
    size_t operator()(std::vector<uint8_t> vec) const 
    {   
        std::size_t seed = vec.size();
        for(auto& i : vec) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

struct MapKeyPointer
{
  typedef spp::sparse_hash_map<std::vector<uint8_t>, std::vector<uint32_t>>::iterator iterator;
  MapKeyPointer(iterator i) : it(i) {}
  MapKeyPointer() {}
  const std::vector<uint8_t>& operator*() const { return it->first; }
  const std::vector<uint8_t>* operator->() const { return &it->first; }
  iterator it;
};

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const { 
        return std::chrono::duration_cast<second_>
            (clock_::now() - beg_).count() * 1000; }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

class LSHForest
{
    public:
        LSHForest(unsigned int d = 128, unsigned int l = 8, bool store = false, bool file_backed = false);
        ~LSHForest() { }
        void AddFromFile(std::string path);
        void Add(std::vector<uint32_t> &vec);
        void BatchAdd(std::vector<std::vector<uint32_t>> &vecs);
        void Index();
        bool IsClean();
        void Store(const std::string &path);
        void Restore(const std::string &path);
        std::vector<uint32_t> GetHash(uint32_t id);
        void GetKNNGraph(std::vector<uint32_t> &from, std::vector<uint32_t> &to, std::vector<float> &weight, unsigned int k, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScan(const std::vector<uint32_t> &vec, unsigned int k, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScanExclude(const std::vector<uint32_t> &vec, unsigned int k, std::vector<uint32_t> &exclude, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScanById(uint32_t id, unsigned int k, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScanExcludeById(uint32_t id, unsigned int k, std::vector<uint32_t> &exclude, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> LinearScan(const std::vector<uint32_t> &vec, std::vector<uint32_t> &indices, unsigned int k = 0, bool weighted = false);
        void FastLinearScan(const std::vector<uint32_t> &vec, std::vector<uint32_t> &indices, std::vector<float> &weights, unsigned int k = 0, bool weighted = false);
        std::vector<uint32_t> Query(const std::vector<uint32_t> &vec, unsigned int k);
        std::vector<uint32_t> QueryExclude(const std::vector<uint32_t> &vec, std::vector<uint32_t> &exclude, unsigned int k);
        std::vector<uint32_t> QueryById(uint32_t id, unsigned int k);
        std::vector<uint32_t> QueryExcludeById(uint32_t id, std::vector<uint32_t> &exclude, unsigned int k);
        std::vector<std::vector<uint32_t>> BatchQuery(const std::vector<std::vector<uint32_t>> &vecs, unsigned int k);
        
        std::vector<uint32_t> GetData(uint32_t id);

        std::vector<float> GetAllDistances(const std::vector<uint32_t> &vec);
        float GetDistance(const std::vector<uint32_t> &vec_a, const std::vector<uint32_t> &vec_b);
        float GetWeightedDistance(const std::vector<uint32_t> &vec_a, const std::vector<uint32_t> &vec_b);
        float GetDistanceById(uint32_t a, uint32_t b);
        float GetWeightedDistanceById(uint32_t a, uint32_t b);

        std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
        GetLayout(LayoutConfiguration config = LayoutConfiguration(), bool create_mst = true, bool mem_dump = true);

        void Clear();
        size_t size();
    private:
        unsigned int d_, l_, k_;
        unsigned long size_;
        bool clean_, store_, file_backed_;
        // std::vector<spp::sparse_hash_map<const std::vector<uint8_t>, std::vector<uint32_t>, MyHash>> hashtables_;
        std::vector<spp::sparse_hash_map<std::vector<uint8_t>, std::vector<uint32_t>, MyHash>> hashtables_;
        std::vector<std::tuple<uint32_t, uint32_t>> hashranges_;
        std::vector<std::vector<MapKeyPointer>> sorted_hashtable_pointers_;
        std::vector<std::vector<uint32_t>> data_;

        uint32_t Swap(uint32_t i);
        std::vector<std::vector<uint8_t>> GetKeysFromHashtable(spp::sparse_hash_map<std::vector<uint8_t>, std::vector<uint32_t>, MyHash> hashtable);
        std::vector<uint32_t> Swap(std::vector<uint32_t> vec);
        std::vector<std::vector<uint32_t>> Swap(std::vector<std::vector<uint32_t>> vecs);
        std::vector<uint8_t> Hash(std::vector<uint32_t> vec);
        std::vector<uint8_t> Hash(std::vector<std::vector<uint32_t>> vecs);
        unsigned int BinarySearch(unsigned int n, const std::function<bool(unsigned int)> &fn);
        void QueryInternal(const std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k);
        void QueryInternalExclude(const std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k, std::vector<uint32_t> &exclude);
};

#endif