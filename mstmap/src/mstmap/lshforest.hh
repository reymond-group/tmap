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

#include <sparsepp/spp.h>

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
        LSHForest(unsigned int d = 128, unsigned int l = 8, bool store = false);
        ~LSHForest() {}
        void AddFromFile(std::string path);
        void Add(std::vector<uint32_t> &vec);
        void BatchAdd(std::vector<std::vector<uint32_t>> &vecs);
        void Index();
        bool IsClean();
        void Store(const std::string &path);
        void Restore(const std::string &path);
        std::vector<uint32_t> GetHash(uint32_t id);
        std::vector<std::tuple<uint32_t, uint32_t, float>> GetKNNGraph(unsigned int k, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScan(std::vector<uint32_t> &vec, unsigned int k, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScanExclude(std::vector<uint32_t> &vec, unsigned int k, std::vector<uint32_t> &exclude, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScanById(uint32_t id, unsigned int k, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> QueryLinearScanExcludeById(uint32_t id, unsigned int k, std::vector<uint32_t> &exclude, unsigned int kc = 10, bool weighted = false);
        std::vector<std::pair<float, uint32_t>> LinearScan(std::vector<uint32_t> &vec, std::vector<uint32_t> indices, unsigned int k = 0, bool weighted = false);
        std::vector<uint32_t> Query(std::vector<uint32_t> &vec, unsigned int k);
        std::vector<uint32_t> QueryExclude(std::vector<uint32_t> &vec, std::vector<uint32_t> &exclude, unsigned int k);
        std::vector<uint32_t> QueryById(uint32_t id, unsigned int k);
        std::vector<uint32_t> QueryExcludeById(uint32_t id, std::vector<uint32_t> &exclude, unsigned int k);
        std::vector<std::vector<uint32_t>> BatchQuery(std::vector<std::vector<uint32_t>> &vecs, unsigned int k);

        float GetDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b);
        float GetWeightedDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b);
        float GetDistanceById(uint32_t a, uint32_t b);
        float GetWeightedDistanceById(uint32_t a, uint32_t b);

        size_t size();
    private:
        unsigned int d_, l_, k_;
        bool clean_, store_;
        
        std::vector<std::vector<std::vector<uint8_t>>> hashtable_keys_;
        std::vector<std::vector<std::vector<uint32_t>>> hashtable_values_;
        std::vector<std::vector<size_t>> sort_maps_;

        std::vector<std::tuple<uint32_t, uint32_t>> hashranges_;
        std::vector<std::vector<std::vector<uint8_t>>> sorted_hashtables_;
        std::vector<std::vector<uint32_t>> data_;

        std::vector<spp::sparse_hash_map<std::vector<uint8_t>, std::vector<uint32_t>>> hashmaps_;

        uint32_t Swap(uint32_t i);
        std::vector<uint32_t> Swap(std::vector<uint32_t> vec);
        std::vector<std::vector<uint32_t>> Swap(std::vector<std::vector<uint32_t>> vecs);
        std::vector<uint8_t> Hash(std::vector<uint32_t> vec);
        std::vector<uint8_t> Hash(std::vector<std::vector<uint32_t>> vecs);
        std::vector<std::vector<uint8_t>> GetKeysFromHashtable(std::map<std::vector<uint8_t>, std::vector<uint32_t>> hashtable);
        
        unsigned int BinarySearch(unsigned int n, const std::function<bool(unsigned int)> &fn);
        void QueryInternal(std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k);
        void QueryInternalExclude(std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k, std::vector<uint32_t> &exclude);
};

#endif