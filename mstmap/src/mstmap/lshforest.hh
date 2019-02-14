#ifndef LSHFOREST_H
#define LSHFOREST_H

#include <vector>
#include <map>
#include <tuple>
#include <set>
#include <stdint.h>
#include <cstring>
#include <algorithm>
#include <functional>
#include <iostream>
#include <chrono>

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
        LSHForest(unsigned int d = 128, unsigned int l = 8);
        ~LSHForest() {}
        void Add(uint32_t key, std::vector<uint32_t> &vec);
        void Index();
        bool IsClean();
        std::set<uint32_t> Query(std::vector<uint32_t> &vec, unsigned int k);
        void QueryInternal(std::vector<uint32_t> &vec, unsigned int r, std::set<uint32_t> &results, unsigned int k);
        uint32_t Swap(uint32_t i);
        std::vector<uint32_t> Swap(std::vector<uint32_t> vec);
        std::vector<uint8_t> Hash(std::vector<uint32_t> vec);
        std::vector<std::vector<uint8_t>> GetKeysFromHashtable(std::map<std::vector<uint8_t>, std::vector<uint32_t>> hashtable);
        unsigned int BinarySearch(unsigned int n, const std::function<bool(unsigned int)> &fn);
    private:
        unsigned int d_, l_, k_;
        bool clean_;
        std::vector<std::map<std::vector<uint8_t>, std::vector<uint32_t>>> hashtables_;
        std::vector<std::tuple<uint32_t, uint32_t>> hashranges_;
        std::vector<std::vector<std::vector<uint8_t>>> sorted_hashtables_;
        std::map<uint32_t, std::vector<std::vector<uint8_t>>> keys_;
        std::map<uint32_t, std::vector<uint32_t>> data_;
};

#endif