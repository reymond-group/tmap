#ifndef EUCLIDEANHASH_H
#define EUCLIDEANHASH_H

#include <stdexcept>
#include <vector>
#include <valarray>
#include <random>
#include <algorithm>
#include <stdint.h>
#include <iostream>
#include "TinySHA1.hh"

class EuclideanHash
{
    public:
        EuclideanHash(unsigned int d = 128, unsigned int seed = 42);
        std::vector<uint32_t> FromWeightArray(std::vector<float> &vec);
        std::vector<std::vector<uint32_t>> BatchFromWeightArray(std::vector<std::vector<float>> &vecs);
        float GetDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b);
        float GetWeightedDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b);
        ~EuclideanHash() {}

    private:
        unsigned int d_;
        uint64_t prime_ = 18446744073709551615UL;
        uint32_t max_hash_ = 4294967295;
        std::valarray<uint32_t> perms_a_;
        std::valarray<uint32_t> perms_b_;
};

#endif