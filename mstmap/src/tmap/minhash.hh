#ifndef MINHASH_H
#define MINHASH_H

#include <stdexcept>
#include <vector>
#include <valarray>
#include <random>
#include <algorithm>
#include <stdint.h>
#include <iostream>
#include "TinySHA1.hh"

class Minhash
{
    public:
        Minhash(unsigned int d = 128, unsigned int seed = 42, unsigned int sample_size = 128);
        std::vector<uint32_t> FromBinaryArray(std::vector<uint8_t> &vec);
        std::vector<std::vector<uint32_t>> BatchFromBinaryArray(std::vector<std::vector<uint8_t>> &vecs);
        std::vector<uint32_t> FromSparseBinaryArray(std::vector<uint32_t> &vec);
        std::vector<std::vector<uint32_t>> BatchFromSparseBinaryArray(std::vector<std::vector<uint32_t>> &vecs);
        std::vector<uint32_t> FromStringArray(std::vector<std::string> &vec);
        std::vector<std::vector<uint32_t>> BatchFromStringArray(std::vector<std::vector<std::string>> &vecs);
        std::vector<uint32_t> FromWeightArray(std::vector<float> &vec);
        std::vector<std::vector<uint32_t>> BatchFromWeightArray(std::vector<std::vector<float>> &vecs);
        float GetDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b);
        float GetWeightedDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b);
        ~Minhash() {}

    private:
        unsigned int d_, sample_size_;
        uint64_t prime_ = 18446744073709551615UL;
        uint32_t max_hash_ = 4294967295;
        std::valarray<uint32_t> perms_a_;
        std::valarray<uint32_t> perms_b_;
        std::vector<std::valarray<float>> rs_;
        std::vector<std::valarray<float>> ln_cs_;
        std::vector<std::valarray<float>> betas_;

};

#endif