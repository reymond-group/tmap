#include "euclideanhash.hh"

EuclideanHash::EuclideanHash(unsigned int d, unsigned int seed) : d_(d),
                                                                  perms_a_((uint32_t)0, d), perms_b_((uint32_t)0, d)
{
    // Permutations for the standard minhash
    std::mt19937 rand;
    rand.seed(seed);
    
    std::uniform_int_distribution<std::mt19937::result_type> dist_a(1, max_hash_);
    std::uniform_int_distribution<std::mt19937::result_type> dist_b(0, max_hash_);
    

    std::vector<uint32_t> perms_a_tmp(d, 0);
    std::vector<uint32_t> perms_b_tmp(d, 0);

    for (unsigned int i = 0; i < d; i++)
    {
        uint32_t a = dist_a(rand);
        uint32_t b = dist_b(rand);

        while (std::find(perms_a_tmp.begin(), perms_a_tmp.end(), a) != perms_a_tmp.end())
        {
            a = dist_a(rand);
        }

        while (std::find(perms_b_tmp.begin(), perms_b_tmp.end(), a) != perms_b_tmp.end())
        {
            b = dist_a(rand);
        }

        perms_a_tmp[i] = a;
        perms_b_tmp[i] = b;

        perms_a_[i] = a;
        perms_b_[i] = b;
    }

    // This is for the weighted minhash
    std::default_random_engine rand_gamma;
    rand_gamma.seed(seed);
    std::gamma_distribution<float> gamma_dist(2.0, 1.0);
    std::uniform_real_distribution<float> dist_beta(0.0f, 1.0f);
}

std::vector<uint32_t> EuclideanHash::FromWeightArray(std::vector<float> &vec)
{
    if (vec.size() != d_)
        throw std::invalid_argument("The length of the weighted input vector has to be the same length as the dimenstionality with which Minhash was initialized");

    std::vector<uint32_t> mh(d_);
    return mh;
}

std::vector<std::vector<uint32_t>> EuclideanHash::BatchFromWeightArray(std::vector<std::vector<float>> &vecs)
{
    std::vector<std::vector<uint32_t>> results(vecs.size());

    #pragma omp parallel for
    for (size_t i = 0; i < vecs.size(); i++)
        results[i] = FromWeightArray(vecs[i]);

    return results;
}

float EuclideanHash::GetDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b)
{
    float intersect = 0;
    size_t length = vec_a.size();

    for (unsigned int i = 0; i < length; i++)
    {
        if (vec_a[i] == vec_b[i])
            intersect++;
    }

    return 1.0f - intersect / length;
}

float EuclideanHash::GetWeightedDistance(std::vector<uint32_t> &vec_a, std::vector<uint32_t> &vec_b)
{
    float intersect = 0;
    size_t length = vec_a.size();

    for (unsigned int i = 0; i < length; i += 2)
    {
        if (vec_a[i] == vec_b[i] && vec_a[i + 1] == vec_b[i + 1])
            intersect++;
    }

    return 1.0f - 2.0f * intersect / length;
}
