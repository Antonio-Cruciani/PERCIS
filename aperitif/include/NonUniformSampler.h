#ifndef NON_UNIFORM_SAMPLER_H
#define NON_UNIFORM_SAMPLER_H
#include <vector>
#include <utility>
#include "SamplingPreprocessing.h"
#include <random>

SamplingPreprocessing preprocessing(const std::vector<double>& X);
std::pair<int, int> non_uniform_sampling_binary_search(const SamplingPreprocessing& kernel,const int n,std::mt19937& rng);

#endif // NON_UNIFORM_SAMPLER_H
