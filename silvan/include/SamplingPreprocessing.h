// SamplingPreprocessing.h
#ifndef SAMPLING_PREPROCESSING_H
#define SAMPLING_PREPROCESSING_H

#include <vector>

struct SamplingPreprocessing {
    std::vector<double> w;
    std::vector<double> c_vals;
    std::vector<double> r;
    double c_total;
    std::vector<double> X;
    std::vector<int> original_index_map;
    //Weights is used in case on linear non uniform sampling
    std::vector<double> weights;
    bool optimized = true;
};

#endif // SAMPLING_PREPROCESSING_H