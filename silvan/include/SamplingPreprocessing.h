// SamplingPreprocessing.h
#ifndef SAMPLING_PREPROCESSING_H
#define SAMPLING_PREPROCESSING_H

#include <vector>

struct SamplingPreprocessing {
    std::vector<double> w;
    std::vector<double> c;
    std::vector<double> r;
    double total_c;
};

#endif // SAMPLING_PREPROCESSING_H