// SamplingPreprocessing.h
#ifndef SAMPLING_PREPROCESSING_H
#define SAMPLING_PREPROCESSING_H

#include <vector>
/*struct TieRange {
    int start;  // inclusive
    int end;    // inclusive
    double value;

    bool contains(int idx) const {
        return idx >= start && idx <= end;
    }
};*/
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
    std::vector<double> uniform_denominator;
    //std::vector<TieRange> tie_ranges;
};

#endif // SAMPLING_PREPROCESSING_H