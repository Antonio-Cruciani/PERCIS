// NonUniformSampler.cpp
#include "NonUniformSampler.h"
#include "SamplingPreprocessing.h"
#include <algorithm>
#include <random>
#include <iostream>
#include <cmath>



std::vector<TieRange> compute_tie_ranges(const std::vector<double>& sorted_X) {
    std::vector<TieRange> tie_ranges;
    int n = sorted_X.size();
    int start = 0;

    for (int i = 1; i <= n; ++i) {
        if (i == n || sorted_X[i] != sorted_X[start]) {
            tie_ranges.push_back({start, i - 1, sorted_X[start]});
            start = i;
        }
    }
    return tie_ranges;
}

int find_tie_group(const std::vector<TieRange>& ranges, int idx) {
    int low = 0, high = ranges.size() - 1;
    while (low <= high) {
        int mid = (low + high) / 2;
        if (ranges[mid].contains(idx)) {
            return mid;
        } else if (idx < ranges[mid].start) {
            high = mid - 1;
        } else {
            low = mid + 1;
        }
    }
    throw std::runtime_error("Index not found in any tie range");
}


SamplingPreprocessing preprocessing(const std::vector<double>& X) {
    std::cout<<"Running Preprocessing"<<std::endl;
    std::random_device rd;
    std::mt19937 rng(rd()); 
    int n = X.size();
    std::vector<std::pair<double, int>> value_index_pairs;
    for (int i = 0; i < n; ++i) {
        value_index_pairs.emplace_back(X[i], i);
    }
    // This is a technical trick, in case of very skewed distribution of X we need to shuffle the indices before sorting them
    // This shuffles the ties and allows the binary search not to fall on the same element of the skewed group.
    // If you do not shuffle, in very skewed scenarios you might pick always the same s and obtain a bit of bias. 
    std::shuffle(value_index_pairs.begin(), value_index_pairs.end(), rng);

    std::sort(value_index_pairs.begin(), value_index_pairs.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    std::vector<double> sorted_X(n);
    std::vector<int> sorted_to_original(n);
    for (int i = 0; i < n; ++i) {
        sorted_X[i] = value_index_pairs[i].first;
        sorted_to_original[i] = value_index_pairs[i].second;
    }

    SamplingPreprocessing kernel;
    kernel.X = sorted_X;
    kernel.original_index_map = sorted_to_original;

    std::vector<double> w(n + 1, 0.0);
    std::vector<double> r(n + 1, 0.0);
    std::vector<double> c_vals(n, 0.0);
    double c_total = 0.0;

    for (int i = n - 1; i >= 0; --i) {
        w[i] = w[i + 1] + sorted_X[i];
        double c_i = (n - i) * sorted_X[i] - w[i];
        c_vals[i] = c_i;
        r[i] = r[i + 1] + c_i;
        c_total += c_i;
    }

    kernel.w = w;
    kernel.r = r;
    kernel.c_vals = c_vals;
    kernel.c_total = c_total;
    kernel.tie_ranges = compute_tie_ranges(kernel.X);
    std::cout<<"Preprocessing Completed"<<std::endl;
    return kernel;
}

std::pair<int, int> non_uniform_sampling_binary_search(const SamplingPreprocessing& kernel,std::mt19937& rng) {
    //std::vector<std::pair<int, int>> S;
    int n = kernel.X.size();

    
    std::uniform_real_distribution<> dis(0.0, 1.0);
  
    int a = 0, b = n - 1;
    while (a < b) {
        int d = ceil((a + b) / 2);
        double k = (kernel.c_total - kernel.r[d + 1]) / kernel.c_total;
        double u = dis(rng);
        if (u <= k) {
            b = d;
        } else {
            a = d + 1;
        }
    }
    int s = b;
    int tie_group_index = find_tie_group(kernel.tie_ranges, s);
    const auto& range = kernel.tie_ranges[tie_group_index];
    std::uniform_int_distribution<int> tie_dist(range.start, range.end);
    s = tie_dist(rng);
    //double x_s = kernel.X[s]; 
    a = s;
    b = n - 1;
    while (a < b) {
        int d = ceil((a + b) / 2);
        double numerator = (d - s + 1) * kernel.X[s] - kernel.w[s] + kernel.w[d + 1];
        double denominator = (n - s) * kernel.X[s] - kernel.w[s];
        double k = denominator == 0 ? 1.0 : numerator / denominator;
        double u = dis(rng);
        if (u <= k) {
            b = d;
        } else {
            a = d + 1;
        }
    }
    int t = b;

    //S.emplace_back(kernel.original_index_map[s], kernel.original_index_map[t]);
    

    return {kernel.original_index_map[s], kernel.original_index_map[t]};
}
