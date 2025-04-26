// NonUniformSampler.cpp
#include "NonUniformSampler.h"
#include <algorithm>
#include <random>

SamplingPreprocessing preprocessing(const std::vector<double>& X) {
    int n = X.size();
    std::vector<std::pair<double, int>> value_index_pairs;
    for (int i = 0; i < n; ++i) {
        value_index_pairs.emplace_back(X[i], i);
    }

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
    return kernel;
}

std::pair<int, int> non_uniform_sampling_binary_search(const SamplingPreprocessing& kernel) {
    //std::vector<std::pair<int, int>> S;
    int n = kernel.X.size();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
  
    int a = 0, b = n - 1;
    while (a < b) {
        int d = (a + b) / 2;
        double k = (kernel.c_total - kernel.r[d + 1]) / kernel.c_total;
        double u = dis(gen);
        if (u <= k) {
            b = d;
        } else {
            a = d + 1;
        }
    }
    int s = b;

    a = s;
    b = n - 1;
    while (a < b) {
        int d = (a + b) / 2;
        double numerator = (d - s + 1) * kernel.X[s] - kernel.w[s] + kernel.w[d + 1];
        double denominator = (n - s) * kernel.X[s] - kernel.w[s];
        double k = denominator == 0 ? 1.0 : numerator / denominator;
        double u = dis(gen);
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
