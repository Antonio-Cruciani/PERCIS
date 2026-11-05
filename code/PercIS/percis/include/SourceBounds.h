#ifndef SOURCE_PC_BOUNDS_H
#define SOURCE_PC_BOUNDS_H

#include <cstdint>

namespace source_pc {

// PercIS pilot bounds. Mean observations have range [0,1]. The diameter
// observations have range [0,internal_bound], and variance is unbiased.
long double percis_mean_upper(long double mean, std::uint64_t samples,
                              long double delta);
long double percis_rho_upper(long double mean, long double unbiased_variance,
                             std::uint64_t samples, long double internal_bound,
                             long double delta);

// Normalized PercIS Bennett bound. max_second_or_var can be either an upper
// bound on the maximum second moment or on the maximum variance. The result
// covers the full theorem domain, including its limit at x=0. Numerical
// subdivision returns an upper enclosure, rather than a sampled maximum.
std::uint64_t percis_sample_budget(long double max_mean,
                                   long double max_second_or_var,
                                   long double rho,
                                   long double epsilon_normalized,
                                   long double delta_production);

}  // namespace source_pc

#endif
