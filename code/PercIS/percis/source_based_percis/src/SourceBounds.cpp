#include "SourceBounds.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>

namespace source_pc {
namespace {

using Real = long double;
constexpr Real numerical_slack = 512 * std::numeric_limits<Real>::epsilon();
constexpr Real relative_tolerance = 1e-8L;
constexpr std::size_t maximum_subdivisions = 262144;

void require_nonnegative(Real value, const char* description) {
    if (!std::isfinite(value) || value < 0)
        throw std::invalid_argument(description);
}

void require_delta(Real delta) {
    if (!std::isfinite(delta) || delta <= 0 || delta >= 1)
        throw std::invalid_argument("confidence delta must lie in (0,1)");
}

Real upward(Real value) {
    return std::nextafter(value, std::numeric_limits<Real>::infinity());
}

// Evaluate g(x) h(epsilon/g(x)) using log(x). The polynomial avoids
// catastrophic cancellation when epsilon/g is small; the logarithmic form
// avoids overflow and permits x smaller than the floating-point range.
Real bennett_denominator(Real log_x, Real epsilon) {
    const Real x = std::exp(log_x);
    const Real log_g = log_x + std::log1p(-x);
    const Real log_u = std::log(epsilon) - log_g;
    if (log_u <= std::log(0.125L)) {
        const Real u = std::exp(log_u);
        // h(u)/u^2 = sum_{j>=0} (-u)^j/((j+1)(j+2)).
        // The omitted tail is positive and smaller than 0.125^32/1122.
        Real polynomial = -1.0L / (32 * 33);
        for (int j = 30; j >= 0; --j) {
            const Real coefficient = (j % 2 ? -1.0L : 1.0L) /
                                     ((j + 1) * (j + 2));
            polynomial = coefficient + u * polynomial;
        }
        return epsilon * u * polynomial;
    }
    if (log_u > 0) {
        const Real inverse_u = std::exp(-log_u);
        const Real log_one_plus_u = log_u + std::log1p(inverse_u);
        return epsilon * ((1 + inverse_u) * log_one_plus_u - 1);
    }
    const Real u = std::exp(log_u);
    return epsilon * ((1 + u) * std::log1p(u) / u - 1);
}

struct Objective {
    Real log_a;
    Real epsilon;

    Real value_lower(Real log_x) const {
        const Real numerator = log_a - log_x;
        const Real denominator = bennett_denominator(log_x, epsilon);
        return std::max(0.0L, numerator / denominator *
                                 (1 - numerical_slack));
    }

    // On x in [left,right], log(A/x) decreases and the denominator
    // decreases, because d[g h(epsilon/g)]/dg = log(1+epsilon/g)-epsilon/g.
    Real interval_upper(Real log_left, Real log_right) const {
        const Real numerator = log_a - log_left;
        const Real denominator = bennett_denominator(log_right, epsilon);
        const Real result = numerator / denominator * (1 + numerical_slack);
        if (!std::isfinite(result) || result <= 0)
            throw std::overflow_error("PercIS Bennett sample bound overflow");
        return upward(result);
    }
};

struct Interval {
    Real left;
    Real right;
    Real upper;
    bool operator<(const Interval& other) const { return upper < other.upper; }
};

// This is the upstream epsilon-only endpoint search. It is retained as a
// subdivision hint, not as a domain restriction: the paper's theorem does
// not justify dropping the remainder of the variance/mean-capped domain.
Real legacy_split(Real epsilon) {
    Real low = std::max(0.0L, 0.5L - std::sqrt(epsilon / 3));
    Real high = 0.5L;
    while (high - low > 0.0001L) {
        const Real middle = (high + low) / 2;
        if (bennett_denominator(std::log(middle), epsilon) <=
            2 * epsilon * epsilon)
            high = middle;
        else
            low = middle;
    }
    return high;
}

std::uint64_t checked_ceiling(Real value) {
    const Real ceiling = std::ceil(value);
    // 2^64 is exactly representable even where long double equals double.
    if (!std::isfinite(ceiling) || ceiling >= std::ldexp(1.0L, 64))
        throw std::overflow_error("PercIS sample budget exceeds uint64_t");
    return static_cast<std::uint64_t>(std::max(1.0L, ceiling));
}

}  // namespace

long double percis_mean_upper(long double mean, std::uint64_t samples,
                              long double delta) {
    require_nonnegative(mean, "pilot mean must lie in [0,1]");
    if (mean > 1 || samples == 0)
        throw std::invalid_argument("pilot mean needs [0,1] observations and samples>0");
    require_delta(delta);
    const Real a = -std::log(delta) / static_cast<Real>(samples);
    const Real bound = mean + a + std::sqrt(a) * std::sqrt(a + 2 * mean);
    return std::min(1.0L, upward(bound * (1 + numerical_slack)));
}

long double percis_rho_upper(long double mean, long double unbiased_variance,
                             std::uint64_t samples, long double internal_bound,
                             long double delta) {
    require_nonnegative(mean, "pilot path-length mean must be finite and nonnegative");
    require_nonnegative(unbiased_variance, "pilot path-length variance must be finite and nonnegative");
    require_nonnegative(internal_bound, "internal-path bound must be finite and nonnegative");
    require_delta(delta);
    if (samples < 2 || mean > internal_bound)
        throw std::invalid_argument("rho pilot needs samples>=2 and mean<=internal bound");
    if (internal_bound == 0) {
        if (unbiased_variance != 0)
            throw std::invalid_argument("zero internal bound requires zero variance");
        return 0;
    }
    // Normalize first: D^2 or C^2 need not be representable.
    const Real m = static_cast<Real>(samples);
    const Real normalized_mean = mean / internal_bound;
    const Real a = -std::log(delta) / m;
    const Real bernstein = normalized_mean + a +
                           std::sqrt(a) * std::sqrt(a + 2 * normalized_mean);
    const Real log_two_delta = std::log(2.0L) - std::log(delta);
    const Real empirical = normalized_mean +
        (std::sqrt(unbiased_variance) / internal_bound) *
            std::sqrt(2 * log_two_delta / m) +
        7 * log_two_delta / (3 * (m - 1));
    // Each stochastic bound has failure probability delta. The caller
    // spends 2*delta on their minimum, as in the upstream pilot allocation.
    const Real normalized_bound = std::min({1.0L, bernstein, empirical});
    return std::min(internal_bound,
                    upward(internal_bound * normalized_bound *
                           (1 + numerical_slack)));
}

std::uint64_t percis_sample_budget(long double max_mean,
                                   long double max_second_or_var,
                                   long double rho,
                                   long double epsilon_normalized,
                                   long double delta_production) {
    require_nonnegative(max_mean, "maximum mean bound must lie in [0,1]");
    require_nonnegative(max_second_or_var, "variance/second-moment bound must lie in [0,1]");
    require_nonnegative(rho, "rho bound must be finite and nonnegative");
    require_delta(delta_production);
    if (max_mean > 1 || max_second_or_var > 1 ||
        !std::isfinite(epsilon_normalized) || epsilon_normalized <= 0 ||
        epsilon_normalized >= 1)
        throw std::invalid_argument("normalized epsilon must lie in (0,1), and mean/moment bounds in [0,1]");
    if (rho == 0 || max_second_or_var == 0 || max_mean == 0)
        return 1;

    const Real epsilon = epsilon_normalized;
    const Real variance = std::min(0.25L, max_second_or_var);
    // Stable equivalent of 1/2-sqrt(1/4-variance).
    const Real variance_endpoint = 2 * variance /
                                  (1 + std::sqrt(std::max(0.0L, 1 - 4 * variance)));
    // A maximum mean is also bounded by the sum of all means. This extra
    // cap uses the same rho event and keeps every logarithm positive.
    const Real endpoint = std::min({0.5L, max_mean, variance_endpoint, rho});
    const Real log_endpoint = std::log(endpoint);
    const Objective objective{std::log(2.0L) + std::log(rho) -
                              std::log(delta_production), epsilon};
    const Real inverse_epsilon = 1 / epsilon;
    (void)checked_ceiling(inverse_epsilon);
    Real lower = inverse_epsilon * (1 - numerical_slack);
    for (int j = 0; j < 4; ++j)
        lower = std::max(lower, objective.value_lower(log_endpoint - j * std::log(2.0L)));
    const Real split = std::min(endpoint, legacy_split(epsilon));
    lower = std::max(lower, objective.value_lower(std::log(split)));

    // Tail enclosure: for x<=a<epsilon/e,
    // J(x)>=epsilon[log(epsilon/x)-1], so sup_{0<x<=a} F(x)
    // <= max(1,log(A/a)/(log(epsilon/a)-1))/epsilon.
    // Selecting a in log space avoids underflow even for a long flat tail.
    const Real log_epsilon = std::log(epsilon);
    const Real target_ratio = (1 + relative_tolerance / 4) * lower * epsilon;
    const Real gap = std::max(relative_tolerance / 8, target_ratio - 1);
    const Real k = objective.log_a - log_epsilon + 1;
    const Real tail_b = std::max(1.0L, std::max(0.0L, k) / gap);
    const Real log_tail = std::min(log_endpoint - 1, log_epsilon - 1 - tail_b);
    const Real tail_denominator = log_epsilon - log_tail - 1;
    const Real tail_upper = upward(inverse_epsilon *
        std::max(1.0L, (objective.log_a - log_tail) / tail_denominator) *
        (1 + numerical_slack));

    std::priority_queue<Interval> intervals;
    auto add_interval = [&](Real left, Real right) {
        if (right > left)
            intervals.push({left, right, objective.interval_upper(left, right)});
    };
    const Real log_split = std::log(split);
    if (log_tail < log_split && log_split < log_endpoint) {
        add_interval(log_tail, log_split);
        add_interval(log_split, log_endpoint);
    } else {
        add_interval(log_tail, log_endpoint);
    }

    for (std::size_t step = 0; step < maximum_subdivisions; ++step) {
        const Interval top = intervals.top();
        if (top.upper <= lower * (1 + relative_tolerance))
            break;
        const Real middle = top.left + (top.right - top.left) / 2;
        if (!(middle > top.left && middle < top.right))
            break;
        intervals.pop();
        lower = std::max(lower, objective.value_lower(middle));
        add_interval(top.left, middle);
        add_interval(middle, top.right);
    }
    // Even if subdivision hits its deterministic work limit, the largest
    // remaining enclosure is valid; the limit never returns a sampled max.
    return checked_ceiling(std::max(tail_upper, intervals.top().upper));
}

}  // namespace source_pc
