#include "SourceOnly.h"
#include "SourceBounds.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace source_pc {
namespace {

struct Sum {
    long double value = 0;
    long double residual = 0;
    void add(long double x) {
        const long double next = value + x;
        if (std::fabs(value) >= std::fabs(x)) residual += (value - next) + x;
        else residual += (x - next) + value;
        value = next;
    }
    long double get() const { return value + residual; }
};

std::uint64_t mixed_seed(std::uint64_t x) {
    x += UINT64_C(0x9e3779b97f4a7c15);
    x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
    return x ^ (x >> 31);
}

void add_count(std::uint64_t& total, std::uint64_t increment) {
    if (increment > std::numeric_limits<std::uint64_t>::max() - total)
        throw std::overflow_error("Source-only sampling counter overflow");
    total += increment;
}

struct Moments {
    std::vector<Sum> sums;
    std::vector<Sum> squares;
    Sum lengths;
    Sum squared_lengths;
    std::uint64_t bags = 0;
    std::uint64_t paths = 0;
    std::uint64_t unreachable = 0;
    explicit Moments(std::size_t n) : sums(n), squares(n) {}
};

struct Worker {
    Sampler sampler;
    Moments moments;
    Worker(const GraphData& graph, const Model& model, std::uint64_t seed,
           std::uint64_t bag_size, long double alpha)
        : sampler(graph, model, seed, bag_size, alpha), moments(graph.size()) {}
};

void collect(Worker& worker, const Model& model, std::uint64_t count,
             const std::atomic<bool>& failed) {
    const long double tolerance = 32 * std::numeric_limits<long double>::epsilon();
    auto& moments = worker.moments;
    for (std::uint64_t i = 0; i < count; ++i) {
        if (failed.load(std::memory_order_relaxed)) return;
        const Bag bag = worker.sampler.sample();
        for (const auto& entry : bag.scores) {
            if (entry.first >= model.correction.size())
                throw std::runtime_error("Sampler returned an invalid vertex");
            // The SAME global range is used for all vertices in PercIS's
            // unit-range sample-size calculation. Normalize epsilon likewise.
            const long double z = entry.second / model.max_correction;
            if (!std::isfinite(z) || z < -tolerance || z > 1 + tolerance)
                throw std::runtime_error("Sampler returned a bag score outside its range");
            const long double value = std::max(0.0L, std::min(1.0L, z));
            moments.sums[entry.first].add(value);
            moments.squares[entry.first].add(value * value);
        }
        const long double length = bag.internal_vertices;
        moments.lengths.add(length);
        moments.squared_lengths.add(length * length);
        add_count(moments.bags, 1);
        add_count(moments.paths, bag.sampled_paths);
        if (!bag.reachable) add_count(moments.unreachable, 1);
        // Empty/unreachable bags stay in all moment denominators, including
        // the zero-length bin used to estimate rho in the pilot.
    }
}

Moments phase(const GraphData& graph, const Model& model, const RunOptions& options,
              int threads, std::uint64_t count, bool pilot, long double alpha) {
    std::vector<std::unique_ptr<Worker>> workers;
    workers.reserve(static_cast<std::size_t>(threads));
    // Production streams agree with a standalone -g run. Pilot streams are
    // separate; restarting the pilot's seed would reuse its samples.
    const auto base = pilot ? mixed_seed(options.seed ^ UINT64_C(0xd1b54a32d192ed03))
                            : options.seed;
    for (int t = 0; t < threads; ++t) {
        const auto seed = t == 0 ? base : mixed_seed(base + static_cast<std::uint64_t>(t));
        workers.emplace_back(new Worker(graph, model, seed, options.bag_size, alpha));
    }
    const auto thread_count = static_cast<std::uint64_t>(threads);
    const auto quotient = count / thread_count;
    const auto remainder = count % thread_count;
    std::atomic<bool> failed(false);
    std::exception_ptr exception;
    std::mutex exception_mutex;
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(static, 1)
#endif
    for (int t = 0; t < threads; ++t) {
        try {
            collect(*workers[static_cast<std::size_t>(t)], model,
                    quotient + (static_cast<std::uint64_t>(t) < remainder ? 1 : 0), failed);
        } catch (...) {
            std::lock_guard<std::mutex> lock(exception_mutex);
            if (!exception) exception = std::current_exception();
            failed.store(true, std::memory_order_relaxed);
        }
    }
    if (exception) std::rethrow_exception(exception);
    Moments combined(graph.size());
    for (const auto& worker : workers) {
        const auto& local = worker->moments;
        add_count(combined.bags, local.bags);
        add_count(combined.paths, local.paths);
        add_count(combined.unreachable, local.unreachable);
        combined.lengths.add(local.lengths.get());
        combined.squared_lengths.add(local.squared_lengths.get());
        for (std::size_t v = 0; v < graph.size(); ++v) {
            combined.sums[v].add(local.sums[v].get());
            combined.squares[v].add(local.squares[v].get());
        }
    }
    if (combined.bags != count)
        throw std::logic_error("Source-only sampler did not execute its exact bag budget");
    return combined;
}

std::uint64_t pilot_budget(long double normalized_epsilon, long double delta,
                           std::uint32_t edge_bound) {
    // Same pilot-size heuristic as Probabilistic::run, in normalized units.
    const long double log_confidence = -std::log(delta);
    const long double count = std::max({1000.0L, log_confidence / normalized_epsilon,
                                      2 * edge_bound * log_confidence});
    if (!std::isfinite(count) || count >=
        static_cast<long double>(std::numeric_limits<std::uint64_t>::max()))
        throw std::overflow_error("PercIS pilot sample size exceeds the counter range");
    return static_cast<std::uint64_t>(count);
}

long double sample_variance(const Sum& sum, const Sum& squares, std::uint64_t count) {
    if (count < 2) return std::numeric_limits<long double>::quiet_NaN();
    const long double m = static_cast<long double>(count);
    return std::max(0.0L, std::fma(-sum.get(), sum.get() / m, squares.get()) / (m - 1));
}

} // namespace

Model::Model(const GraphData& graph) {
    const std::size_t n = graph.size();
    if (n < 3) throw std::invalid_argument("Source-only centrality requires at least three vertices");
    if (graph.states.size() != n)
        throw std::invalid_argument("The number of states must equal the number of vertices");
    if (n > static_cast<std::size_t>(std::numeric_limits<Vertex>::max()))
        throw std::invalid_argument("Too many vertices for the source-only vertex type");

    std::size_t positive = 0;
    for (long double state : graph.states) {
        if (!std::isfinite(state) || state < 0 || state > 1)
            throw std::invalid_argument("Percolation states must be finite and in [0,1]");
        if (state > 0) ++positive;
    }
    if (positive < 2)
        throw std::invalid_argument("Source-only centrality for every vertex requires at least two positive states");

    std::vector<long double> prefix(n + 1, 0);
    std::vector<long double> suffix(n + 1, 0);
    Sum accumulated;
    for (std::size_t v = 0; v < n; ++v) {
        accumulated.add(graph.states[v]);
        prefix[v + 1] = accumulated.get();
    }
    accumulated = Sum();
    for (std::size_t v = n; v-- > 0;) {
        accumulated.add(graph.states[v]);
        suffix[v] = accumulated.get();
    }
    total_state = prefix[n];
    if (!std::isfinite(total_state) || !(total_state > 0))
        throw std::invalid_argument("Total percolation state is not finite and positive");

    correction.resize(n);
    cumulative_states.resize(n);
    const long double endpoint_factor =
        static_cast<long double>(n - 1) / static_cast<long double>(n - 2);
    for (std::size_t v = 0; v < n; ++v) {
        // Do not form total_state - x_v: a dominant state would erase tiny
        // positive states by cancellation. Prefix/suffix sums exclude v.
        Sum excluded;
        excluded.add(prefix[v]);
        excluded.add(suffix[v + 1]);
        const long double remaining = excluded.get();
        if (!std::isfinite(remaining) || !(remaining > 0))
            throw std::invalid_argument("Every vertex must have positive source mass outside itself");
        correction[v] = endpoint_factor *
            (1 + graph.states[v] / remaining);
        if (!std::isfinite(correction[v]) || !(correction[v] > 0))
            throw std::overflow_error("Source-only importance correction is not representable");
        max_correction = std::max(max_correction, correction[v]);
        // Unnormalized long-double prefix masses; Sampler draws in [0,S).
        cumulative_states[v] = prefix[v + 1];
    }
    cumulative_states.back() = total_state;
}

RunResult estimate(const GraphData& graph, const Model& model,
                   const RunOptions& options) {
    const std::size_t n = graph.size();
    if (model.correction.size() != n || !(model.total_state > 0) ||
        !std::isfinite(model.max_correction) || !(model.max_correction > 0))
        throw std::invalid_argument("Source-only model does not match the graph");
    if (options.threads < 0) throw std::invalid_argument("Thread count cannot be negative");
    int threads = options.threads;
#ifdef _OPENMP
    if (threads == 0) threads = omp_get_max_threads();
#else
    if (threads == 0) threads = 1;
    if (threads != 1)
        throw std::invalid_argument("This build has no OpenMP support; use --threads 1");
#endif
    if (!std::isfinite(options.alpha))
        throw std::invalid_argument("Shortest-path sampling alpha must be finite");
    if (!std::isfinite(options.epsilon) || !(options.epsilon > 0 && options.epsilon < 1))
        throw std::invalid_argument("Epsilon must lie strictly between zero and one");
    if (!std::isfinite(options.delta) || !(options.delta > 0 && options.delta < 1))
        throw std::invalid_argument("Delta must lie strictly between zero and one");
    if (options.pilot_samples == 1)
        throw std::invalid_argument("The pilot must contain at least two bags");
    if (options.samples && options.pilot_samples)
        throw std::invalid_argument("A fixed -g budget skips the pilot; omit --pilot-samples");

    RunResult result;
    result.threads = threads;
    result.estimates.assign(n, 0);
    result.sample_variances.assign(n, 0);
    result.production_alpha = options.alpha;
    result.planned_bags = options.samples;
    if (options.samples == 0) {
        // PercIS: (1) pilot, (2) choose ONE sample size, (3) independent
        // production run. No production checkpoints or early stopping.
        result.edge_diameter_bound = percis_edge_diameter_bound(graph);
        const auto internal_bound = result.edge_diameter_bound > 0
                                      ? result.edge_diameter_bound - 1 : 0;
        if (internal_bound == 0) {
            if (!options.quiet)
                std::cerr << "percis: no shortest path can have internal vertices\n";
            return result;
        }
        const long double normalized_epsilon = options.epsilon / model.max_correction;
        if (!(normalized_epsilon > 0))
            throw std::overflow_error("Normalized epsilon is not representable");
        const auto tau = options.pilot_samples ? options.pilot_samples
            : pilot_budget(normalized_epsilon, options.delta, result.edge_diameter_bound);
        if (!options.quiet)
            std::cerr << "percis: phase=pilot bags=" << tau
                      << " edge_diameter_bound=" << result.edge_diameter_bound << '\n';
        {
            const Moments pilot = phase(graph, model, options, threads, tau, true, options.alpha);
            result.pilot_bags = pilot.bags;
            result.pilot_paths = pilot.paths;
            const long double m = static_cast<long double>(pilot.bags);
            long double max_mean = 0;
            long double max_second = 0;
            for (std::size_t v = 0; v < n; ++v) {
                max_mean = std::max(max_mean, pilot.sums[v].get() / m);
                max_second = std::max(max_second, pilot.squares[v].get() / m);
            }
            result.max_mean_upper = percis_mean_upper(std::min(1.0L, max_mean), tau,
                                                      options.delta / 8);
            result.second_moment_upper = percis_mean_upper(std::min(1.0L, max_second), tau,
                                                           options.delta / 8);
            result.rho_upper = percis_rho_upper(pilot.lengths.get() / m,
                sample_variance(pilot.lengths, pilot.squared_lengths, tau),
                tau, internal_bound, options.delta / 8);
            // Preserve PercIS's optional post-pilot single-path heuristic.
            // If it changes the bag rule, use E[Z_path^2] <= E[Z_path]
            // rather than a second moment estimated for larger pilot bags.
            if (options.bag_size == 0 && !options.alpha_given &&
                result.second_moment_upper > 0.9L * result.max_mean_upper) {
                result.production_alpha = -1;
                result.second_moment_upper = result.max_mean_upper;
            }
            result.variance_upper = std::min(0.25L, result.second_moment_upper);
            result.planned_bags = percis_sample_budget(result.max_mean_upper,
                result.variance_upper, result.rho_upper, normalized_epsilon, options.delta / 2);
        } // The pilot is discarded before collecting the output sample.
        if (!options.quiet)
            std::cerr << "percis: rho_upper=" << result.rho_upper
                      << " max_mean_upper_normalized=" << result.max_mean_upper
                      << " variance_upper_normalized=" << result.variance_upper
                      << " normalized_epsilon=" << normalized_epsilon
                      << " planned_bags=" << result.planned_bags << '\n';
    }
    if (!options.quiet)
        std::cerr << "percis: phase=production bags=" << result.planned_bags << '\n';
    const Moments production = phase(graph, model, options, threads,
                                     result.planned_bags, false, result.production_alpha);
    result.bags = production.bags;
    result.sampled_paths = production.paths;
    result.unreachable_bags = production.unreachable;
    const long double m = static_cast<long double>(result.bags);
    for (std::size_t v = 0; v < n; ++v) {
        result.estimates[v] = model.max_correction * (production.sums[v].get() / m);
        const long double variance = sample_variance(production.sums[v], production.squares[v],
                                                     production.bags);
        result.sample_variances[v] = model.max_correction * (model.max_correction * variance);
        if (!std::isfinite(result.estimates[v]) ||
            (result.bags > 1 && !std::isfinite(result.sample_variances[v])))
            throw std::overflow_error("Reported source-only moments exceed the numeric range");
    }
    return result;
}

} // namespace source_pc
