#include "SourceOnly.h"

#include <boost/multiprecision/cpp_int.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

namespace source_pc {
namespace {

using Count = boost::multiprecision::cpp_int;

// Uniform on [0, bound), including bounds beyond the machine word size.
// Rejection preserves the distribution; neither counts nor ranks are rounded.
Count random_below(const Count& bound, std::mt19937_64& rng) {
    if (bound <= 0) {
        throw std::logic_error("A random integer bound must be positive");
    }
    if (bound == 1) return 0;
    if (bound <= std::numeric_limits<std::uint64_t>::max()) {
        const auto upper = bound.convert_to<std::uint64_t>();
        return std::uniform_int_distribution<std::uint64_t>(0, upper - 1)(rng);
    }

    const Count largest = bound - 1;
    const std::size_t bits =
        static_cast<std::size_t>(boost::multiprecision::msb(largest)) + 1;
    for (;;) {
        Count rank = 0;
        std::size_t remaining = bits;
        while (remaining >= 64) {
            rank <<= 64;
            rank += rng();
            remaining -= 64;
        }
        if (remaining != 0) {
            rank <<= remaining;
            rank += rng() >> (64 - remaining);
        }
        if (rank < bound) return rank;
    }
}

// x = mantissa * 2^exponent, exactly for the represented binary state.
// Removing powers of two keeps the common-integer CDF small for simple states.
struct BinaryWeight {
    Count mantissa = 0;
    int exponent = 0;
};

BinaryWeight binary_weight(long double x) {
    if (x == 0) return {};
    int exponent = 0;
    const long double fraction = std::frexp(x, &exponent);
    constexpr int digits = std::numeric_limits<long double>::digits;
    Count mantissa(std::ldexp(fraction, digits));
    const unsigned zeros = boost::multiprecision::lsb(mantissa);
    mantissa >>= zeros;
    return {std::move(mantissa), exponent - digits + static_cast<int>(zeros)};
}

} // namespace

struct Sampler::Impl {
    const GraphData& graph;
    const Model& model;
    const std::uint64_t fixed_bag_size;
    const long double alpha;
    BinaryWeight alpha_weight;
    std::mt19937_64 rng;
    std::vector<Count> source_cdf;
    std::vector<std::vector<Vertex>> incoming;

    static constexpr std::size_t unseen = std::numeric_limits<std::size_t>::max();
    std::vector<unsigned char> owner;
    std::vector<std::size_t> distance;
    std::vector<Count> paths;
    std::vector<std::vector<Vertex>> predecessors;
    // Each vertex belongs to one search ball and enters this reset list once.
    std::vector<Vertex> queue;
    std::vector<Vertex> source_frontier;
    std::vector<Vertex> target_frontier;
    std::vector<Vertex> next_frontier;
    std::vector<std::uint64_t> hits;
    std::vector<Vertex> hit_vertices;

    Impl(const GraphData& graph_, const Model& model_, std::uint64_t seed,
         std::uint64_t bag_size_, long double alpha_)
        : graph(graph_), model(model_), fixed_bag_size(bag_size_), alpha(alpha_), rng(seed) {
        const auto n = graph.size();
        if (!std::isfinite(alpha)) throw std::invalid_argument("Path-sampling alpha must be finite");
        if (alpha > 0) alpha_weight = binary_weight(alpha);
        if (n < 3) throw std::invalid_argument("Source-only centrality needs at least three vertices");
        if (n - 1 > std::numeric_limits<Vertex>::max()) {
            throw std::invalid_argument("Graph vertex indices exceed the supported range");
        }
        if (graph.states.size() != n || model.correction.size() != n) {
            throw std::invalid_argument("Graph states and model must match the graph size");
        }
        for (std::size_t v = 0; v < n; ++v) {
            if (!std::isfinite(graph.states[v]) || graph.states[v] < 0 || graph.states[v] > 1) {
                throw std::invalid_argument("Percolation states must be finite and in [0,1]");
            }
            if (!std::isfinite(model.correction[v]) || model.correction[v] <= 0) {
                throw std::invalid_argument("Importance corrections must be finite and positive");
            }
            for (Vertex w : graph.adjacency[v]) {
                if (static_cast<std::size_t>(w) >= n) {
                    throw std::invalid_argument("Graph contains an out-of-range neighbor");
                }
            }
        }
        if (graph.directed) {
            incoming.resize(n);
            for (std::size_t v = 0; v < n; ++v) {
                for (Vertex w : graph.adjacency[v]) {
                    incoming[w].push_back(static_cast<Vertex>(v));
                }
            }
        }

        // A floating prefix sum can silently erase a small positive source
        // following a large state. An exact integer CDF retains its support.
        std::vector<BinaryWeight> weights;
        weights.reserve(n);
        int minimum_exponent = std::numeric_limits<int>::max();
        for (long double state : graph.states) {
            weights.push_back(binary_weight(state));
            if (weights.back().mantissa != 0) {
                minimum_exponent = std::min(minimum_exponent, weights.back().exponent);
            }
        }
        if (minimum_exponent == std::numeric_limits<int>::max()) {
            throw std::invalid_argument("At least one source must have positive state");
        }
        source_cdf.reserve(n);
        Count cumulative = 0;
        for (const auto& weight : weights) {
            if (weight.mantissa != 0) {
                const auto shift = static_cast<unsigned>(weight.exponent - minimum_exponent);
                cumulative += weight.mantissa << shift;
            }
            source_cdf.push_back(cumulative);
        }

        owner.assign(n, 0);
        distance.assign(n, unseen);
        paths.resize(n);
        predecessors.resize(n);
        hits.assign(n, 0);
    }

    void reset_workspace() {
        for (Vertex v : queue) {
            owner[v] = 0;
            distance[v] = unseen;
            paths[v] = 0;
            predecessors[v].clear();
        }
        queue.clear();
        source_frontier.clear();
        target_frontier.clear();
        next_frontier.clear();
        for (Vertex v : hit_vertices) hits[v] = 0;
        hit_vertices.clear();
    }

    Bag sample() {
        const Count rank = random_below(source_cdf.back(), rng);
        const auto it = std::upper_bound(source_cdf.begin(), source_cdf.end(), rank);
        if (it == source_cdf.end()) throw std::logic_error("Invalid source sampling rank");
        const Vertex source = static_cast<Vertex>(it - source_cdf.begin());
        const auto reduced_target = std::uniform_int_distribution<std::uint64_t>(
            0, static_cast<std::uint64_t>(graph.size() - 2))(rng);
        const Vertex target = static_cast<Vertex>(
            reduced_target + (reduced_target >= source ? 1 : 0));
        return sample_pair(source, target);
    }

    void add_hit(Vertex v, std::uint64_t increment) {
        if (hits[v] == 0) hit_vertices.push_back(v);
        if (hits[v] > std::numeric_limits<std::uint64_t>::max() - increment) {
            throw std::overflow_error("Bag hit count overflow");
        }
        hits[v] += increment;
    }

    std::uint64_t bag_size_for(const Count& count) const {
        if (fixed_bag_size != 0) return fixed_bag_size;
        if (alpha <= 0 || count == 1) return 1;
        // Preserve upstream's floating multiply followed by truncation on its
        // ordinary uint64_t count domain (e.g. double(2.3)*10 rounds to 23).
        // Check before conversion; the legacy int cast could overflow or give
        // an empty bag. Wider counts use the exact binary product below.
        if (count <= std::numeric_limits<std::uint64_t>::max()) {
            const double product = static_cast<double>(alpha) *
                static_cast<double>(count.convert_to<std::uint64_t>());
            if (std::isfinite(product) && product < std::ldexp(1.0, 64)) {
                return product < 1 ? 1 : static_cast<std::uint64_t>(product);
            }
        }
        Count scaled = count * alpha_weight.mantissa;
        if (alpha_weight.exponent >= 0) scaled <<= alpha_weight.exponent;
        else scaled >>= -alpha_weight.exponent;
        if (scaled < 1) return 1;
        if (scaled > std::numeric_limits<std::uint64_t>::max()) {
            throw std::overflow_error(
                "PercIS alpha times the shortest-path count exceeds uint64_t; use --bag-size to set a fixed bag size");
        }
        return scaled.convert_to<std::uint64_t>();
    }

    const std::vector<Vertex>& backward_neighbors(Vertex v) const {
        return graph.directed ? incoming[v] : graph.adjacency[v];
    }

    static std::uint64_t add_degree(std::uint64_t total, std::size_t degree) {
        // Degree sums only choose which frontier to expand. Saturation cannot
        // change the sampled law; frontier emptiness is checked separately.
        const auto maximum = std::numeric_limits<std::uint64_t>::max();
        return degree > maximum - total ? maximum : total + degree;
    }

    // Each rank identifies a path from one root to a bridge endpoint. At each predecessor,
    // subtract the sizes of the preceding prefix-path blocks and keep the
    // residual rank. This is exact unranking, with no recursive stack.
    void accumulate_side(Vertex root, Vertex start, Count rank,
                         std::uint64_t increment) {
        Vertex current = start;
        if (current != root) add_hit(current, increment);
        while (current != root) {
            bool selected = false;
            for (Vertex previous : predecessors[current]) {
                if (rank < paths[previous]) {
                    current = previous;
                    selected = true;
                    break;
                }
                rank -= paths[previous];
            }
            if (!selected) throw std::logic_error("Shortest-path rank exceeds predecessor counts");
            if (current != root) add_hit(current, increment);
        }
    }

    Bag sample_pair(Vertex source, Vertex target) {
        if (static_cast<std::size_t>(source) >= graph.size() ||
            static_cast<std::size_t>(target) >= graph.size()) {
            throw std::out_of_range("Sample endpoints must be graph vertices");
        }
        if (source == target) throw std::invalid_argument("Sample endpoints must be distinct");
        reset_workspace();
        Bag bag;
        bag.source = source;
        bag.target = target;

        owner[source] = 1;
        owner[target] = 2;
        distance[source] = distance[target] = 0;
        paths[source] = paths[target] = 1;
        queue.push_back(source);
        queue.push_back(target);
        source_frontier.push_back(source);
        target_frontier.push_back(target);
        std::uint64_t source_degree = graph.adjacency[source].size();
        std::uint64_t target_degree = backward_neighbors(target).size();

        struct Bridge {
            Vertex left;
            Vertex right;
            Count count;
        };
        std::vector<Bridge> bridges;
        std::size_t shortest_distance = unseen;

        // As in PercIS, expand a whole frontier from the ball with the smaller
        // degree sum. Reverse edges are used for the target ball. Finish the
        // meeting level so every shortest path is represented by one bridge.
        while (!source_frontier.empty() && !target_frontier.empty()) {
            const bool forward = source_degree <= target_degree;
            auto& frontier = forward ? source_frontier : target_frontier;
            next_frontier.clear();
            std::uint64_t next_degree = 0;
            for (Vertex v : frontier) {
                const auto& neighbors = forward ? graph.adjacency[v] : backward_neighbors(v);
                for (Vertex w : neighbors) {
                    if (owner[w] == 0) {
                        owner[w] = owner[v];
                        distance[w] = distance[v] + 1;
                        paths[w] = paths[v];
                        predecessors[w].push_back(v);
                        queue.push_back(w);
                        next_frontier.push_back(w);
                        next_degree = add_degree(next_degree, forward
                            ? graph.adjacency[w].size() : backward_neighbors(w).size());
                    } else if (owner[w] == owner[v]) {
                        if (distance[w] == distance[v] + 1) {
                            paths[w] += paths[v];
                            predecessors[w].push_back(v);
                        }
                    } else {
                        const Vertex left = forward ? v : w;
                        const Vertex right = forward ? w : v;
                        const std::size_t candidate = distance[left] + distance[right] + 1;
                        if (candidate < shortest_distance) {
                            shortest_distance = candidate;
                            bridges.clear();
                        }
                        if (candidate == shortest_distance) {
                            bridges.push_back({left, right, 0});
                        }
                    }
                }
            }
            if (!bridges.empty()) break;
            frontier.swap(next_frontier);
            (forward ? source_degree : target_degree) = next_degree;
        }
        if (bridges.empty()) return bag;

        Count total_paths = 0;
        for (auto& bridge : bridges) {
            bridge.count = paths[bridge.left] * paths[bridge.right];
            total_paths += bridge.count;
        }
        const std::uint64_t size = bag_size_for(total_paths);

        bag.reachable = true;
        bag.sampled_paths = size;
        if (shortest_distance - 1 > std::numeric_limits<std::uint32_t>::max()) {
            throw std::overflow_error("Shortest-path length exceeds the supported range");
        }
        bag.internal_vertices = static_cast<std::uint32_t>(shortest_distance - 1);
        // Direct paths have no internal vertices. A unique path has a
        // deterministic bag mean, irrespective of the positive bag size.
        if (bag.internal_vertices == 0) return bag;
        const auto accumulate = [&](Count rank, std::uint64_t increment) {
            for (const auto& bridge : bridges) {
                if (rank < bridge.count) {
                    const Count prefix_rank = rank / paths[bridge.right];
                    const Count suffix_rank = rank % paths[bridge.right];
                    accumulate_side(source, bridge.left, prefix_rank, increment);
                    accumulate_side(target, bridge.right, suffix_rank, increment);
                    return;
                }
                rank -= bridge.count;
            }
            throw std::logic_error("Shortest-path rank exceeds meeting-edge counts");
        };
        if (total_paths == 1) {
            accumulate(0, size);
        } else {
            for (std::uint64_t draw = 0; draw < size; ++draw) {
                accumulate(random_below(total_paths, rng), 1);
            }
        }
        bag.scores.reserve(hit_vertices.size());
        for (Vertex v : hit_vertices) {
            const long double fraction = static_cast<long double>(hits[v]) /
                                         static_cast<long double>(size);
            bag.scores.emplace_back(v, model.correction[v] * fraction);
        }
        return bag;
    }
};

Sampler::Sampler(const GraphData& graph, const Model& model,
                 std::uint64_t seed, std::uint64_t bag_size, long double alpha)
    : impl_(std::make_unique<Impl>(graph, model, seed, bag_size, alpha)) {}

Sampler::~Sampler() = default;
Sampler::Sampler(Sampler&&) noexcept = default;
Sampler& Sampler::operator=(Sampler&&) noexcept = default;

Bag Sampler::sample() { return impl_->sample(); }
Bag Sampler::sample_pair(Vertex source, Vertex target) {
    return impl_->sample_pair(source, target);
}

} // namespace source_pc
