#ifndef PERCIS_SOURCE_ONLY_H
#define PERCIS_SOURCE_ONLY_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace source_pc {

using Vertex = std::uint32_t;

struct GraphData {
    std::vector<std::vector<Vertex>> adjacency;
    std::vector<std::uint64_t> original_ids;
    std::vector<long double> states;
    bool directed = false;
    std::size_t size() const { return adjacency.size(); }
    // Default ID order matches upstream PercIS: first appearance in the edge list.
    // Indexed mode uses state-line indices as IDs, including isolated vertices.
    static GraphData load(const std::string& graph_file,
                          const std::string& states_file,
                          bool directed, bool indexed_ids = false);
};

struct Model {
    std::vector<long double> correction;
    std::vector<long double> cumulative_states;
    long double total_state = 0;
    long double max_correction = 0;
    explicit Model(const GraphData& graph);
};

struct Bag {
    Vertex source = 0;
    Vertex target = 0;
    bool reachable = false;
    std::uint64_t sampled_paths = 0;
    std::uint32_t internal_vertices = 0;
    // Sparse importance-weighted BAG means, with both endpoints excluded.
    std::vector<std::pair<Vertex, long double>> scores;
};

class Sampler {
public:
    Sampler(const GraphData& graph, const Model& model,
            std::uint64_t seed, std::uint64_t bag_size = 0,
            long double alpha = 2.3L);
    ~Sampler();
    Sampler(Sampler&&) noexcept;
    Sampler& operator=(Sampler&&) noexcept;
    Sampler(const Sampler&) = delete;
    Sampler& operator=(const Sampler&) = delete;
    Bag sample();
    // Also exposed for checking uniform shortest-path sampling independently.
    Bag sample_pair(Vertex source, Vertex target);
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

struct RunOptions {
    std::uint64_t samples = 0; // Zero selects PercIS pilot -> fixed production.
    std::uint64_t pilot_samples = 0; // Zero uses PercIS's pilot-size heuristic.
    std::uint64_t bag_size = 0; // Zero uses PercIS's alpha * shortest-path count.
    long double alpha = 2.3L;
    bool alpha_given = false;
    std::uint64_t seed = 42;
    int threads = 0; // Zero uses all available OpenMP workers (one in serial).
    long double epsilon = 0.01L;
    long double delta = 0.05L;
    bool quiet = false;
};

struct RunResult {
    std::vector<long double> estimates;
    std::vector<long double> sample_variances;
    std::uint64_t bags = 0;
    std::uint64_t sampled_paths = 0;
    std::uint64_t unreachable_bags = 0;
    std::uint64_t pilot_bags = 0;
    std::uint64_t pilot_paths = 0;
    std::uint64_t planned_bags = 0;
    std::uint32_t edge_diameter_bound = 0;
    long double rho_upper = 0;
    long double max_mean_upper = 0;
    long double second_moment_upper = 0;
    long double variance_upper = 0;
    long double production_alpha = 2.3L;
    int threads = 1;
};

// PercIS's SCC/pivot upper bound on shortest-path EDGE distance.
std::uint32_t percis_edge_diameter_bound(const GraphData& graph);

RunResult estimate(const GraphData& graph, const Model& model,
                   const RunOptions& options);

} // namespace source_pc
#endif
