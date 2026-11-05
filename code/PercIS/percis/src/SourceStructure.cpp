#include "SourceOnly.h"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace source_pc {

std::uint32_t percis_edge_diameter_bound(const GraphData& graph) {
    const std::size_t n = graph.size();
    if (n == 0) return 0;
    const Vertex unseen = std::numeric_limits<Vertex>::max();
    if (n > static_cast<std::size_t>(unseen))
        throw std::invalid_argument("Too many vertices for the diameter bound");

    std::vector<std::vector<Vertex>> reverse(n);
    for (std::size_t u = 0; u < n; ++u) {
        for (Vertex v : graph.adjacency[u]) {
            if (v >= n) throw std::invalid_argument("Invalid adjacency in diameter bound");
            reverse[v].push_back(static_cast<Vertex>(u));
        }
    }

    // Kosaraju's two passes are iterative, including the finishing-order DFS.
    // A long directed path therefore does not consume the process call stack.
    struct Frame {
        Vertex vertex;
        std::size_t next_neighbor;
    };
    std::vector<unsigned char> visited(n, 0);
    std::vector<Vertex> finish;
    finish.reserve(n);
    std::vector<Frame> frames;
    for (std::size_t start = 0; start < n; ++start) {
        if (visited[start]) continue;
        visited[start] = 1;
        frames.push_back({static_cast<Vertex>(start), 0});
        while (!frames.empty()) {
            Frame& frame = frames.back();
            const auto& neighbors = graph.adjacency[frame.vertex];
            if (frame.next_neighbor == neighbors.size()) {
                finish.push_back(frame.vertex);
                frames.pop_back();
            } else {
                const Vertex next = neighbors[frame.next_neighbor++];
                if (!visited[next]) {
                    visited[next] = 1;
                    frames.push_back({next, 0});
                }
            }
        }
    }

    std::vector<Vertex> component(n, unseen);
    std::vector<Vertex> stack;
    std::size_t components = 0;
    for (auto it = finish.rbegin(); it != finish.rend(); ++it) {
        if (component[*it] != unseen) continue;
        const Vertex id = static_cast<Vertex>(components++);
        component[*it] = id;
        stack.push_back(*it);
        while (!stack.empty()) {
            const Vertex vertex = stack.back();
            stack.pop_back();
            for (Vertex next : reverse[vertex]) {
                if (component[next] == unseen) {
                    component[next] = id;
                    stack.push_back(next);
                }
            }
        }
    }

    // Match upstream: among maximum total-degree vertices, use the first
    // internal vertex ID. Degrees include edges that leave or enter the SCC.
    std::vector<Vertex> pivot(components, unseen);
    for (std::size_t vertex = 0; vertex < n; ++vertex) {
        Vertex& selected = pivot[component[vertex]];
        if (selected == unseen ||
            graph.adjacency[vertex].size() + reverse[vertex].size() >
                graph.adjacency[selected].size() + reverse[selected].size()) {
            selected = static_cast<Vertex>(vertex);
        }
    }

    std::vector<Vertex> queue;
    queue.reserve(n);
    const auto eccentricities = [&](const std::vector<std::vector<Vertex>>& adjacency) {
        std::vector<std::uint32_t> distance(n, unseen);
        std::vector<std::uint32_t> eccentricity(components, 0);
        for (std::size_t id = 0; id < components; ++id) {
            queue.clear();
            const Vertex start = pivot[id];
            distance[start] = 0;
            queue.push_back(start);
            for (std::size_t head = 0; head < queue.size(); ++head) {
                const Vertex vertex = queue[head];
                for (Vertex next : adjacency[vertex]) {
                    if (component[next] == id && distance[next] == unseen) {
                        distance[next] = distance[vertex] + 1;
                        queue.push_back(next);
                    }
                }
            }
            eccentricity[id] = distance[queue.back()];
        }
        return eccentricity;
    };
    const auto forward = eccentricities(graph.adjacency);
    const auto backward = eccentricities(reverse);

    // Parallel edges in the component DAG are harmless: count them both in
    // its indegrees and in the corresponding decrements. This keeps O(n+m)
    // construction time without requiring sets or sorting.
    std::vector<std::vector<Vertex>> dag(components);
    std::vector<std::size_t> indegree(components, 0);
    for (std::size_t u = 0; u < n; ++u) {
        for (Vertex v : graph.adjacency[u]) {
            if (component[u] != component[v]) {
                dag[component[u]].push_back(component[v]);
                ++indegree[component[v]];
            }
        }
    }
    std::vector<Vertex> topological;
    topological.reserve(components);
    for (std::size_t id = 0; id < components; ++id) {
        if (indegree[id] == 0) topological.push_back(static_cast<Vertex>(id));
    }
    for (std::size_t head = 0; head < topological.size(); ++head) {
        for (Vertex next : dag[topological[head]]) {
            if (--indegree[next] == 0) topological.push_back(next);
        }
    }
    if (topological.size() != components)
        throw std::logic_error("The SCC condensation graph contains a cycle");

    // The original AllCCUpperBound recurrence, in reverse topological order.
    // Wider arithmetic avoids overflow while combining component bounds.
    std::vector<std::uint64_t> reach(components, 0);
    std::uint64_t diameter = 0;
    for (auto it = topological.rbegin(); it != topological.rend(); ++it) {
        const Vertex id = *it;
        reach[id] = forward[id];
        for (Vertex next : dag[id]) {
            reach[id] = std::max(reach[id], static_cast<std::uint64_t>(forward[id]) +
                1 + backward[next] + reach[next]);
        }
        diameter = std::max(diameter, reach[id] + backward[id]);
    }
    // Even when the recurrence exceeds uint32_t, every finite shortest path
    // has at most n-1 edges, so saturation remains a valid upper bound.
    return static_cast<std::uint32_t>(std::min<std::uint64_t>(diameter, unseen));
}

} // namespace source_pc
