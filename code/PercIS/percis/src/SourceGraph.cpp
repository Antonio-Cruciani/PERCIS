#include "SourceOnly.h"

#include <algorithm>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace source_pc {
namespace {

std::string location(const std::string& file, std::size_t line) {
    return file + ":" + std::to_string(line) + ": ";
}

std::string uncomment(std::string line) {
    const auto comment = line.find('#');
    if (comment != std::string::npos) line.erase(comment);
    return line;
}

std::uint64_t parse_id(const std::string& token, const std::string& where) {
    std::uint64_t value = 0;
    const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
    if (result.ec != std::errc() || result.ptr != token.data() + token.size()) {
        throw std::runtime_error(where + "vertex IDs must be unsigned decimal integers");
    }
    return value;
}

std::vector<long double> read_states(const std::string& file) {
    std::ifstream input(file);
    if (!input) throw std::runtime_error("cannot open percolation states file: " + file);
    std::vector<long double> states;
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream fields(uncomment(line));
        std::string token, extra;
        if (!(fields >> token)) continue;
        if (fields >> extra) {
            throw std::runtime_error(location(file, line_number) +
                                     "expected exactly one percolation state");
        }
        errno = 0;
        char* end = nullptr;
        const long double value = std::strtold(token.c_str(), &end);
        if (end != token.c_str() + token.size() || end == token.c_str() ||
            errno == ERANGE || !std::isfinite(value) || value < 0 || value > 1) {
            throw std::runtime_error(location(file, line_number) +
                                     "percolation state must be a finite number in [0,1]");
        }
        if (states.size() > std::numeric_limits<Vertex>::max()) {
            throw std::runtime_error("too many vertices in percolation states file: " + file);
        }
        states.push_back(value);
    }
    if (input.bad()) throw std::runtime_error("error reading percolation states file: " + file);
    if (states.empty()) throw std::runtime_error("percolation states file contains no states: " + file);
    return states;
}

} // namespace

GraphData GraphData::load(const std::string& graph_file,
                          const std::string& states_file,
                          bool directed, bool indexed_ids) {
    GraphData graph;
    graph.directed = directed;
    graph.states = read_states(states_file);
    std::ifstream input(graph_file);
    if (!input) throw std::runtime_error("cannot open graph edge-list file: " + graph_file);

    std::unordered_map<std::uint64_t, Vertex> internal_ids;
    if (indexed_ids) {
        graph.adjacency.resize(graph.states.size());
        graph.original_ids.resize(graph.states.size());
        for (std::size_t i = 0; i < graph.states.size(); ++i) graph.original_ids[i] = i;
    }
    const auto map_id = [&](std::uint64_t original) -> Vertex {
        const auto found = internal_ids.find(original);
        if (found != internal_ids.end()) return found->second;
        if (graph.original_ids.size() > std::numeric_limits<Vertex>::max()) {
            throw std::runtime_error("too many distinct vertices in graph: " + graph_file);
        }
        const Vertex id = static_cast<Vertex>(graph.original_ids.size());
        internal_ids.emplace(original, id);
        graph.original_ids.push_back(original);
        graph.adjacency.emplace_back();
        return id;
    };

    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream fields(uncomment(line));
        std::string first, second, extra;
        if (!(fields >> first)) continue;
        if (!(fields >> second) || (fields >> extra)) {
            throw std::runtime_error(location(graph_file, line_number) +
                                     "expected exactly two vertex IDs");
        }
        const auto where = location(graph_file, line_number);
        const auto source = parse_id(first, where);
        const auto target = parse_id(second, where);
        if (indexed_ids && (source >= graph.size() || target >= graph.size())) {
            throw std::runtime_error(where + "vertex ID exceeds the state-file index range");
        }
        // Upstream ignores a self-loop before assigning either endpoint an ID.
        if (source == target) continue;
        const Vertex u = indexed_ids ? static_cast<Vertex>(source) : map_id(source);
        const Vertex v = indexed_ids ? static_cast<Vertex>(target) : map_id(target);
        graph.adjacency[u].push_back(v);
        if (!directed) graph.adjacency[v].push_back(u);
    }
    if (input.bad()) throw std::runtime_error("error reading graph edge-list file: " + graph_file);
    if (graph.size() != graph.states.size()) {
        throw std::runtime_error("state count (" + std::to_string(graph.states.size()) +
                                 ") does not match graph vertex count (" +
                                 std::to_string(graph.size()) +
                                 "); default IDs follow first appearance in the edge list; "
                                 "use --indexed-ids to include isolated vertices");
    }
    for (auto& neighbors : graph.adjacency) {
        std::sort(neighbors.begin(), neighbors.end());
        neighbors.erase(std::unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
    return graph;
}

} // namespace source_pc
