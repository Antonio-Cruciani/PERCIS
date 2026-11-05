// Source-only percolation centrality with independent shortest-path bags.
#include "SourceOnly.h"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <getopt.h>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <system_error>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace {

void usage(std::ostream& out, const char* program) {
    out << "Usage: " << program
        << " [options] epsilon delta states_file graph_file\n\n"
        << "Estimate SOURCE-ONLY percolation centrality. Each independent sample\n"
        << "draws s proportional to x_s, t uniformly excluding s, and a bag of\n"
        << "uniform shortest paths with replacement. Unreachable pairs count as zero.\n\n"
        << "  -d, --directed           Read directed edges (default: undirected)\n"
        << "  -g, --samples N          Use exactly N independent bags\n"
        << "  -s, --alpha A            PercIS path-sampling factor (default: 2.3)\n"
        << "                          A<=0 selects one shortest path per pair\n"
        << "      --bag-size B         Override with B paths per reachable bag\n"
        << "      --pilot-samples N    Override PercIS's pilot size (N>=2)\n"
        << "      --seed N             Random seed (default: 42)\n"
        << "  -t, --threads N          Worker threads (0/default: all available)\n"
        << "  -o, --output FILE        Write internal_id,estimate CSV (default: stdout)\n"
        << "      --diagnostics FILE   Write estimates and production bag variances\n"
        << "      --id-map FILE        Write internal_id,original_id CSV\n"
        << "      --indexed-ids        Vertex IDs are state-row indices 0,...,n-1;\n"
        << "                          includes isolated vertices\n"
        << "  -q, --quiet              Suppress phase progress\n"
        << "  -h, --help               Show this help\n\n"
        << "Without -g, use PercIS: pilot -> one computed production budget ->\n"
        << "fresh production samples. No adaptive stopping or confidence checkpoints.\n"
        << "Both epsilon and delta must be in (0,1). With -g they do not select the\n"
        << "budget. The pilot is discarded; output uses production bags only.\n"
        << "Default bags contain max(1,floor(alpha*sigma_st)) paths when sigma_st>1,\n"
        << "and one path otherwise. --bag-size and --alpha are mutually exclusive.\n"
        << "States: one number in [0,1] per vertex, at least two positive states.\n"
        << "Default vertex order: first appearance in the edge list, as in upstream.\n";
}

std::uint64_t integer(const char* text, const std::string& name, bool allow_zero = false) {
    const std::string value(text);
    std::uint64_t parsed = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), parsed);
    if (result.ec != std::errc() || result.ptr != value.data() + value.size() ||
        value.empty() || (!allow_zero && parsed == 0)) {
        throw std::invalid_argument(name + " must be " +
            (allow_zero ? "a nonnegative" : "a positive") + " 64-bit integer");
    }
    return parsed;
}

long double probability(const char* text, const std::string& name) {
    const std::string value(text);
    std::size_t consumed = 0;
    long double parsed = 0;
    try {
        parsed = std::stold(value, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument(name + " must be a finite number in (0,1)");
    }
    if (consumed != value.size() || !std::isfinite(parsed) || parsed <= 0 || parsed >= 1) {
        throw std::invalid_argument(name + " must be a finite number in (0,1)");
    }
    return parsed;
}

long double alpha_value(const char* text) {
    std::size_t consumed = 0;
    const std::string value(text);
    long double parsed;
    try {
        parsed = std::stold(value, &consumed);
    } catch (const std::exception&) {
        throw std::invalid_argument("alpha must be a finite number");
    }
    if (consumed != value.size() || !std::isfinite(parsed))
        throw std::invalid_argument("alpha must be a finite number");
    return parsed;
}

std::ofstream output_file(const std::string& path) {
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("Cannot open output file: " + path);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << std::setprecision(std::numeric_limits<long double>::max_digits10);
    return stream;
}

bool same_file(const std::string& left, const std::string& right) {
    namespace fs = std::filesystem;
    if (fs::weakly_canonical(left) == fs::weakly_canonical(right)) return true;
    return fs::exists(left) && fs::exists(right) && fs::equivalent(left, right);
}

} // namespace

int main(int argc, char** argv) {
    try {
        source_pc::RunOptions options;
        bool directed = false;
        bool indexed_ids = false;
        std::string output, diagnostics, id_map;
        enum { BAG_SIZE = 1000, SEED, INDEXED_IDS, PILOT_SAMPLES, DIAGNOSTICS, ID_MAP };
        static const option long_options[] = {
            {"directed", no_argument, nullptr, 'd'},
            {"samples", required_argument, nullptr, 'g'},
            {"bag-size", required_argument, nullptr, BAG_SIZE},
            {"alpha", required_argument, nullptr, 's'},
            {"pilot-samples", required_argument, nullptr, PILOT_SAMPLES},
            {"seed", required_argument, nullptr, SEED},
            {"threads", required_argument, nullptr, 't'},
            {"output", required_argument, nullptr, 'o'},
            {"diagnostics", required_argument, nullptr, DIAGNOSTICS},
            {"id-map", required_argument, nullptr, ID_MAP},
            {"indexed-ids", no_argument, nullptr, INDEXED_IDS},
            {"quiet", no_argument, nullptr, 'q'},
            {"help", no_argument, nullptr, 'h'},
            {nullptr, 0, nullptr, 0}
        };
        opterr = 0;
        int flag;
        while ((flag = getopt_long(argc, argv, "dg:t:o:s:qh", long_options, nullptr)) != -1) {
            switch (flag) {
                case 'd': directed = true; break;
                case 'g': options.samples = integer(optarg, "samples"); break;
                case 's': options.alpha = alpha_value(optarg); options.alpha_given = true; break;
                case 't': {
                    const auto threads = integer(optarg, "threads", true);
                    if (threads > static_cast<std::uint64_t>(std::numeric_limits<int>::max()))
                        throw std::invalid_argument("threads is too large");
                    options.threads = static_cast<int>(threads);
                    break;
                }
                case 'o': output = optarg; break;
                case 'q': options.quiet = true; break;
                case 'h': usage(std::cout, argv[0]); return 0;
                case BAG_SIZE: options.bag_size = integer(optarg, "bag-size"); break;
                case SEED: options.seed = integer(optarg, "seed", true); break;
                case INDEXED_IDS: indexed_ids = true; break;
                case PILOT_SAMPLES: options.pilot_samples = integer(optarg, "pilot-samples"); break;
                case DIAGNOSTICS: diagnostics = optarg; break;
                case ID_MAP: id_map = optarg; break;
                default: throw std::invalid_argument("Unknown option or missing argument; use --help");
            }
        }
        if (argc - optind != 4) {
            usage(std::cerr, argv[0]);
            return 1;
        }
        options.epsilon = probability(argv[optind], "epsilon");
        options.delta = probability(argv[optind + 1], "delta");
        if (options.samples && options.pilot_samples)
            throw std::invalid_argument("A fixed -g budget skips the pilot; omit --pilot-samples");
        if (options.pilot_samples == 1)
            throw std::invalid_argument("The pilot must contain at least two bags");
        if (options.bag_size && options.alpha_given)
            throw std::invalid_argument("Use either --bag-size or --alpha, not both");
        if (options.threads == 0) {
#ifdef _OPENMP
            options.threads = omp_get_max_threads();
#else
            options.threads = 1;
#endif
        }
        // Separate result files prevent one requested output from overwriting another.
        const std::string states_path = argv[optind + 2];
        const std::string graph_path = argv[optind + 3];
        const std::vector<std::string> outputs = {output, diagnostics, id_map};
        for (std::size_t i = 0; i < outputs.size(); ++i) {
            if (outputs[i].empty()) continue;
            if (same_file(outputs[i], states_path) || same_file(outputs[i], graph_path))
                throw std::invalid_argument("Output files must differ from the input files");
            for (std::size_t j = i + 1; j < outputs.size(); ++j) {
                if (!outputs[j].empty() && same_file(outputs[i], outputs[j]))
                    throw std::invalid_argument("Output, diagnostics, and ID map must use distinct files");
            }
        }

        const auto start = std::chrono::steady_clock::now();
        const auto graph = source_pc::GraphData::load(graph_path, states_path, directed, indexed_ids);
        const source_pc::Model model(graph);
        std::cerr << std::setprecision(10)
                  << "measure=source-only vertices=" << graph.size()
                  << " directed=" << (directed ? "yes" : "no")
                  << " state_sum=" << model.total_state
                  << " max_importance_weight=" << model.max_correction
                  << " bag_size=" << (options.bag_size ? std::to_string(options.bag_size) : "percis")
                  << " alpha=" << options.alpha << " seed=" << options.seed
                  << " threads=" << options.threads << '\n';
        const auto result = source_pc::estimate(graph, model, options);

        std::ofstream file;
        std::ostream* estimates = &std::cout;
        if (!output.empty()) {
            file = output_file(output);
            estimates = &file;
        }
        *estimates << std::setprecision(std::numeric_limits<long double>::max_digits10);
        for (std::size_t v = 0; v < graph.size(); ++v)
            *estimates << v << ',' << result.estimates[v] << '\n';
        estimates->flush();
        if (!*estimates) throw std::runtime_error("Failed to write estimates");
        if (!diagnostics.empty()) {
            auto detail = output_file(diagnostics);
            detail << "internal_id,estimate,bag_sample_variance\n";
            for (std::size_t v = 0; v < graph.size(); ++v)
                detail << v << ',' << result.estimates[v] << ','
                       << result.sample_variances[v] << '\n';
            detail.close();
        }
        if (!id_map.empty()) {
            auto mapping = output_file(id_map);
            for (std::size_t v = 0; v < graph.size(); ++v)
                mapping << v << ',' << graph.original_ids[v] << '\n';
            mapping.close();
        }
        if (file.is_open()) file.close();
        const double elapsed = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - start).count();
        std::cerr << "mode=" << (options.samples ? "fixed" : "percis")
                  << " pilot_bags=" << result.pilot_bags
                  << " pilot_paths=" << result.pilot_paths
                  << " planned_bags=" << result.planned_bags
                  << " bags=" << result.bags << " inner_paths=" << result.sampled_paths
                  << " unreachable_bags=" << result.unreachable_bags
                  << " production_alpha=" << result.production_alpha
                  << " threads=" << result.threads
                  << " elapsed_seconds=" << elapsed << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "percis-source: " << error.what() << '\n';
        return 1;
    }
}
