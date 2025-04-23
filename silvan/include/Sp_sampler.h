#ifndef SP_SAMPLER_H
#define SP_SAMPLER_H

#include <vector>
#include <map>
#include <ctype.h>
#include <random>

#include "Rand_gen.h"
#include "Graph.h"
#include "SamplingPreprocessing.h"

class Sp_sampler
{
    public:
        Sp_sampler( const Graph *g, const uint32_t seed ,const double sum_perc, const bool uniform,const SamplingPreprocessing& sampling_kernel,std::vector<double>& sorted_X);
        virtual ~Sp_sampler();
        std::map<uint32_t, double>/*vector<uint32_t>*/ random_path(int &path_length , int &num_paths , double alpha_sp_sampling);
        uint64_t vis_edges;
        std::vector<double> percolation_states;
        std::vector<double> sorted_X;        //Sorted percolation states
        SamplingPreprocessing sampling_kernel;
    protected:
    private:
        void backtrack_path( const uint32_t u, const uint32_t v, const uint32_t start, std::vector<uint32_t> &path );
        void backtrack_all_paths( const uint32_t u, const uint32_t v, const uint32_t start, std::vector<uint32_t> &path );
        inline uint32_t random_node() const;
        //std::vector<double> build_outgoing_weights(const std::vector<double>& X);
        std::pair<int, int> non_uniform_sampling(const std::vector<double>& X, const SamplingPreprocessing& preproc, std::mt19937& rng);
        double denominator_kernel;
        Graph *pred;
        uint32_t *ball_indicator;
        uint32_t *dist;
        uint32_t *q;
        uint64_t *n_paths;
        // Kernel used for the non uniform sampling
        //std::vector<double> sampling_kernel;
        //SamplingPreprocessing sampling_kernel;
        bool uniform;
        Rand_gen *randgen;
        const Graph *g;
};

#endif // SP_SAMPLER_H
