#include <iostream>
#include <dirent.h>
#include <limits.h>
#include <math.h>
#include <ctime>
#include <cstddef>
#include <unistd.h>
#include <ctype.h>
#include <fstream>

#include "Rand_gen.h"
#include "Graph.h"
#include "utilities.h"
#include "Probabilistic.h"

extern char *optarg;
static const std::string ERROR_HEADER = "ERROR: ";
using namespace std;

bool directed = false;
double verb = 60;
double delta;
double err;
char *graph_file;
char *percolation_file;

std::string output_file;
int64_t k = 0;
double sampling_rate = 2.3;
bool alpha_given = false;
double empirical_peeling_param = 2.0;
bool m_hat_enabled = true;
bool uniform = false;
// mcrade
int num_mc = 10;
bool optimized = true;
bool fixed_ss = false;
uint32_t sample_size;
bool vc_dim = false;
/**
 * Print usage on stderr.
 */
void usage(const char *binary_name) {
    std::cerr << binary_name
        << ": compute percolation centrality approximations for all nodes"
        << std::endl;
    std::cerr << "USAGE: " << binary_name << " [-fudhm] [-f vc_dimension ] [-u uniform sampling] [-v verbosity] [-k k_value] [-o output] [-a a_emp_peeling] [-s alpha]  [-g sample_size] [-b linear_sampler] epsilon delta percolation_states graph"
        << std::endl;
    std::cerr << "\t-f: use the vc dimension upper bound on the sample size" << std::endl;
    std::cerr << "\t-u: use uniform sampling for the approximation" << std::endl;
    std::cerr << "\t-d: consider the graph as directed" << std::endl;
    std::cerr << "\t-k: compute the top-k betweenness centralities (if 0, compute all of them with absolute error) " << std::endl;
    std::cerr << "\t-h: print this help message" << std::endl;
    std::cerr << "\t-v: print additional messages (verbosity is the time in second between subsequent outputs)" << std::endl;
    std::cerr << "\t-o: path for the output file (if empty, do not write the output)" << std::endl;
    std::cerr << "\t-a: parameter a for empirical peeling (def. = 2)" << std::endl;
    std::cerr << "\t-s: parameter alpha for sampling shortest paths (def. = 2.3)" << std::endl;
    std::cerr << "\t-m: disable the computation of m_hat" << std::endl;
    std::cerr << "\t-g: sample size in case of fixed sample size" << std::endl;
    std::cerr << "\t-b: use the linear time non uniform sampler" << std::endl;
    std::cerr << "\terr: accuracy (0 < epsilon < 1), relative accuracy if k > 0" << std::endl;
    std::cerr << "\tdelta: confidence (0 < delta < 1)" << std::endl;
    std::cerr << "\tpercolation states file" << std::endl;
    std::cerr << "\tgraph: graph edge list file" << std::endl;
}

/**
 * Parse command line options.
 * Return 0 if everything went well, 1 if there were errors, 2 if -h was specified.
 */
int parse_command_line(int& argc, char *argv[]) {
    int opt;
    while ((opt = getopt(argc, argv,"bfudhmk:o:s:a:v:g:")) != -1) {
        switch (opt) {
        case 'b':
            optimized = false;
            break;
        case 'f':
            vc_dim = true;
            break;
        case 'u':
            uniform = true;
            cout<<"Selected Uniform Sampling"<<endl;
            break;
        case 'd':
            directed = true;
            break;
        case 'h':
            return 2;
            break;
        case 'm':
            m_hat_enabled = false;
            break;
        case 'o':
            std::cerr << "Writing output to " << optarg << std::endl;
            output_file = optarg;
            break;
        case 's':
            sampling_rate = std::strtod(optarg, NULL);
            alpha_given = true;
            break;
        case 'a':
            empirical_peeling_param = std::strtod(optarg, NULL);
            if (errno == ERANGE || empirical_peeling_param <= 1) {
                std::cerr << ERROR_HEADER
                    << "The value a should be >= 1. Empirical peeling disabled."
                    << std::endl;
            }
            alpha_given = true;
            break;
        case 'k':
            k = std::strtod(optarg, NULL);
            if (errno == ERANGE || k < 0 || k > UINT_MAX) {
                std::cerr << ERROR_HEADER
                    << "The value k should be between 0 and 2^32-1."
                    << std::endl;
                return 1;
            }
            break;
        case 'v':
            verb = std::strtod(optarg, NULL);
            if (errno == ERANGE || verb < 0) {
                std::cerr << ERROR_HEADER
                    << "The verbosity should be a positive number, or 0 to produce no output."
                    << std::endl;
                return 1;
            }
            break;
        case 'g':
            sample_size = std::strtod(optarg, NULL);
            fixed_ss = true;
            if (errno == ERANGE || sample_size < 0) {
                std::cerr << ERROR_HEADER
                        << "The sample size should be a positive number"
                        << std::endl;
                return 1;
            }
        }
       
    }


    if (optind != argc - 4) {
        std::cerr << ERROR_HEADER << "Wrong number of arguments" << std::endl;
        return 1;
    } else {
        err = std::strtod(argv[argc - 4], NULL);
        if (errno == ERANGE || err >= 1.0 || err <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "The error err should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        delta = std::strtod(argv[argc - 3], NULL);
        if (errno == ERANGE || delta >= 1.0 || delta <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "Delta should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        
        percolation_file = argv[argc - 2];
        graph_file = argv[argc - 1];
    }
    
    /*
    if (optind != argc - 3) {
        std::cerr << ERROR_HEADER << "Wrong number of arguments" << std::endl;
        return 1;
    } else {
        err = std::strtod(argv[argc - 3], NULL);
        if (errno == ERANGE || err >= 1.0 || err <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "The error err should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        delta = std::strtod(argv[argc - 2], NULL);
        if (errno == ERANGE || delta >= 1.0 || delta <= 0.0) {
            std::cerr << ERROR_HEADER <<
                "Delta should be greater than 0 and smaller than 1"
                << std::endl;
            return 1;
        }
        graph_file = argv[argc - 1];
        // test if input graph file exists
        std::ifstream infile(graph_file);
        if(!infile.good()){
          std::cerr << "Problems with input graph file" << std::endl;
          return 1;
        }
        percolation_file = argv[argc - 2];
        
        
    }
    */
    return 0;
}

int main(int argc, char *argv[]){
    int correct_parse = parse_command_line(argc, argv);

    if (correct_parse != 0) {
        usage(argv[0]);
        return correct_parse!=2;
    }
    if (!fixed_ss){
        Probabilistic G( graph_file,percolation_file, uniform,optimized, directed, verb , sampling_rate , alpha_given , empirical_peeling_param , m_hat_enabled , vc_dim,output_file);
        G.run((uint32_t) k, delta, err,uniform,optimized,vc_dim);
    }else{
        Probabilistic G( graph_file,percolation_file,sample_size, uniform,optimized,directed, verb ,sampling_rate , alpha_given, empirical_peeling_param , m_hat_enabled , output_file);
        G.run_fixed_sample_size(k,delta,err,sample_size,uniform,optimized);
    }
    std::cout << "run finished" << std::endl;
    return 0;
}
