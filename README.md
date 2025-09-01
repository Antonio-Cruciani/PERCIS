# Estimating Percolation Centrality with Progressive Sampling and Non-uniform Sampling #

This repository contains the implementation of our approximation algorithm PercIS, supporting our paper ["Fast Percolation Centrality Approximation
with Importance Sampling"]().

PercIS is an algorithm to compute approximations of the Percolation Centralities from graphs with non-uniform sampling. More details can be found in the paper, 


## Installation

The software requires the [OpenMP API](http://openmp.org/wp/). After cloning this repository,
build the software by issuing the `make` command inside the silvan folder.

To run the scripts you need to install ``Julia``.

### Installing Julia
Download the ``julia`` lastest version at: https://julialang.org/



### Running PercIS and reproduce the experiments ###

The network has to be provided as a file containing two space-separated
integers `u v` per line, for each edge `(u,v)` of the network. The labels of
the vertices are assumed to be consecutive. The percolation states has to be provided as a single column file were line `i` represents the percolation state of node `i` in the graph.

To reproduce the experiments in the paper:

	1) Unzip the ``datasets.zip`` folder (keep it in the same level as this README file)
	2) Compile PerIS and Silvan:
		a) ``cd code/PercIS/percis/`` and run ``g++ -fopenmp -std=c++17 -Ofast -Wall -Iinclude main.cpp src/* -lm -o percis``
		b) ``cd code/PercIS/silvan/`` and run ``g++ -fopenmp -std=c++11 -Ofast -Wall -g -Iinclude main.cpp src/* -lm -o silvan``
	2) Move to ``code/PercIS/julia_scripts/`` folder (``cd code/PercIS/julia_scripts/``).
	3) Run the scripts with the command ``julia <script_name>.jl`` (for example ``julia run_percis.jl``).

### Running the Exact Algorithm ###

You will need to install some dependencies to run the exact algorithm. To do so, type ``julia`` in your terminal and execute the following commands:

	1) using Pkg
	2) Pkg.add(["Graphs","StatsBase","DataStructures","Logging","Random"])
	3) exit(1)

To run the exact algorithm, move to the folder: `exact_algorithm/`, modify the ``julia`` file called ``run_algorithm.jl`` and set the dataset for which you want to compute the exact percolation centrality. The dataset name and the percolation states files must have the same name. 

Once you modified the file, run ``julia --threads <number_of_threads_to_use> run_algorithm.jl``. 

## Running PercIS
To run PercIS on other instances, move to PercIS's folder and run 
``./percis -h``. This command will dislpay the helper function's output.

``
compute percolation centrality approximations for all nodes
USAGE: 

	./percis [-fudhm] [-f vc_dimension ] [-i data independent] [-u uniform sampling] [-v verbosity] [-k k_value] [-e exact] [-w window] [-o output] [-a a_emp_peeling] [-s alpha]  [-g sample_size] [-b linear_sampler] [-t thread_number] epsilon delta percolation_states graph
		-f: use the pseudo dimension upper bound on the sample size
		-i: use the data independent upper bound on the sample size
		-u: use uniform sampling for the approximation
		-d: consider the graph as directed
		-k: compute the top-k betweenness centralities (if 0, compute all of them with absolute error) 
		-h: print this help message
		-v: print additional messages (verbosity is the time in second between subsequent outputs)
		-o: path for the output file (if empty, do not write the output)
		-e: path for the exact percolation file (can be empty)
		-w: Sampling window (use only if -e)
		-a: parameter a for empirical peeling (def. = 2)
		-s: parameter alpha for sampling shortest paths (def. = 2.3)
		-m: disable the computation of m_hat
		-g: sample size in case of fixed sample size (or upper bound on samples for -e)
		-b: use the linear time non uniform sampler
		-t: set the thread number to use (max threads as default)
		err: accuracy (0 < epsilon < 1), relative accuracy if k > 0
		delta: confidence (0 < delta < 1)
		percolation states file
		graph: graph edge list file

``

As an example consider running PercIS on the ``01_musae_facebook_edges_e_log`` (undirected) graph with epsilon 0.001 and delta 0.05 and to print the output in a folder called ``path_to_folder/tmp/``:


``./percis -o "path_to_folder/tmp/" 0.001 0.05 "../../datasets/percolation_states/`01_musae_facebook_edges_e_log.txt ../../datasets/graphs/`01_musae_facebook_edges_e_log.txt" ``
