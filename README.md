# Estimating Percolation Centrality with Progressive Sampling and Non-uniform Sampling #

This repository contains the implementation of our approximation algorithm PercIS, supporting our paper ["Fast Percolation Centrality Approximation
with Importance Sampling"]().

PercIS is an algorithm to compute approximations of the Percolation Centralities from graphs with non-uniform sampling. More details can be found in the paper, 


## Installation

The software requires the [OpenMP API](http://openmp.org/wp/). After cloning this repository,
build the software by issuing the `make` command inside the silvan folder.

### Running PercIS ###

The network has to be provided as a file containing two space-separated
integers `u v` per line, for each edge `(u,v)` of the network. The labels of
the vertices are assumed to be consecutive. The percolation states has to be provided as a single column file were line `i` represents the percolation state of node `i` in the graph.

To run PercIS, you can use the julia scripts found in the `julia_scripts` folder.



