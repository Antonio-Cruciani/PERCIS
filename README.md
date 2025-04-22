# Estimating Percolation Centrality with Progressive Sampling and Non-uniform Sampling #

This repository contains the implementation of our approximation algorithm APERITIF, supporting our paper ["APERITIF: Estimating Percolation Centrality with Progressive Sampling and
Non-uniform Sampling"]().

APERITIF is an algorithm to compute approximations of the Percolation Centralities from graphs with non uniform progressive sampling. More details can be found in the paper, 


Part of the underlying implementation of APERITIF is based on the sampling algorithm [SILVAN](https://github.com/VandinLab/SILVAN) by Leonardo Pellegrina and Fabio Vandin. Therefore, it has the same compilation dependencies (described below) and it is distributed with the same license (Apache License 2.0).

## Installation

The software requires the [OpenMP API](http://openmp.org/wp/). After cloning this repository,
build the software by issuing the `make` command inside the silvan folder.

### Running APERITIF ###

The network has to be provided as a file containing two space-separated
integers `u v` per line, for each edge `(u,v)` of the network. The labels of
the vertices are assumed to be consecutive.

To run APERITIF, you can use the `run_aperitif.py` python script found in the `scripts` folder. It takes the following input parameters:


(TODO)



For example, to approximate the betweenness centrality of all nodes of the undirected graph `graph.txt` with absolute accuracy 0.01 and with probability at least 0.95, you can use

`python run_aperitif.py -db graph.txt -e 0.01 -d 0.05`

or to approximate the top-10 most central nodes of the directed graph `digraph.txt` with relative accuracy 0.1, you can use

`python run_aperitif.py -db digraph.txt -e 0.1 -k 10 -t 1`

### Reproducing experiments
To reproduce the experiments described in the paper, follow the instructions listed below.


### Contacts ###
You can contact us at antonio.cruciani@aalto.fi and leonardo.pellegrina@unipd.it for any questions and for reporting bugs.

### Aknowledgments ###
