include("header.jl")



graphs_path = "../../../datasets/graphs/"
percolation_path = "../../../datasets/percolation_states/"
output_path = ""
# Insert the name of the dataset header
# example: datasets = ["name.txt"]

datasets = ["01_musae_facebook_edges_e_log.txt"]
# If undirected otherwise set it to true
directed = false


# File separator
separator = "\t"

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])

    x = read_percolation(percolation_path*ds_name*".txt")

    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_new_target(g,x)
    save_results_new(ds_name,"exact",y[1],y[6],y[5])
end