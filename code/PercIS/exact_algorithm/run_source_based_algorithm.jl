include("header.jl")



graphs_path = "../../../datasets/graphs/"
percolation_path = "../../../datasets/percolation_states/"
output_path = ""
# Insert the name of the dataset header
# example: datasets = ["name.txt"]

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]
# If undirected otherwise set it to true
directed = false

perc_state_list = ["_unif","_rnd_init_50","_e_log"]

# File separator
separator = "\t"

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1]) * ps
    for ps in perc_state_list
        x = read_percolation(percolation_path*ds_name*".txt")
        @info("Computing Ground Truth percolation centrality for "*ds_name)
        flush(stderr)
        y = parallel_percolation_centrality_source_based(g,x)
        save_results_new(ds_name,"exact_source_based",y[1],y[5],y[4])
    end
end


# Only for _lcc_in_50
datasets = ["01_musae_facebook_edges_lcc_in_50.txt","02_email_enron_lcc_in_50.txt","03_ca_astroph_lcc_in_50.txt"]

for ds in datasets
    gf = graphs_path*ds
    g = load_graph(gf,directed,separator)
    ds_name = string(split(ds,".txt")[1])
    x = read_percolation(percolation_path*ds_name*".txt")
    @info("Computing Ground Truth percolation centrality for "*ds_name)
    flush(stderr)
    y = parallel_percolation_centrality_source_based(g,x)
    save_results_new(ds_name,"exact_source_based",y[1],y[5],y[4])
    
end