function create_folder(nn)
    mkpath("scores/" * nn * "/")
    mkpath("times/" * nn * "/")
end

function check_file_existence(fp)
    if !isfile(fp)
        @info("File $fp does not exist!!!")
        exit(1)
    end
end
function _catch_and_update!(line,results)
    dmax = match(r"d_max:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if dmax!== nothing
        d_max = parse(Float64, dmax.captures[1])
        results["d_max"] = d_max
    end

    kern = match(r"Sampling Kernel Built in\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if kern !== nothing
        kernel_bulding_time = parse(Float64, kern.captures[1])
        results["kernel_bulding_time"] = kernel_bulding_time
    end
    voids = match(r"void_samples:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if voids !== nothing
        void_samples = parse(Float64, voids.captures[1])
        results["void_samples"] = void_samples
    end
    ns = match(r"num_samples:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if ns !== nothing
        num_samples = parse(Float64, ns.captures[1])
        results["num_samples"] = num_samples
    end
    m = match(r"Total time:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if m !== nothing
        total_time = parse(Float64, m.captures[1])
        results["total_time"] = total_time
    end
    m = match(r"Final Supremum Deviation:\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if m !== nothing
        final_sd = parse(Float64, m.captures[1])
        results["sd"] = final_sd
    end
    return nothing
end

function save_results(results,output_path,nn,fn)
    #mkpath("results_greedy/"*graph_name*"/")
    header = true
    if isfile(output_path*"times/"*nn*"/"*fn)
        header = false
    end
    open(output_path*"times/"*nn*"/"*fn,"a") do file
        if header
            write(file,"total_time num_samples kernel_bulding_time void_samples d_max sd epsilon run\n")
        end
    end
    f = open(output_path*"times/" * nn * "/"*fn, "a")
    write(f, string(results["total_time"]) * " " * string(results["num_samples"]) *" "*string(results["kernel_bulding_time"])*" "*string(results["void_samples"])*" "*string(results["d_max"])*" "*string(results["sd"])*" "*string(results["epsilon"])*" "*string(results["run"])   *"\n")
    close(f)
end

function read_centrality_values(file_name::String)::Array{Float64}
    @assert isfile(file_name) "The centrality value file does not exist "*file_name
    f::IOStream = open(file_name, "r")
    centrality::Array{Float64} = []
    value::Float64 = 0.0
    for l in eachline(f)
        value = parse(Float64, l)
        if (value < -0.1)
            println("ERROR. There are negative values with absolute big values")
            return Array{Float64}([])
        end
        if (value < 0)
            value = 0
        end
        push!(centrality, value)
    end
    close(f)
    return centrality
end


sample_size_list = [1000,5000,10000,50000,100000]
#sample_size_list = [1000]
epsilon_list = [0.05,0.01,0.005,0.001,0.0005]
ss_save = [1,2,3,4,5]
delta = 0.05
sampling_window = 1000
#epsilon = 0.1
runs = 10
#graphs_path = "../../percolation_centrality_via_non_uniform_sampling/julia_scripts/graphs/"
#percolation_path = "../../percolation_centrality_via_non_uniform_sampling/julia_scripts/percolation_states/"
#exact_scores_path = "../../percolation_centrality/scores/"
epsilon_list = [0.0003]
#graphs_path = "/home/antonio/Desktop/RES_PERCOLATION/EXACT/graphs/"
#percolation_path = "/home/antonio/Desktop/RES_PERCOLATION/EXACT/percolation_states/"
#exact_scores_path = "/home/antonio/Desktop/RES_PERCOLATION/EXACT/scores/"
graphs_path = "../../percolation_centrality/graphs/"
percolation_path = "../../percolation_centrality/percolation_states/"
exact_scores_path = "../../percolation_centrality/scores/"
tn = 64
directed = false
output = ""

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]

#datasets = ["10_flickr.txt"]

@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Supremum Deviation Experiment")
@info("Running Experiments for Random Initiator Experiment")
# Undirected
# Non Uniform 
global j = 1
fractional = [2,4,5,10]
for denom in fractional
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)  
                flush(stderr)              
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end


# Directed

datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
#datasets = ["08_web_berkstan.txt"]

global j = 1
for denom in fractional
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"

            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -d -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Random Spread Experiment")


datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]
#datasets = ["10_flickr.txt"]

global j = 1

for denom in fractional
     for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target_e_log.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end





# Directed

global  j = 1
datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
#datasets = ["08_web_berkstan.txt"]




for denom in fractional
     for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target_e_log.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -d -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end
    end
    global j+=1
end

@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Uniform Distributed Percolation States Experiment")

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]
#datasets = ["10_flickr.txt"]

global j = 1

for denom in fractional
     for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_unif"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end





# Directed

global  j = 1
datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
#datasets = ["08_web_berkstan.txt"]




for denom in fractional
     for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_unif"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -d -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end
    end
    global j+=1
end

#=
@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Worst Case Experiment")







datasets = ["02_email_enron_lcc_in_50.txt","03_ca_astroph_lcc_in_50.txt"]
global j = 1


for eps in epsilon_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end
=#

global  j=1 
#datasets = ["15_cit_hepph_lcc_in_50.txt" ,"14_p2p_gnutella31_lcc_in_50.txt","11_soc_epinions_lcc_in_50.txt","12_soc_slashdot_lcc_in_50.txt","04_web_notredame_lcc_in_50.txt","06_web_google_lcc_in_50.txt"]
#=
datasets = ["14_p2p_gnutella31_lcc_in_50.txt"]
for eps in epsilon_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -d -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end
=#
#=
global  j=1 
datasets = ["14_p2p_gnutella31_lcc_in_50.txt","11_soc_epinions_lcc_in_50.txt","15_cit_hepph_lcc_in_50.txt","12_soc_slashdot_lcc_in_50.txt","04_web_notredame_lcc_in_50.txt","06_web_google_lcc_in_50.txt"]
fractional = [2,4,5,10]
for denom in fractional
    for ds in datasets
            ds_name = string(split(ds,".txt")[1])
            gf = graphs_path*ds
            create_folder(ds_name)
            ps = percolation_path*ds_name*".txt"
            outpath = "../julia_scripts/scores/"*ds_name*"/"
            es = exact_scores_path*ds_name*"/exact_target.txt"
            check_file_existence(gf)
            check_file_existence(ps)
            check_file_existence(es)
            @info("Input Graph Path: $gf")
            @info("Input Percolation States Path: $ps")
            @info("Input Exact Values Path: $es")
            @info("Running experiements for "*gf)
            @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
            perc_cent = read_centrality_values(es)
            mp = maximum(perc_cent)
            eps = mp/denom
            @info("Maximum Exact Percolation = $mp , Target ε = $eps")
            for i in 1:runs
                op = outpath *"non_uniform_bs_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
                #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

                #@info("Running Run Number "*string(i))
                results = Dict(
                    "total_time" =>0.0,
                    "d_max" => 0.0,
                    "kernel_bulding_time"=>0.0,
                    "void_samples"=> 0.0,
                    "time_bfs"=>0.0,
                    "num_samples"=>0.0,
                    "sd" => 0.0,
                    "run" => i,
                    "epsilon"=>eps
                )
                
                
                #println(ps)
                #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
                args = `-v 10 -d -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
                #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
                @info("----------------------------------------------------------------------------------")
                @info("Run Number $i")
                for line in eachline(`../aperitif/aperitif $args`)
                    @info("$line")
                    _catch_and_update!(line,results)      
                    flush(stderr)                        
                end
                op_times = "non_uniform_bs_sd_"*string(j)*".txt"
                #op_times = "uniform_ss_"*string(ss)*".txt"

                save_results(results,"../julia_scripts/",ds_name,op_times)
                @info("Completed")
                flush(stderr)
            end

    end
    global j+=1
end
global  j=1 

#for eps in epsilon_list
for denom in fractional

    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        es = exact_scores_path*ds_name*"/exact_target.txt"
        check_file_existence(gf)
        check_file_existence(ps)
        check_file_existence(es)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Input Exact Values Path: $es")
        @info("Running experiements for "*gf)
        @info("Reading Centrality and setting ε = (1/2)⋅max p(v)")
        perc_cent = read_centrality_values(es)
        mp = maximum(perc_cent)
        eps = mp/denom
        @info("Maximum Exact Percolation = $mp , Target ε = $eps")
        for i in 1:runs
            op = outpath *"uniform_SD_ss_"*string(j)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "sd" => 0.0,
                 "run" => i,
                 "epsilon"=>eps
            )
            
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -d -u -o $op -e $es -w $sampling_window -t $tn $eps $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "uniform_sd_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
#end
end

=#