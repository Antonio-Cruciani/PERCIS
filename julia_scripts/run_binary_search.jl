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
            write(file,"total_time num_samples kernel_bulding_time void_samples d_max run\n")
        end
    end
    f = open(output_path*"times/" * nn * "/"*fn, "a")
    write(f, string(results["total_time"]) * " " * string(results["num_samples"]) *" "*string(results["kernel_bulding_time"])*" "*string(results["void_samples"])*" "*string(results["d_max"])*" "*string(results["run"])   *"\n")
    close(f)
end

#sample_size_list = [1000,5000,10000,50000,100000,500000,1000000]
#sample_size_list = [1000]
#sample_size_list = [500000,1000000]
sample_size_list = [1000,5000,10000,50000,100000,500000,1000000]

ss_save = [1,2,3,4,5,6,7]
delta = 0.05
epsilon = 0.1
runs = 10

#graphs_path = "../julia_scripts/graphs/"
#percolation_path = "../julia_scripts/percolation_states/"
graphs_path = "../../percolation_centrality/graphs/"
percolation_path = "../../percolation_centrality/percolation_states/"
#graphs_path = "/home/antonio/Desktop/RES_PERCOLATION/EXACT/graphs/"
#percolation_path = "/home/antonio/Desktop/percolation_states/"
tn = 64
directed = false
output = ""

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]
#datasets = ["02_email_enron.txt"]
#datasets = ["10_flickr.txt"]

@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Random Initiator Experiment")
# Undirected
# Non Uniform 

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)  
                flush(stderr)              
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end




# Directed

datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
#datasets = ["08_web_berkstan.txt"]

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-d -v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end




@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Random Spread Experiment")


datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]
#datasets = ["10_flickr.txt"]


for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_e_log"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end


# Directed


datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]
#datasets = ["14_p2p_gnutella31.txt"]

#datasets = ["08_web_berkstan.txt"]



for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_e_log"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-d -v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)   
                flush(stderr)                           
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end


graphs_path = "../../percolation_centrality/components/"


@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Worst Case Setting")

datasets = ["01_musae_facebook_edges_lcc_in_50.txt","02_email_enron_lcc_in_50.txt","03_ca_astroph_lcc_in_50.txt"]
#datasets = ["10_flickr.txt"]

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)         
                flush(stderr)                     
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end





# Directed

datasets = ["15_cit_hepph_lcc_in_50.txt" ,"14_p2p_gnutella31_lcc_in_50.txt","11_soc_epinions_lcc_in_50.txt","12_soc_slashdot_lcc_in_50.txt","04_web_notredame_lcc_in_50.txt","06_web_google_lcc_in_50.txt"]
#datasets = ["08_web_berkstan_lcc_in_50.txt"]

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-d -v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)                
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end



graphs_path = "../../percolation_centrality/graphs/"



@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Uniform Percolation States")
# Undirected
# Non Uniform 
datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_unif"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)  
                flush(stderr)              
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end




# Directed
graphs_path = "../../percolation_centrality/graphs/"
datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]

#datasets = ["04_web_notredame.txt" ]
#datasets = ["08_web_berkstan.txt"]

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_unif"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-d -v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)      
                flush(stderr)                        
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end


datasets = ["24_uselections.txt","23_twitter_pol.txt","22_obamacare.txt","21_brexit.txt","20_abortion.txt"]
@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Real-World Instances Experiment")
# Undirected
# Non Uniform 

for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        check_file_existence(ps)
        @info("Input Graph Path: $gf")
        @info("Input Percolation States Path: $ps")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"non_uniform_bs_ss_"*string(ss)*"_run_"*string(i)*".txt"
            #op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

            #@info("Running Run Number "*string(i))
            results = Dict(
                "total_time" =>0.0,
                "d_max" => 0.0,
                "kernel_bulding_time"=>0.0,
                 "void_samples"=> 0.0,
                 "time_bfs"=>0.0,
                 "num_samples"=>0.0,
                 "run" => i
            )
            
            #println(ps)
            #//output = read(`./aperitif -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`, String)
            args = `-v 10 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)  
                flush(stderr)              
            end
            op_times = "non_uniform_bs_ss_"*string(ss)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end