function create_folder(nn)
    mkpath("scores/" * nn * "/")
    mkpath("times/" * nn * "/")
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

sample_size_list = [1000,5000,10000,50000,100000]
#sample_size_list = [1000]

ss_save = [1,2,3,4,5]
delta = 0.05
epsilon = 0.1
runs = 1
graphs_path = "../julia_scripts/graphs/"
percolation_path = "../julia_scripts/percolation_states/"
tn = 20
directed = false
output = ""
datasets = ["02_email_enron.txt"]
# Non Uniform 
for ss in sample_size_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])*"_rnd_init_50"
        gf = graphs_path*ds
        create_folder(ds_name)
        ps = percolation_path*ds_name*".txt"
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        
        @info("Running experiements for "*gf)
        for i in 1:runs
            #op = outpath *"non_uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"
            op = outpath *"uniform_ss_"*string(ss)*"_run_"*string(i)*".txt"

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
            #args = `-v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../aperitif/aperitif $args`)
                @info("$line")
                _catch_and_update!(line,results)                
            end
            #op_times = "non_uniform_ss_"*string(ss)*".txt"
            op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
end