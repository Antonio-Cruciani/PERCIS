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

epsilon_list = [0.05,0.01,0.005,0.001,0.0005,0.0001]
ss_save = [1,2,3,4,5,6]
delta = 0.05
runs = 10
graphs_path = "../../../datasets/graphs/"

#sampling_rate_ = 2.3
sampling_rate_ = 0.0
tn = 64
directed = false
output = ""

@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for SILVAN on Labeled Graphs")


datasets = ["25_combined_edges.txt","27_guncontrol_edges.txt"]




global j = 1
for eps in epsilon_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        @info("Input Graph Path: $gf")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"silvan_ss_"*string(j)*"_run_"*string(i)*".txt"
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
            args = `-v 10 -s $sampling_rate_ -o $op $eps $delta $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../silvan/silvan $args`)
                @info("$line")
                _catch_and_update!(line,results)  
                flush(stderr)              
            end
            op_times = "silvan_ss_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end



datasets = ["32_youtube_edges.txt"]




global j = 1
for eps in epsilon_list
    for ds in datasets
        ds_name = string(split(ds,".txt")[1])
        gf = graphs_path*ds
        create_folder(ds_name)
        outpath = "../julia_scripts/scores/"*ds_name*"/"
        check_file_existence(gf)
        @info("Input Graph Path: $gf")
        @info("Running experiements for "*gf)
        for i in 1:runs
            op = outpath *"silvan_ss_"*string(j)*"_run_"*string(i)*".txt"
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
            args = `-d -v 10 -s $sampling_rate_ -o $op $eps $delta $gf`
            #args = `-u -v 1 -g $ss -o $op -t $tn $epsilon $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            for line in eachline(`../silvan/silvan $args`)
                @info("$line")
                _catch_and_update!(line,results)  
                flush(stderr)              
            end
            op_times = "silvan_ss_"*string(j)*".txt"
            #op_times = "uniform_ss_"*string(ss)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end


