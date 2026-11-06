# Source-only rho experiments. Run from your existing julia_scripts directory.
# Executable: ../percis/percis_source
# num_samples counts production bags; void_samples counts unreachable pairs.
# kernel_bulding_time is NaN because it is not reported separately.

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
    dmax = match(r"max_importance_weight=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if dmax !== nothing
        results["d_max"] = parse(Float64, dmax.captures[1])
    end

    # The current source executable does not report kernel construction time.
    # Its column remains NaN unless a future executable supplies this message.
    kern = match(r"Sampling Kernel Built in\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
    if kern !== nothing
        results["kernel_bulding_time"] = parse(Float64, kern.captures[1])
    end

    # Read counts from the final summary, not the pilot's progress messages.
    if startswith(line, "mode=")
        voids = match(r"(?:^|\s)unreachable_bags=(\d+)", line)
        if voids !== nothing
            results["void_samples"] = parse(Float64, voids.captures[1])
        end
        ns = match(r"(?:^|\s)bags=(\d+)", line)
        if ns !== nothing
            results["num_samples"] = parse(Float64, ns.captures[1])
        end
        m = match(r"elapsed_seconds=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", line)
        if m !== nothing
            results["total_time"] = parse(Float64, m.captures[1])
        end
    end
    return nothing
end

function save_results(results,output_path,nn,fn)
    for key in ("total_time", "num_samples", "void_samples", "d_max")
        isfinite(results[key]) || error("Missing statistic: $key")
    end
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
epsilon_list = [0.05,0.01,0.005,0.001,0.0005,0.0001]
ss_save = [1,2,3,4,5,6]
delta = 0.05
#epsilon = 0.1
runs = 10
graphs_path = "../../../datasets/graphs/"
percolation_path = "../../../datasets/percolation_states/"

sampling_rate_ = 2.3
#sampling_rate_ = 0.0
tn = 64
directed = false
output = ""

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]
@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Random Initiator Experiment")

# Undirected
# Non Uniform

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
                flush(stderr)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

# Directed

datasets = ["15_cit_hepph.txt" ,"14_p2p_gnutella31.txt","11_soc_epinions.txt","12_soc_slashdot.txt","04_web_notredame.txt","06_web_google.txt"]

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-d -s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
                flush(stderr)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis_source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

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

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
                flush(stderr)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis_source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

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

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-d -s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
                flush(stderr)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Worst Case Setting")
#graphs_path = "../../../percolation_centrality/components/"

graphs_path = "../../../percolation_centrality/components/"
datasets = ["01_musae_facebook_edges_lcc_in_50.txt","02_email_enron_lcc_in_50.txt","03_ca_astroph_lcc_in_50.txt"]

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
                flush(stderr)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

# Directed
#graphs_path = "../../../percolation_centrality/components/"

global  j=1
datasets = ["15_cit_hepph_lcc_in_50.txt" ,"14_p2p_gnutella31_lcc_in_50.txt","11_soc_epinions_lcc_in_50.txt","12_soc_slashdot_lcc_in_50.txt","04_web_notredame_lcc_in_50.txt","06_web_google_lcc_in_50.txt"]

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-d -s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

# Uniform Percolation States
# Restore the graph directory after the component experiments.
graphs_path = "../../../datasets/graphs/"
@info("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
@info("Running Experiments for Uniform Percolation States")

datasets = ["01_musae_facebook_edges.txt","02_email_enron.txt","03_ca_astroph.txt"]

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
                flush(stderr)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

datasets = ["25_combined_edges.txt","27_guncontrol_edges.txt"]

global j = 1
for eps in epsilon_list
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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

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
            op = outpath *"non_uniform_rho_source_ss_"*string(j)*"_run_"*string(i)*".txt"

            results = Dict(
                "total_time" => NaN,
                "d_max" => NaN,
                "kernel_bulding_time"=> NaN,
                 "void_samples"=> NaN,
                 "time_bfs"=>0.0,
                 "num_samples"=> NaN,
                 "run" => i
            )

            args = `-d -s $sampling_rate_ --seed $i -o $op -t $tn $eps $delta $ps $gf`
            @info("----------------------------------------------------------------------------------")
            @info("Run Number $i")
            # Merge stdout and stderr: source statistics are written to stderr.
            log_pipe = Pipe()
            process = run(pipeline(`../percis/percis-source $args`,
                stdout=log_pipe, stderr=log_pipe); wait=false)
            close(log_pipe.in)
            for line in eachline(log_pipe)
                @info("$line")
                _catch_and_update!(line,results)
            end
            wait(process)
            close(log_pipe)
            success(process) || error("percis-source failed")
            op_times = "non_uniform_rho_source_ss_"*string(j)*".txt"

            save_results(results,"../julia_scripts/",ds_name,op_times)
            @info("Completed")
            flush(stderr)
        end

    end
    global j+=1
end

