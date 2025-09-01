
struct static_graph
    adjacency::Array{Array{Int64}}
    incidency::Array{Array{Int64}}
    degrees_adj::Array{Int64}
    degrees_idj::Array{Int64}
    function static_graph(adj::Array{Array{Int64}},idj::Array{Array{Int64}})
        return new(adj,idj,[lastindex(adj[u]) for u in 1:lastindex(adj)],[lastindex(idj[u]) for u in 1:lastindex(idj)])
    end
end


function load_graph(file_name::String,directed::Bool = true,sep::String = " ")
    @assert isfile(file_name) "The edge list file " * file_name * " does not exist"
    start_time = time()
    file_id_to_graph_id::Dict{Int64,Int64} = Dict{Int64,Int64}()
    current_node_id::Int64 = 1
    f::IOStream = open(file_name, "r")
	edges::Set{Tuple{Int64,Int64}} = Set{Tuple{Int64,Int64}}()
    nodes::Set{Int64} = Set{Int64}()
    for line in eachline(f)
        split_line::Vector{String} = split(line, sep)
        @assert length(split_line) == 2 "Bad line format: " * line
        u = parse(Int64,split_line[1])
		v = parse(Int64,split_line[2])
        if u != v
            if (!haskey(file_id_to_graph_id, u))
                file_id_to_graph_id[u] = current_node_id
                #push!(file_id, u)
                current_node_id = current_node_id + 1
            end
            if (!haskey(file_id_to_graph_id, v))
                file_id_to_graph_id[v] = current_node_id
                #push!(file_id, v)
                current_node_id = current_node_id + 1
            end 
            
            if directed
                push!(edges,(u,v))
            else
                push!(edges,(min(u,v),max(u,v)))
            end
            push!(nodes,u)
            push!(nodes,v)
        end
    end
    close(f)

    g = SimpleGraph()
    #n_vertices::Int64 = maximum(collect(keys(file_id_to_graph_id)))
    n_vertices::Int64 = maximum(collect(values(file_id_to_graph_id)))
    if directed
        g = SimpleDiGraph(n_vertices)
    else
        g = SimpleGraph(n_vertices)
    end
    for e in collect(edges)
        res = add_edge!(g,file_id_to_graph_id[e[1]],file_id_to_graph_id[e[2]])
        if !res 
			@info("Edge not added! "*string(e[1])*" , "*string(e[2]))
			@info("Edge not added! "*string(file_id_to_graph_id[e[1]])*" , "*string(file_id_to_graph_id[e[2]]))
		end
    end
    loading_time = time() - start_time
    @info("Loaded Graph ")
	@info("#Nodes "*string(nv(g)))
	@info("#Edges "*string(ne(g)))
	@info("Directed ? "*string(directed))
    @info("Loading time "*string(loading_time)*" seconds")
    flush(stderr)
	return g
end


function print_stats(g; graph_name="anonymous")
    @info("====================================================")
    @info("Network: " * graph_name)
    @info("====================================================")
    @info("Number of nodes " * string(nv(g)))
    @info("Number of edges " * string(length(ne(g))))
    @info("Is the graph directed? " * string(is_directed(g)))
    @info("====================================================")
    flush(stderr)
end


function adjacency_list(g)::Array{Array{Int64}}
    adj_list::Array{Array{Int64}} = Array{Array{Int64}}([])
    for u in 1:nv(g)
        push!(adj_list,outneighbors(g,u))
    end
    return adj_list
end

function incidency_list(g)::Array{Array{Int64}}
    in_list::Array{Array{Int64}} = Array{Array{Int64}}([])
    for u in 1:nv(g)
        push!(in_list,inneighbors(g,u))
    end
    return in_list
end


function save_graph(g,nn,sep::String = " ")
    mkpath("components/")
    f = open("components/" * nn * ".txt", "w")
    for e in collect(edges(g))
        write(f,string(src(e))*sep*string(dst(e))*"\n")
    end
    close(f)
end


function load_and_normalize_edgelist(path::String)
    # Containers
    edges = Tuple{String, String}[]
    label_to_id = Dict{String, Int}()
    id_to_label = String[]
    
    # Read file
    open(path, "r") do file
        for line in eachline(file)
            if isempty(strip(line)) || startswith(line, "#")
                continue  # skip comments and empty lines
            end
            tokens = split(strip(line))
            if length(tokens) < 2
                continue  # skip malformed lines
            end
            push!(edges, (tokens[1], tokens[2]))
        end
    end

    # Map labels to consecutive IDs
    next_id = 0
    for (u, v) in edges
        for node in (u, v)
            if !haskey(label_to_id, node)
                label_to_id[node] = next_id
                push!(id_to_label, node)
                next_id += 1
            end
        end
    end

    # Normalize edges
    normalized_edges = [(label_to_id[u], label_to_id[v]) for (u, v) in edges]

    return normalized_edges, label_to_id, id_to_label
end


function load_and_normalize_edgelist_percolation_states(path::String,path_percs::String)
    # Containers
    edges = Tuple{String, String}[]
    label_to_id = Dict{String, Int}()
    id_to_label = String[]
    
    # Read file
    open(path, "r") do file
        for line in eachline(file)
            if isempty(strip(line)) || startswith(line, "#")
                continue  # skip comments and empty lines
            end
            tokens = split(strip(line))
            if length(tokens) < 2
                continue  # skip malformed lines
            end
            push!(edges, (tokens[1], tokens[2]))
        end
    end

    # Map labels to consecutive IDs
    next_id = 0
    for (u, v) in edges
        for node in (u, v)
            if !haskey(label_to_id, node)
                label_to_id[node] = next_id
                push!(id_to_label, node)
                next_id += 1
            end
        end
    end

    # Normalize edges
    normalized_edges = [(label_to_id[u], label_to_id[v]) for (u, v) in edges]

    normalized_states::Array{Float64} = zeros(Float64,next_id)
    # Read file
    open(path_percs, "r") do file
        for line in eachline(file)
            if isempty(strip(line)) || startswith(line, "#")
                continue  # skip comments and empty lines
            end
            tokens = split(strip(line))
            if length(tokens) < 2
                continue  # skip malformed lines
            end
            normalized_states[label_to_id[tokens[1]]+1] = parse(Float64,tokens[2])
        end
    end

    return normalized_edges,normalized_states ,label_to_id, id_to_label
end

function save_edge_list(el,nn,sep = " ")
    mkpath("normalized/")
    f = open("normalized/" * nn * ".txt", "w")
    for e in el
        write(f,string(e[1])*sep*string(e[2])*"\n")
    end
    close(f)

end


@inline function compute_d_max(n::Int64,X::Array{Float64})
    start_time::Float64 = time()
    @info("Computing d_v = max_{sv} κ/̃κ")
    tmp_perc_states::Dict{Int64,Float64} = Dict(v => X[v] for v in 1:n)
    sorted_dict = OrderedDict(sort(collect(tmp_perc_states), by = kv -> kv[2]))
    percolation_data::Tuple{Float64,Dict{Int64,Float64}} = percolation_differences(sorted_dict,n)



    d::Array{Float64} = ones(Float64,n)

    sorted_indices::Vector{Int64} = sortperm(X)            # Indices that would sort X
    #println("TIPO  ",typeof(sorted_indices))
    sorted_X::Array{Float64} = sort(X)                      # Sort X values
    prefix_sum::Array{Float64} = zeros(Float64, n + 1)      # Prefix sum array

    # Compute prefix sums
    for i in 1:n
        prefix_sum[i + 1] = prefix_sum[i] + sorted_X[i]
    end

    Y::Array{Float64} = zeros(Float64, n)                   # Result array

    # Compute Y values for sorted_X
    for i in 1:n
        left_sum = (i - 1) * sorted_X[i] - prefix_sum[i]
        right_sum = (prefix_sum[end] - prefix_sum[i]) - (n - i) * sorted_X[i]
        Y[sorted_indices[i]] = left_sum + right_sum
    end
   
    for v in 1:n
        d[v] += Y[v]/percolation_data[2][v]
        if Y[v]/percolation_data[2][v] < 0 
            println("Negative "*string(Y[v])*"/"*string(percolation_data[2][v])*" = "*string(Y[v]/percolation_data[2][v]))
        end
    end
    filtered_d::Array{Float64} = Array{Float64}([])
    for v in 1:n
        if !isnan(d[v])
            push!(filtered_d,d[v])
        end
    end

    end_time::Float64 = time()-start_time
    @info("minimum d_v = "*string(minimum(filtered_d))*" maximum d_v = "*string(maximum(filtered_d))*" computed in "*string(end_time)*" seconds")
    return filtered_d,end_time,Y,percolation_data
end


function percolation_differences(percolation_states,n::Int64)::Tuple{Float64,Dict{Int64,Float64}}
    sum::Float64 = 0.0
    minus_sum::Dict{Int64,Float64} = Dict(v => 0.0 for v in 1:n)
    svp::Array{Float64} = zeros(Float64, n + 1)

    j::Int64 = 0
    k::Float64 = 0.0
    for (key, value) in percolation_states
        if j > 0
            svp[j + 1] = svp[j] + k
            sum += j * value - svp[j + 1]
        end
        k = value
        j += 1
    end
    svp[end] = svp[n] + k

    j = 0
    for (key, value) in percolation_states
        minus_sum[key] = sum - value * (2 * j - n) - svp[end] + 2 * svp[j + 1]
        j += 1
    end

    return sum,minus_sum
end

function ramp(x::Float64,y::Float64)::Float64
    if x-y>0
        return (x-y)
    else
        return 0
    end
end