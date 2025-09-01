function read_percolation(file_name::String)
    @assert isfile(file_name) "The percolation values file does not exist"
    f = open(file_name, "r")
    p::Array{Float64} = []
    for line in eachline(f)
        l = parse(Float64, line)
        if l > 1.0 || l < 0.0
            println("Error, percolation state "*string(l)*" is not a valid state!")
            return nothing
        end 
        push!(p,parse(Float64, line))
    end
    close(f)
    return p
end


function save_results_new(nn::String, cn::String, c::Array{Float64}, t::Float64,dv::Float64)
    if (length(c) > 0)
        mkpath("scores/" * nn * "/")
        save_centrality_values("scores/" * nn * "/" * cn * ".txt", c)
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, string(t))
        close(f)
        mkpath("dv/"*nn*"/")
        f = open("dv/" * nn * "/dv_" * cn * ".txt", "w")
        write(f, string(dv))
        close(f)
    else
        mkpath("times/" * nn * "/")
        f = open("times/" * nn * "/time_" * cn * ".txt", "w")
        write(f, "-1.0\n")
        close(f)
    end
end



function save_centrality_values(file_name::String, centrality::Array{Float64})::Nothing
    f::IOStream = open(file_name, "w")
    for u in 1:lastindex(centrality)
        write(f, string(centrality[u]) * "\n")
    end
    close(f)
end


