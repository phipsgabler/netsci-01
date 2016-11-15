module NetSci01

using LightGraphs, GraphPlot
using Distributions

export readnetwork, samplenetwork

const line_regex = r"^(\d+)\s(\d+)"

"Read a space separated file into an (undirected) Graph in an efficient way."
function readnetwork(filename::String, limit::Number = Inf; fromzero::Bool = false)
    graph = Graph()
    vertices = 0
    correction = Int(fromzero)

    # using grep to exclude comment lines
    open(filename) do file
        for (l, line) in enumerate(eachline(file))
            if l > limit
                break
            end

            if (m = match(line_regex, line)) !== nothing
                raw_v1, raw_v2 = m.captures
                v1 = parse(Int, raw_v1) + correction
                v2 = parse(Int, raw_v2) + correction
            
                if v1 > vertices || v2 > vertices
                    @assert add_vertices!(graph, max(v1, v2) - vertices)
                    vertices = max(v1, v2)
                end
                
                # we explicitely ignore inverse directions, if there are any
                @assert has_edge(graph, v1, v2) || add_edge!(graph, v1, v2)
            end
        end
    end

    return graph
end



"Choose a vertex which has not already been visited"
@inline function choosestart(translation::Dict{Int, Int}, vertices::Int)
    next = rand(1:vertices)
    while haskey(translation, next)
        next = rand(1:vertices)
    end
    
    return next
end

"Randomly sample n nodes from a given network"
function samplenetwork(graph::Graph, n::Integer, method::Symbol = :rw; options...)
    sampling_method = @eval $(Symbol("samplenetwork_", method))
    sampling_method(graph, n; options...)
end

"Sample using a random walk"
samplenetwork_rw(graph::Graph, n::Integer; options...) =
    samplenetwork_rj(graph, n; c = 0.0, options...)


"Samples from a network graph using a random walk until n nodes are visited, with
some stochastic jumping"
function samplenetwork_rj(graph::Graph, n::Integer;
                          c::Float64 = 0.15, limit_factor::Integer = 100)
    @assert n >= 0
    @assert 0 <= c <= 1
    @assert limit_factor > 0

    if n == 0
        return Graph()
    end
    
    step = 0
    walk_limit = round(Int, limit_factor / log10(n + 1) * n)
    vertices = nv(graph)
    
    last = rand(1:vertices)
    
    translation = Dict(last => 1)
    sizehint!(translation, n)
    
    node_count = 1
    sample = Graph(n)
    
    while node_count < n
        if rand() > c && !isempty(neighbors(graph, last)) && step < walk_limit
            # choose random neighbour with probability (1-c)
            next = rand(neighbors(graph, last))
            
            if !haskey(translation, next)
                node_count += 1
                translation[next] = node_count
            end
            
            add_edge!(sample, translation[last], translation[next])
            last = next
            step += 1
        else
            # reinitialize; jump, or got stuck in component, or neighbours were empty
            last = choosestart(translation, vertices)
            node_count += 1
            translation[last] = node_count
            step = 0
        end
    end
    
    return sample
end

"Forest fire: samples from a network graph using a randomized bread-first search,
with occasional restarts if necessary"
function samplenetwork_ff(graph::Graph, n::Integer;
                          p::Float64 = 0.15, limit_factor::Integer = 100)
    
#    Leskovec, J. Kleinberg, and C. Faloutsos. Graphs over time: Densification laws, shrinking
#    diamaters and possible explanations. In ACM SIGKDD , 2005.
    
    @assert n >= 0
    @assert 0 <= p <= 1
    @assert limit_factor > 0

    if n == 0
        return Graph()
    end
    
    vertices = nv(graph)
    
    translation = Dict{Int, Int}()
    sizehint!(translation, n)
    
    node_count = 0
    sample = Graph(n)
    
    # we can't add too many nodes, because the number of nodes of `sample` is fixed
    while node_count < n
        initial = choosestart(translation, vertices)
        translation[initial] = (node_count += 1)

        burning = [initial]
        finished = false
        
        while !isempty(burning)
            ignited = Int[]
            
            for b in burning
                for n in neighbors(graph, b)
                    if rand() > p
                        continue
                    elseif !haskey(translation, n)
                        translation[n] = (node_count += 1)
                    end
                    
                    if !add_edge!(sample, translation[b], translation[n])
                        finished = true
                        break
                    end
                    push!(ignited, n)
                end

                if finished
                    break
                end
            end
            
            burning = ignited
        end        
    end
    
    return sample
end


end
