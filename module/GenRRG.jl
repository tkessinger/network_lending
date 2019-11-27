#!/usr/bin/env julia

## GenRRG.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Generate connected random regular graphs.
## Heavily borrowed from LightGraphs' random_regular_graph function.

module GenRRG

    using Random, StatsBase

    export generate_RRG, generate_connected_RRG, is_connected

    function is_suitable(
        edges::Array{Array{Int64,1}},
        potential_edges::Dict{Int64,Int64},
        extant_edges::Array{Array{Int64,1}} = [Int64[]]
        )
        # given a Dict of potential edges,
        # see if there's any way to construct a new edge

        # if the Dict is empty, everything is hunky dory
        if isempty(potential_edges)
            return true
        end
        list = keys(potential_edges)
        # if it's not empty, check if a new edge can be created
        for s1 in list, s2 in list
            s1 >= s2 && continue
            # if an edge doesn't already exist, that's good
            if [s1, s2] ∉ edges && [s1, s2] ∉ extant_edges
                return true
            end
        end
        # if there's no way to create a new edge, return false
        return false
    end

    function try_random_regular(
        N::Int64,
        k::Int64,
        extant_edges::Array{Array{Int64,1},1} = Array{Int64,1}[]
        )
        # attempt to generate a random regular graph
        # based on a degree and a list of current edges
        # if this fails, return an empty list: the ``main'' function will try again


        # list of ``unfilled'' edges
        stubs = vcat([fill(i,k) for i in 1:N] ...)
        # the empty edge list of the current type, which we will populate
        edges = Array{Int64,1}[]
        # if the stub list isn't empty, let's try to empty it
        while !isempty(stubs)
            # a Dict() containing the number of possible edges
            # that each stub could be part of
            potential_edges = Dict{Int64,Int64}()
            shuffle!(stubs)
            # consider random pairs of stubs
            for i in 1:2:length(stubs)
                s1, s2 = stubs[i:(i+1)]
                # order the stubs
                if s1 > s2
                    s1, s2 = s2, s1
                end
                # if the edge isn't a loop and doesn't already appear anywhere,
                # then add it to the list of edges
                if s1 != s2 && ∉([s1, s2], edges) && ∉([s1, s2], extant_edges)
                    # println("adding edge [$s1, $s2]")
                    push!(edges, [s1, s2])
                # if this condition fails, add 1 to potential_edges
                else
                    potential_edges[s1] = get(potential_edges, s1, 0) + 1
                    potential_edges[s2] = get(potential_edges, s2, 0) + 1
                end
            end
            # see if there's any way to construct a new edge based on what we've got
            if !(is_suitable(edges, potential_edges, extant_edges))
                # println("not suitable!")
                # if not, return an empty edge list
                return Array{Int64, 1}[]
            end

            # repopulate the stub list based on what remains, then try again
            stubs = Int64[]
            for (e, ct) in potential_edges
                append!(stubs, fill(e, ct))
            end
            #println("$stubs, $edges, $potential_edges")
        end
        return edges
    end

    function generate_RRG(
        N::Int64,
        k::Int64
        )

        num_attempts = 0

        neighbors = [Int64[] for i in 1:N]

        edges = []
        while isempty(edges)
            num_attempts += 1
            edges = try_random_regular(N, k)
        end
        for (s1, s2) in edges
            push!(neighbors[s1], s2)
            push!(neighbors[s2], s1)
        end
        #println("$num_attempts")
        adjacency = falses(N,N)
        for edge in edges
            adjacency[edge[1], edge[2]] = true
            adjacency[edge[2], edge[1]] = true
        end
        graph = (neighbors, edges, adjacency)
        return graph
    end

    function generate_connected_RRG(
        N::Int64, # number of individuals
        k::Int64  # degree of graph
        )
        # generate a random multitype regular graph,
        # then check to see if it's connected
        # if not, keep going until it is
        num_attempts = 1
        let connected = false
            while !(connected)
                global graph = generate_RRG(N, k)
                neighbors, edges = graph
                connected = is_connected(neighbors)
                num_attempts += 1
            end
        end
        # graph[1] is neighbors
        # graph[2] is edges
        # graph[3] is the adjacency matrix
        return graph
    end

    function assemble_connected_list(
        neighbors::Array{Array{Int64,1},1}, # list of neighbors of each individual
        curr_connections::Array{Int64,1}, # list containing individuals that have already been counted
        indv::Int64, # current individual
        )
        # create a list of all individuals the current individual is connected to
        # this is used to determine whether the graph is connected or not

        # if the current individual isn't on the list, put them there
        if indv ∉ curr_connections
            push!(curr_connections, indv)
        end
        # look at all their neighbors
        for neighb in neighbors[indv]
            if neighb ∉ curr_connections
                # add the neighbor to the list
                push!(curr_connections, neighb)
                # add all their connections, and all their connections' connections (recursively)
                new_connections = assemble_connected_list(neighbors, curr_connections, neighb)
                [push!(curr_connections, x) for x in new_connections if x ∉ curr_connections]
            end
        end
        # concatenate the list so ∈ will work
        return vcat(curr_connections ...)
    end

    function is_connected(
        neighbors::Array{Array{Int64,1},1},
        )
        # check to see if our graph is connected
        N = length(neighbors)
        # assemble a list of everyone connected to individual 1
        connected_list = assemble_connected_list(neighbors, Int64[], 1)
        connected_list = sort(unique(connected_list))
        # if it's as long as the entire network, there must be only one connected component
        if length(connected_list) == N
            return true
        else
            return false
        end
    end
end
