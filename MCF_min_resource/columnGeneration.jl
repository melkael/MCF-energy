include("labelCorrectingForbiddenExperimental.jl")
using Suppressor
using JLD2

function loop_erasure(p)
    i = 1
    out = []
    while i <= length(p)
        j = i+1
        while j <= length(p)
            if p[i] == p[j]
                i = j
            end
            j = j+1
        end
        push!(out, p[i]) 
        i += 1
     end
    out
end

function findImprovingPathsBrum(problem::Problem)
    paths_to_add = Dict()
    for k in problem.commodities
        sn_with_dual_costs = makeBrumGraphEnergy(problem, k)

        forbidden_paths = []#problem.paths[k.id]
        forbidden_paths = edgePathToNodePath.(values(forbidden_paths), k.source)
        
        paths_to_add[k.id] = nothing
        paths, delays = solvePricingSingleObjective(sn_with_dual_costs, k.source, k.meta_destination, [], k.dests)

        for i in 1:length(paths)
            p = p = nodePathToEdgePath(paths[i])

            if delays[i] * 2 < k.max_delay
                reduced_cost = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
                if reduced_cost < 0 && !isapprox(reduced_cost, 0.0, atol=1e-8) && !(p ∈ problem.paths[k.id]) #&& !(paths[i] in paths_seen)# && !has_already_been_added(p, original_paths[k.id])
                    if paths_to_add[k.id] == nothing || reduced_cost < paths_to_add[k.id]["reduced cost"]
                        paths_to_add[k.id] = Dict("reduced cost" => reduced_cost, "path" => p)
                            #break
                    end
                end
            end
        end
        
        if paths_to_add[k.id] == nothing
            paths, delays = solvePricingExperimental(sn_with_dual_costs, k.source, k.meta_destination, forbidden_paths, k.dests)

            for i in 1:length(paths)
                p = nodePathToEdgePath(paths[i])
                
                if delays[i] * 2 < k.max_delay
                    reduced_cost = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
                    
                    if reduced_cost < 0 && !isapprox(reduced_cost, 0.0, atol=1e-8) && !(p ∈ problem.paths[k.id]) #&& !(paths[i] in paths_seen)# && !has_already_been_added(p, original_paths[k.id])
                        if paths_to_add[k.id] == nothing || reduced_cost < paths_to_add[k.id]["reduced cost"]
                            paths_to_add[k.id] = Dict("reduced cost" => reduced_cost, "path" => p)
                            #break
                        end
                    end
                end
            end
        end
    end
    paths_to_add
end

function makeBrumGraphEnergy(problem::Problem, commodity::Commodity)
    sn_modified = deepcopy(problem.sn)

    for e in edges(sn_modified)
        if isForbidden(e, problem.forbidden_arcs[commodity.id])
            set_prop!(sn_modified, e, :wTimesBW, 0)
            set_prop!(sn_modified, e, :propagationDelay, 9999999)
            #println("here")
            #stop()
        else
            BW = commodity.BW_demanded
            # set w_ij
            if get_prop(sn_modified, e, :BW_max) < 100000
                #println("here")
                set_prop!(sn_modified, e, :wTimesBW, ((-shadow_price(problem.c_edge[e])) * BW +BW))
            else
                set_prop!(sn_modified, e, :wTimesBW, ((-shadow_price(problem.c_edge[e])) * BW))
            end
            #set_prop!(sn_modified, e, :sumDualCuts, 0)
            # time is simply propagation delay
            #set_prop!(sn_modified, e, :propagationDelay, get_prop(sn_modified, e, :weight))
        end
    end

    # we add abs(smallest edge) +1 to all edges in order to ensure there are no negative arcs (and hence no negative cycle)
    # the +1 ensures we do not have edges that have (0, 0) costs, which would cause issues when retrieving paths (see retrieve_paths function in labelCorrectingSP.jl)
    min_value = Inf
    for e in edges(sn_modified)
        #if get_prop(sn_modified, e, :wTimesBW) + get_prop(sn_modified, e, :sumDualCuts) < min_value 
        if get_prop(sn_modified, e, :wTimesBW) < min_value 
            min_value = get_prop(sn_modified, e, :wTimesBW)# + get_prop(sn_modified, e, :sumDualCuts)
        end
    end
    set_prop!(sn_modified, :normalization_wTimesBW, abs(min_value)+1)
    for e in edges(sn_modified)
        set_prop!(sn_modified, e, :wTimesBW, get_prop(sn_modified, e, :wTimesBW) + abs(min_value)+1)
    end

    for n in 1:nv(problem.sn)
        for n2 in inneighbors(problem.sn, n)
            set_prop!(sn_modified, n2, n, :sumWandH, get_prop(sn_modified, n2, n, :wTimesBW))
        end
    end
    return sn_modified
end

function getReducedCost(sn::MetaDiGraph, p::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}, commodity::Commodity, c_edge, c_commodity, c_node, forbidden_arcs, knapsackCutsNodes)
    #println("bad version")
    s_wij = 0
    for e in p
        # even when setting costs to large values, we might still have a path with forbidden arc returned here, if its the only possible path in graph
        # in this case we return infinity
        if e in forbidden_arcs
            return Inf
        end
        #BW = commodity.BW_demanded
        #s_wij += BW * -shadow_price(c_edge[e])
        #println(e, " ", get_prop(sn, e, :wTimesBW)- get_prop(sn, :normalization_wTimesBW))
        s_wij += get_prop(sn, e, :wTimesBW) - get_prop(sn, :normalization_wTimesBW)
    end
    #println()

    bw = length_non_dummy(sn, p) * commodity.BW_demanded
    #s_wij -= commodity.BW_demanded

    σ = shadow_price(c_commodity[commodity.id])

    cpu = getCpu(sn, commodity, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p))

    # we look at the SOURCE of last edge (and not destination) because the last node is the meta-destination which we do not care about
    z_end = shadow_price(c_node[p[end].src])
    for cut in knapsackCutsNodes[p[end].src]
        z_end += shadow_price(cut)
    end
    #@show (s_wij)
    #@show σ
    #@show cpu * (z_end)
    #@show cpu
    #@show commodity.BW_demanded
    # this reduced cost is only valid if the path is NOT in the RMP
    rc = -s_wij + σ - cpu * (z_end-1)
    #rc = cpu + -s_wij + σ - cpu * (z_end)
    #@show rc
    #println()

    return rc
end


function columnGeneration(problem::Problem, improvingPaths::Dict{Any, Any}, tree)
    objectives = []
    while !noPathsToAdd(improvingPaths)
        for K_id in keys(improvingPaths)
            if improvingPaths[K_id] != nothing 
                p = improvingPaths[K_id]["path"]
                com = retrieve_commodity(problem.commodities, K_id)
                v = @variable(problem.model, lower_bound=0)
                cpu = getCpu(problem.sn, com, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p))
                push!(problem.cpus[K_id], cpu)

                # add the new path to path list
                push!(problem.paths[K_id], p)
                # add corresponding variable to vars list
                push!(problem.vars[K_id], v)

                problem.map_vars_paths[problem.vars[K_id][end]] = p

                final_node = p[end-1].dst

                # update constraints (2)
                for e in p
                    set_normalized_coefficient(problem.c_edge[e], v, com.BW_demanded)
                    push!(problem.edge_contains_paths[e], (v, com.BW_demanded))
                end
                # update constraints (3)
                set_normalized_coefficient(problem.c_commodity[K_id], v, 1.0)
                # update constraints (4)
                set_normalized_coefficient(problem.c_node[final_node], v, cpu)

                # update objective
                set_objective_coefficient(problem.model, v, com.BW_demanded * length_non_dummy(problem.sn, p) + cpu)
            end
        end

        @suppress optimize!(problem.model)

        #println(objective_value(problem.model))
        #println(onlyHasIntegerPathVars(problem))
        push!(objectives, objective_value(problem.model))
        #println(objective_value(problem.model))
        if length(objectives) > 100 && length(unique(objectives[end-100:end])) == 1
            tree.column_generation_has_broken = true
            println("breakin")
            break
        end
        #println("finding improving paths")
        improvingPaths = findImprovingPathsBrum(problem)
        #println("oopsie")
    end
    return problem.vars, problem.map_vars_paths, problem.paths
end
