include("labelCorrectingForbiddenExperimental.jl")
include("utils.jl")
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

        forbidden_paths = problem.paths[k.id]
        forbidden_paths = edgePathToNodePath.(values(forbidden_paths), k.source)
        
        paths_to_add[k.id] = nothing
        paths, delays = solvePricingSingleObjective(sn_with_dual_costs, k.source, k.meta_destination, [], k.dests)

        for i in 1:length(paths)
            p = p = nodePathToEdgePath(paths[i])

            if delays[i] * 2 < k.max_delay
                reduced_cost = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.c_energy_l, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
                if reduced_cost < 0 && !isapprox(reduced_cost, 0.0, atol=1e-8) && !(p ∈ problem.paths[k.id]) #&& !(paths[i] in paths_seen)# && !has_already_been_added(p, original_paths[k.id])
                    if paths_to_add[k.id] == nothing || reduced_cost < paths_to_add[k.id]["reduced cost"]
                        paths_to_add[k.id] = Dict("reduced cost" => reduced_cost, "path" => p)
                            #break
                    end
                end
            end
        end
        # try without forbidden paths
        if paths_to_add[k.id] == nothing
            paths, delays = solvePricingExperimental(sn_with_dual_costs, k.source, k.meta_destination, [], k.dests)
            for i in 1:length(paths)
                p = p = nodePathToEdgePath(paths[i])
                
                if delays[i] * 2 < k.max_delay
                    reduced_cost = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.c_energy_l, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
                    if reduced_cost < 0 && !isapprox(reduced_cost, 0.0, atol=1e-8) && !(p ∈ problem.paths[k.id]) #&& !(paths[i] in paths_seen)# && !has_already_been_added(p, original_paths[k.id])
                        if paths_to_add[k.id] == nothing || reduced_cost < paths_to_add[k.id]["reduced cost"]
                            paths_to_add[k.id] = Dict("reduced cost" => reduced_cost, "path" => p)
                            #break
                        end
                    end
                end
            end
        end

        # if it didn't work, try with forbidden paths
        if paths_to_add[k.id] == nothing
            paths, delays = solvePricingExperimental(sn_with_dual_costs, k.source, k.meta_destination, forbidden_paths, k.dests)
            for i in 1:length(paths)
                p = p = nodePathToEdgePath(paths[i])
                
                if delays[i] * 2 < k.max_delay
                    reduced_cost = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.c_energy_l, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
                    if reduced_cost < 0 && !isapprox(reduced_cost, 0.0, atol=1e-8)  #&& !(paths[i] in paths_seen)# && !has_already_been_added(p, original_paths[k.id])
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
            set_prop!(sn_modified, e, :wTimesBW, -shadow_price(problem.c_edge[e]) * BW)
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

    min_value = Inf
    # we set energy-related dual variable in nodes
    # note we do not set cut dual variables for nodes because shortest paths are computed with respect to the destination node
    # hence it is equivalent to adding a constant to all paths with the same destination (which does not modify the shortest path we find)
    for n in keys(problem.vars_nodes_a)  #1:nv(sn_modified)
        if shadow_price(problem.c_energy_normalized_a[n]) * 1/length(problem.commodities) < min_value
            min_value = shadow_price(problem.c_energy_normalized_a[n]) * 1/length(problem.commodities)
        end
    end

    set_prop!(sn_modified, :normalization_h, abs(min_value)+1)
    # then we update the value to account for normalization
    for n in 1:nv(problem.sn)
        if n in keys(problem.c_energy_normalized_a)
            set_prop!(sn_modified, n, :h, shadow_price(problem.c_energy_normalized_a[n]) * 1/length(problem.commodities) + abs(min_value)+1)
        else
            set_prop!(sn_modified, n, :h, 0)
        end
    end

    for n in 1:nv(problem.sn)
        for n2 in inneighbors(problem.sn, n)
            set_prop!(sn_modified, n2, n, :sumWandH, get_prop(sn_modified, n2, n, :wTimesBW) + get_prop(sn_modified, n, :h))
        end
    end
    return sn_modified
end

function getReducedCost(sn::MetaDiGraph, p::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}, commodity::Commodity, c_edge, c_commodity, c_node, c_energy_l, forbidden_arcs, knapsackCutsNodes)
    s_wij = 0
    for e in p
        # even when setting costs to large values, we might still have a path with forbidden arc returned here, if its the only possible path in graph
        # in this case we return infinity
        if e in forbidden_arcs
            return Inf
        end
        BW = commodity.BW_demanded
        s_wij += get_prop(sn, e, :wTimesBW) - get_prop(sn, :normalization_wTimesBW)
    end

    s_h = 0
    for e in p
        s_h += get_prop(sn, e.dst, :h) - get_prop(sn, :normalization_h)
    end
    σ = shadow_price(c_commodity[commodity.id])

    cpu = getCpu(sn, commodity, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p))

    # we look at the SOURCE of last edge (and not destination) because the last node is the meta-destination which we do not care about
    z_end = shadow_price(c_node[p[end].src])
    for cut in knapsackCutsNodes[p[end].src]
        z_end += shadow_price(cut)
    end


    q_end = shadow_price(c_energy_l[p[end].src])

    lamda = 1/mean(commodity.arr_mean)
    cpurq = commodity.cpu_request

    # this reduced cost is only valid if the path is NOT in the RMP
    rc = -s_wij - s_h + σ - cpu * z_end - q_end * ((lamda * cpurq) / (get_prop(sn, p[end].src, :cpu_max) - get_prop(sn, p[end].src, :cpu_used)))
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

                problem.map_vars_paths[problem.vars[K_id][length(problem.vars[K_id])]] = p

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
                
                # update constraints (5)
                for i in keys(problem.vars_nodes_a)
                    if isInP(p, i) && (i != p[1].src)
                        push!(problem.c_energy_a[i], @constraint(problem.model, problem.vars_nodes_a[i] >= v))
                    end
                end
                # update constraints (6)

                lamda = 1/mean(com.arr_mean)
                cpurq = com.cpu_request

                set_normalized_coefficient(problem.c_energy_l[final_node], v, -(lamda * cpurq)/(get_prop(problem.sn, final_node, :cpu_max) - get_prop(problem.sn, final_node, :cpu_used)))

                # update constraints (17)
                for i in keys(problem.vars_nodes_a)
                    if isInP(p, i) && (i != p[1].src)
                        set_normalized_coefficient(problem.c_energy_normalized_a[i], v, -1/length(problem.commodities))
                    end
                end
            end
        end

        @suppress optimize!(problem.model)

        if onlyHasIntegerPathVars(problem)
            obj_integer = objective_value(problem.model)
            if obj_integer < tree.incumbent
                tree.incumbent = obj_integer

                pb_cpy = copyProblem(problem)
                optimize!(pb_cpy.model)
                tree.incumbent_solution = NodeBnB(obj_integer, obj_integer, typemax(Int64), nothing, pb_cpy, nothing, false, (nothing, nothing))
                println("Found new tree incumbent in CG: ", tree.incumbent)
                println("Gap is ", (tree.incumbent - tree.ub) / tree.ub)
                if (tree.incumbent - tree.ub) / tree.ub <= tree.threshold
                    tree.node_queue = PriorityQueue{Int,Tuple{Float64, Int}}()
                end
            end
        end
        #println(objective_value(problem.model))
        #println(onlyHasIntegerPathVars(problem))
        push!(objectives, objective_value(problem.model))
        #println(objective_value(problem.model))
        #if length(objectives) > 100 && length(unique(objectives[end-100:end])) == 1
        #    tree.column_generation_has_broken = true
        #    println("breakin")
        #    break
        #end
        #println("finding improving paths")
        improvingPaths = findImprovingPathsBrum(problem)
        #println("oopsie")
    end
    return problem.vars, problem.map_vars_paths, problem.paths
end
