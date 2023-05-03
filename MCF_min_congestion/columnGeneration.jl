include("labelCorrectingForbidden.jl")
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
        
        #println("pre-pricing")
        paths, delays = solvePricing(sn_with_dual_costs, k.source, k.meta_destination, forbidden_paths, k.dests)
        #println(problem.paths[k.id])
        #println(paths)
        #println(k.dests)
        for i in 1:length(paths)
            p = nodePathToEdgePath(paths[i])
            if delays[i] * 2 < k.max_delay
                reduced_cost = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.c_cpu_cong, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
                if reduced_cost < 0 && !isapprox(reduced_cost, 0.0, atol=1e-18)  #&& !(paths[i] in paths_seen)# && !has_already_been_added(p, original_paths[k.id])
                    if paths_to_add[k.id] == nothing || reduced_cost < paths_to_add[k.id]["reduced cost"]
                        #println(p)
                        paths_to_add[k.id] = Dict("reduced cost" => reduced_cost, "path" => p)
                        #break
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
            set_prop!(sn_modified, e, :w2Ratio, 0)
            #println("here")
            #stop()
        else
            BW = commodity.BW_demanded
            # set w_ij
            set_prop!(sn_modified, e, :wTimesBW, shadow_price(problem.c_edge[e]) * BW)
            #println(shadow_price(problem.c_edge[e]) * BW)
            #set_prop!(sn_modified, e, :sumDualCuts, 0)
            # time is simply propagation delay
            #set_prop!(sn_modified, e, :propagationDelay, get_prop(sn_modified, e, :weight))
            if get_prop(sn_modified, e, :BW_max) <= 100000
                set_prop!(sn_modified, e, :w2Ratio, -(-shadow_price(problem.c_BW_cong[e])  * BW / (get_prop(sn_modified, e, :BW_max) - get_prop(sn_modified, e, :BW_used))))
                #println((-shadow_price(problem.c_BW_cong[e])  * BW / (get_prop(sn_modified, e, :BW_max) - get_prop(sn_modified, e, :BW_used))))
            else
                set_prop!(sn_modified, e, :w2Ratio, 0)
            end
        end
    end

    # we add abs(smallest edge) +1 to all edges in order to ensure there are no negative arcs (and hence no negative cycle)
    # the +1 ensures we do not have edges that have (0, 0) costs, which would cause issues when retrieving paths (see retrieve_paths function in labelCorrectingSP.jl)
    min_value = Inf
    for e in edges(sn_modified)
        #if get_prop(sn_modified, e, :wTimesBW) + get_prop(sn_modified, e, :sumDualCuts) < min_value 
        if get_prop(sn_modified, e, :wTimesBW) < min_value 
            min_value = get_prop(sn_modified, e, :wTimesBW) + get_prop(sn_modified, e, :w2Ratio)
        end
    end
    set_prop!(sn_modified, :normalization_wTimesBW, abs(min_value)+1)
    for e in edges(sn_modified)
        set_prop!(sn_modified, e, :wTimesBW, get_prop(sn_modified, e, :wTimesBW) + abs(min_value)+1)
    end

    for n in 1:nv(problem.sn)
        for n2 in inneighbors(problem.sn, n)
            set_prop!(sn_modified, n2, n, :sumWandH, get_prop(sn_modified, n2, n, :wTimesBW) + get_prop(sn_modified, n2, n, :w2Ratio))
        end
    end
    return sn_modified
end

"""
function getReducedCostOK(sn::MetaDiGraph, p::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}, commodity::Commodity, c_edge, c_commodity, c_node, c_cpu_cong, forbidden_arcs, knapsackCutsNodes, c_BW_cong)
    s_wij = 0
    s_w2 = 0
    for e in p
        # even when setting costs to large values, we might still have a path with forbidden arc returned here, if its the only possible path in graph
        # in this case we return infinity
        if e in forbidden_arcs
            return Inf
        end
        BW = commodity.BW_demanded
        s_wij += shadow_price(c_edge[e]) * BW
        #println(s_wij)
        if get_prop(sn, e, :BW_max) < 100000
            s_w2  +=  shadow_price(c_BW_cong[e])  * BW / (get_prop(sn, e, :BW_max) - get_prop(sn, e, :BW_used))
            #println(s_w2)
        end
    end

#    for e in edges(sn)
#        println(e, " ", shadow_price(c_BW_cong[e]), " ", shadow_price(c_edge[e]))
#    end
#
#    println()
#    for n in keys(c_node)
#        println(n, " ", shadow_price(c_node[n]), " ", shadow_price(c_cpu_cong[n]))
#    end

#    s_h = 0
#    for e in p
#        s_h += get_prop(sn, e.dst, :h) - get_prop(sn, :normalization_h)
#    end
    σ = shadow_price(c_commodity[commodity.id])

    cpu = getCpu(sn, commodity, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p))

    # we look at the SOURCE of last edge (and not destination) because the last node is the meta-destination which we do not care about
    z_end = shadow_price(c_node[p[end].src])
    for cut in knapsackCutsNodes[p[end].src]
        z_end += shadow_price(cut)
    end


    w_1 = shadow_price(c_cpu_cong[p[end].src])
    println()
    println(w_1 * (cpu / get_prop(sn, p[end].src, :cpu_max)))
    println(cpu * z_end)



    # this reduced cost is only valid if the path is NOT in the RMP
    rc = -s_wij -s_w2 + σ - cpu * z_end - w_1 * (cpu / (get_prop(sn, p[end].src, :cpu_max) - get_prop(sn, p[end].src, :cpu_used)))
    @show (rc)
    return rc
end
"""

function getReducedCost(sn::MetaDiGraph, p::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}, commodity::Commodity, c_edge, c_commodity, c_node, c_cpu_cong, forbidden_arcs, knapsackCutsNodes)
    s_wij = 0
    s_w2 = 0
    for e in p
        # even when setting costs to large values, we might still have a path with forbidden arc returned here, if its the only possible path in graph
        # in this case we return infinity
        if e in forbidden_arcs
            return Inf
        end
        BW = commodity.BW_demanded
        s_wij += get_prop(sn, e, :wTimesBW) - get_prop(sn, :normalization_wTimesBW)
        if get_prop(sn, e, :BW_max) <= 100000
            s_w2  += get_prop(sn, e, :w2Ratio)
        end
    end

#    s_h = 0
#    for e in p
#        s_h += get_prop(sn, e.dst, :h) - get_prop(sn, :normalization_h)
#    end
    σ = shadow_price(c_commodity[commodity.id])

    cpu = getCpu(sn, commodity, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p))

    # we look at the SOURCE of last edge (and not destination) because the last node is the meta-destination which we do not care about
    z_end = shadow_price(c_node[p[end].src])
    for cut in knapsackCutsNodes[p[end].src]
        z_end += shadow_price(cut)
    end


    w_1 = shadow_price(c_cpu_cong[p[end].src])

    #@show s_wij
    #@show s_w2
    #@show σ
    #@show cpu * z_end
    #@show w_1 * (cpu / (get_prop(sn, p[end].src, :cpu_max) - get_prop(sn, p[end].src, :cpu_used)))
#
    #println()

    # this reduced cost is only valid if the path is NOT in the RMP
    rc = -s_wij -s_w2 + σ - cpu * z_end - w_1 * (cpu / (get_prop(sn, p[end].src, :cpu_max) - get_prop(sn, p[end].src, :cpu_used)))
    return rc
end

function columnGeneration(problem::Problem, improvingPaths::Dict{Any, Any})
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
                
                # update constraints (6)
                set_normalized_coefficient(problem.c_cpu_cong[final_node], v, -cpu/(get_prop(problem.sn, final_node, :cpu_max) - get_prop(problem.sn, final_node, :cpu_used)))

                for e in p
                    if get_prop(problem.sn, e, :BW_max) < 100000
                        BW_d = com.BW_demanded
                        BW_max = get_prop(problem.sn, e, :BW_max) - get_prop(problem.sn, e, :BW_used)
                        set_normalized_coefficient(problem.c_BW_cong[e], v, -BW_d/BW_max)
                    end
                end
            end
        end

        @suppress optimize!(problem.model)

        #println(objective_value(problem.model))
        #println(onlyHasIntegerPathVars(problem))
        push!(objectives, objective_value(problem.model))
        if length(objectives) > 10 && length(unique(objectives[end-10:end])) == 1
#            tree.column_generation_has_broken = true
            println("breakin")
            break
        end
        #println("finding improving paths")
        improvingPaths = findImprovingPathsBrum(problem)
        #println("oopsie")
    end
    return problem.vars, problem.map_vars_paths, problem.paths
end
