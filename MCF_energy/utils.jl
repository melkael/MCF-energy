using Statistics

function retrieve_commodity(commodities::Vector{Commodity}, id::Int64)
    for k in commodities
        if k.id == id
            return k
        end
    end
    return nothing
end

function arc_path_delay(sn::MetaDiGraph, path::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}})
    delay = 0
    for i in 1:length(path)
        delay += get_prop(sn, path[i].src, path[i].dst, :delay)
    end
    return delay
end

function is_in_arc_path(p::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}, i::Int64)
    for e in p
        if i == e.src || i == e.dst
            return true
        end
    end
    return false
end

function edgePathToNodePath(p_edges, source)
    p_nodes = []
    if p_edges[1].src == source
        push!(p_nodes, p_edges[1].src)
        push!(p_nodes, p_edges[1].dst)
    else
        push!(p_nodes, p_edges[1].dst)
        push!(p_nodes, p_edges[1].src)
    end

    for e in p_edges[2:length(p_edges)]
        if e.src == p_nodes[length(p_nodes)]
            push!(p_nodes, e.dst)
        elseif e.dst == p_nodes[length(p_nodes)]
            push!(p_nodes, e.src)
        end
    end
    p_nodes
end

function nodePathToEdgePath(p_nodes)
    p = []
    for i = 1:length(p_nodes)-1
        push!(p, Edge(p_nodes[i], p_nodes[i+1]))
    end
    return p
end

function isForbidden(e, forbidden_arcs)
    return e in forbidden_arcs
end

function noPathsToAdd(paths_to_add::Dict{Any, Any})
    for k in keys(paths_to_add)
        if paths_to_add[k] != nothing
            return false
        end
    end
    return true
end

function isInP(p, i)
    for e in p
        if i == e.src || i == e.dst
            return true
        end
    end
    return false
end

function toSimplePath(p::Vector{Int64})
    #return p
    i = 1
    output = []
    while i <= length(p)
        j = i+1
        while j<= length(p)
            if p[i] == p[j]
                i = j
            end
            j = j + 1
        end
        push!(output, p[i])
        i+=1
    end 
    output
end

function copyProblem(problem::Problem)
    model_new, rf = copy_model(problem.model)
    c_edge_new = Dict()
    edge_contains_paths_new = deepcopy(problem.edge_contains_paths)
    for k in keys(problem.c_edge)
        c_edge_new[k] = rf[problem.c_edge[k]]
    end

    c_commodity_new = Dict()
    for k in keys(problem.c_commodity)
        c_commodity_new[k] = rf[problem.c_commodity[k]]
    end
    c_node_new = Dict()
    for k in keys(problem.c_node)
        c_node_new[k] = rf[problem.c_node[k]]
    end
    c_energy_a_new = Dict()
    for k in keys(problem.c_energy_a)
        if problem.c_energy_a[k] == nothing
            c_energy_a_new[k] = nothing
        else
            c_energy_a_new[k] = rf[problem.c_energy_a[k]]
        end
    end

    c_energy_l_new = Dict()
    for k in keys(problem.c_energy_l)
        if problem.c_energy_l[k] == nothing
            c_energy_l_new[k] = nothing
        else
            c_energy_l_new[k] = rf[problem.c_energy_l[k]]
        end
    end

    c_energy_normalized_a_new = Dict()

    for k in keys(problem.c_energy_normalized_a)
        c_energy_normalized_a_new[k] = rf[problem.c_energy_normalized_a[k]]
    end


    vars_new = Dict()
    for k in keys(problem.vars)
        vars_new[k] = map(x->rf[x], problem.vars[k])
        for i in 1:length(vars_new[k])
            set_start_value(vars_new[k][i], value(problem.vars[k][i]))
        end
    end

    vars_nodes_a_new = Dict()
    for k in keys(problem.vars_nodes_a)
        vars_nodes_a_new[k] = rf[problem.vars_nodes_a[k]]
    end

    vars_nodes_l_new = Dict()
    for k in keys(problem.vars_nodes_a)
        vars_nodes_l_new[k] = rf[problem.vars_nodes_l[k]]
    end


    #vars_nodes_a_new = map(x->rf[x], problem.vars_nodes_a)
    #vars_nodes_l_new = map(x->rf[x], problem.vars_nodes_l)
    for v in keys(vars_nodes_a_new)
        set_start_value(vars_nodes_a_new[v], value(problem.vars_nodes_a[v]))
        set_start_value(vars_nodes_l_new[v], value(problem.vars_nodes_l[v]))
    end

    
    artificial_variables_new = map(x->rf[x], problem.artificial_variables)
    for i in 1:length(artificial_variables_new)
        set_start_value(artificial_variables_new[i], value(problem.artificial_variables[i]))
    end
    
    cpus_new = deepcopy(problem.cpus)

    map_vars_paths_new = Dict()
    for k in keys(problem.map_vars_paths)
        map_vars_paths_new[rf[k]] = deepcopy(problem.map_vars_paths[k])
    end

    set_optimizer(model_new, OPT.Optimizer)
    set_silent(model_new)

    forbidden_arcs_constraints_new = rf[problem.forbidden_arcs_constraints]
    #for c in problem.forbidden_arcs_constraints
    #    push!(forbidden_arcs_constraints_new, rf[c])
    #end

    c_energy_a_paths_map_new = Dict()
    for i in 1:nv(problem.sn)
        c_energy_a_paths_map_new[i] = []
        for K_id in keys(problem.paths)
            for (index, p) in enumerate(problem.paths[K_id])
                if is_in_arc_path(Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), i) && (i != p[1].src)
                    push!(c_energy_a_paths_map_new[i], p)
                end
            end
        end
    end

    Problem(deepcopy(problem.sn), deepcopy(problem.commodities), deepcopy(problem.paths), model_new, c_edge_new, edge_contains_paths_new ,c_commodity_new, c_node_new, c_energy_a_new, c_energy_normalized_a_new, c_energy_l_new, vars_new, vars_nodes_a_new, vars_nodes_l_new, cpus_new, artificial_variables_new, map_vars_paths_new, c_energy_a_paths_map_new, forbidden_arcs_constraints_new, deepcopy(problem.forbidden_arcs), deepcopy(problem.knapsackCoverCutsEdges), deepcopy(problem.knapsackCoverCutsNodes))
    #return model_new, c_edge_new, c_commodity_new, c_node_new, c_energy_a_new, c_energy_normalized_a_new, c_energy_l_new, vars_new, vars_nodes_a_new, vars_nodes_l_new, cpus_new, min_cpus_new, artificial_variables_new, map_vars_paths_new, deepcopy(paths), deepcopy(forbidden_arcs), rf
end

function alreadySeen(forbidden_arcs, tree)
    for seen in tree.seen_nodes
        if seen == forbidden_arcs
            return true
        end
    end
    return false
end

function addNode!(tree::BnBTree, n::NodeBnB)
    if n.lb < tree.incumbent
        tree.nodes[n.id] = n
        tree.node_queue[n.id] = (n.lb, n.id)
        tree.num_nodes += 1
    end
end

function get_resource_consumption(problem)
    consumption = 0.0
    for (K_id, vars) in problem.vars
        for (index, v) in enumerate(vars)
            if !isapprox(value(v), 0, atol=1e-8)
                p = problem.map_vars_paths[v]
                for e in p
                    if get_prop(problem.sn, e, :BW_max) <= 100000
                        consumption += retrieve_commodity(problem.commodities, K_id).BW_demanded
                    end
                end
                consumption += problem.cpus[K_id][index]
            end
        end
    end
    consumption
end


function get_energy_consumption(problem)
    on_nodes = Set()
    usage = Dict(i => 0.0 for i in 1:nv(problem.sn))
    for vars in problem.vars
        for (index, v) in enumerate(vars[2])
            if !isapprox(value(v), 0, atol=1e-8)
                p = problem.map_vars_paths[v]
                for e in p
                    push!(on_nodes, e.dst)
                end
                com = retrieve_commodity(problem.commodities, vars[1])
                usage[p[end].src] += com.cpu_request * (1/mean(com.arr_mean))
            end
        end
    end
    consumption = 0.0
    for n in on_nodes
        consumption += get_prop(problem.sn, n, :energy_idle)
    end
    for n in keys(usage)
        if usage[n] != 0
            println(usage[n])
            consumption += (get_prop(problem.sn, n, :energy_max) - get_prop(problem.sn, n, :energy_idle)) * usage[n]/get_prop(problem.sn, n, :cpu_max)
        end
    end
    return consumption
end

function get_congestion(problem)
    congestion = 0.0
    usage = Dict(i => 0 for i in 1:nv(problem.sn))
    for vars in problem.vars
        for (index, v) in enumerate(vars[2])
            if !isapprox(value(v), 0, atol=1e-8)
                p = problem.map_vars_paths[v]
                usage[p[end].src] += problem.cpus[vars[1]][index]
            end
        end
    end
    for n in keys(usage)
        if usage[n] != 0
            node_cong = usage[n]/get_prop(problem.sn, n, :cpu_max)
            if node_cong > congestion
                congestion = node_cong
            end
        end
    end

    for e in edges(problem.sn)
        edge_cong = 0
        if get_prop(problem.sn, e, :BW_max) < 100000
            for (K_id, vars) in problem.vars
                for (index, v) in enumerate(vars)
                    if !(isapprox(value(v), 0, atol=1e-8)) && (e in problem.map_vars_paths[v])
                        edge_cong += retrieve_commodity(problem.commodities, K_id).BW_demanded / (get_prop(problem.sn, e, :BW_max) - get_prop(problem.sn, e, :BW_used))
                    end
                end
            end
        end
        if edge_cong > congestion
            congestion = edge_cong
        end
    end
    congestion
end

function addMetaDestinations(G::MetaDiGraph, commodities::Vector{Commodity})
    G2::MetaDiGraph = MetaDiGraph(nv(G) + length(commodities))
    weightfield!(G2, :delay)
    
    for n in 1:nv(G)
        set_props!(G2, n, props(G, n))
    end
    for e in edges(G)
        add_edge!(G2, e)
        set_props!(G2, e, props(G, e))
    end

    for k in 1:length(commodities)
        metaD = k + nv(G)
        commodities[k].meta_destination = metaD
        for d in commodities[k].dests
            add_edge!(G2, d, metaD)
            set_props!(G2, d, metaD, Dict(:BW_max=>2*BIG_M, :BW_used=>0, :delay=>0))
        end
    end

    for n in nv(G)+1:nv(G2)
        set_props!(G2, n, Dict(:cpu_max=>0, :cpu_used=>0, :energy_idle=>0, :energy_max=>0, :type=>"dummy"))
    end
    G2
end