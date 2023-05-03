using JuMP, Statistics

function getInitialPaths(sn, commodities)
    paths = Dict()
    num_good = 0
    for k in commodities
        paths[k.id] = []
        sp = a_star(sn, k.source, k.meta_destination)
        dist = 0
        for e in sp
            dist += get_prop(sn, e, :delay)
        end
        if dist*2 < k.max_delay
            push!(paths[k.id], sp)
        end
        if paths[k.id] == []
            error("Commodity" * String(k.id) * "has no suitable path at all. Problem is infeasible, try with looser delay")
        end 
    end
    paths
end


function makeModel(sn::MetaDiGraph, paths::Dict, commodities::Vector{Commodity})
    artificial_variables = []

    vars_nodes_a = Dict()
    vars_nodes_l = Dict()

    vars = Dict()
    cpus = Dict()
    map_vars_paths = Dict()
    model = Model(OPT.Optimizer)
    #set_optimizer_attribute(model, "presolve", "off")
    # declaring LP variables
    println(1)
    for K_id in keys(paths)
        vars[K_id] = []
        cpus[K_id] = []
        for p in paths[K_id]
            push!(vars[K_id], @variable(model, lower_bound=0))
            map_vars_paths[vars[K_id][length(vars[K_id])]] = p
            
            cpu = getCpu(sn, retrieve_commodity(commodities, K_id), p)
            #@show cpu
            push!(cpus[K_id], cpu)
        end
    end

    for i in 1:nv(sn)
        if get_prop(sn, i, :type) != "dummy" && get_prop(sn, i, :type) != "anchor" 
            vars_nodes_a[i] = @variable(model, lower_bound=0, upper_bound=1)
            vars_nodes_l[i] = @variable(model, lower_bound=0, upper_bound=1)
        end
    end
    
    c_edge = Dict()
    edge_contains_paths = Dict()
    println("Defining constraints (2)")
    for e in edges(sn)
        edge_contains_paths[e] = []
        s = 0
        for K_id in keys(paths)
            for (index, p) in enumerate(paths[K_id])
                for i in 1:length(p)
                    if (p[i].src == e.src && p[i].dst == e.dst)
                        s += vars[K_id][index] * retrieve_commodity(commodities, K_id).BW_demanded
                        push!(edge_contains_paths[e], (vars[K_id][index], retrieve_commodity(commodities, K_id).BW_demanded))
                    end
                end
            end
        end
        av = @variable(model, lower_bound=0)
        push!(artificial_variables, av)
        c_edge[e] = @constraint(model, s <= (get_prop(sn, e, :BW_max) - get_prop(sn, e, :BW_used))+av)
    end

    c_commodity = Dict()
    println("Defining constraints (3)")
    for K_id in keys(paths)
        s = 0
        for (index, p) in enumerate(paths[K_id])
            s += vars[K_id][index]
        end
        av = @variable(model, lower_bound=0)
        push!(artificial_variables, av)
        c_commodity[K_id] = @constraint(model, s == 1-av)
    end

    # defining constraints (4)
    c_node = Dict()
    println("Defining constraints (4)")
    for i in 1:nv(sn)
        s = 0
        for K_id in keys(paths)
            for (index, p) in enumerate(paths[K_id])
                # we look at the n-1th node of path because the last one is the meta destination
                if i == p[end-1].dst
                    s += cpus[K_id][index] * vars[K_id][index]
                end
            end
        end
        av = @variable(model, lower_bound=0)
        push!(artificial_variables, av)
        c_node[i] = @constraint(model, s <= get_prop(sn, i, :cpu_max) - get_prop(sn, i, :cpu_used)+av)
    end



    c_energy_a = Dict()
    c_energy_a_paths_map = Dict()
    println("Defining constraints (5)")
    for i in keys(vars_nodes_a)
        c_energy_a[i] = []
        c_energy_a_paths_map[i] = []
        for K_id in keys(paths)
            for (index, p) in enumerate(paths[K_id])
                if is_in_arc_path(p, i) && (i != p[1].src)
                    push!(c_energy_a[i], @constraint(model, vars_nodes_a[i] >= vars[K_id][index]))
                    push!(c_energy_a_paths_map[i], p)
                end
            end
        end
    end

    c_energy_l = Dict()
    println("Defining constraints (6)")
    for i in keys(vars_nodes_l)
        s = 0
        for K_id in keys(paths)
            for (index, p) in enumerate(paths[K_id])
                # don't forget we have a dummy meta destination so use n-1th node as final for cpu stuff
                final_node = p[end-1].dst
                if final_node == i
                    lamda = 1/mean(retrieve_commodity(commodities, K_id).arr_mean)
                    cpurq = retrieve_commodity(commodities, K_id).cpu_request
                    println(lamda*cpurq/ (get_prop(sn, i, :cpu_max) - get_prop(sn, i, :cpu_used)))
                    s += vars[K_id][index] * lamda * cpurq / (get_prop(sn, i, :cpu_max) - get_prop(sn, i, :cpu_used))
                end
            end
        end
        c_energy_l[i] = @constraint(model, vars_nodes_l[i] >= s)
    end
    obj = 0

    c_energy_normalized_a = Dict()
    println("Defining constraints (17)")
    for i in keys(vars_nodes_a)
        s = 0
        for K_id in keys(paths)
            for (index, p) in enumerate(paths[K_id])
                if is_in_arc_path(p, i) && (i != p[1].src)
                    s += vars[K_id][index]
                end
            end
        end
        c_energy_normalized_a[i] = @constraint(model, vars_nodes_a[i] >= s * (1/length(commodities)))
    end

    for i in 1:nv(sn)
        if i in keys(vars_nodes_a)
            obj += vars_nodes_a[i] * get_prop(sn, i, :energy_idle) + vars_nodes_l[i] * (get_prop(sn, i, :energy_max) - get_prop(sn, i, :energy_idle))
        end
    end

    for av in artificial_variables
        obj += BIG_M * av
    end
    @objective(model, Min, obj)
    constraint_forbidden = @constraint(model, 0 == 0)
    println("done")
    return model, c_edge, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, edge_contains_paths
end