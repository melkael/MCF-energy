using JuMP

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
            println("Commodity ", k.id, " has no suitable path at all. Problem is infeasible, try with looser delay. Delay: ", k.max_delay)
            error()
        end 
    end
    paths
end


function makeModel(sn::MetaDiGraph, paths::Dict, commodities::Vector{Commodity})
    artificial_variables = []
    
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
        c_commodity[K_id] = @constraint(model, s >= 1-av)
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

    obj = 0
    for K_id in keys(paths)
        for (index, p) in enumerate(paths[K_id])
            com = retrieve_commodity(commodities, K_id)
            obj += vars[K_id][index] * (com.BW_demanded * length_non_dummy(sn, p) + cpus[K_id][index])
        end
    end


    for av in artificial_variables
        obj += BIG_M * av
    end
    @objective(model, Min, obj)
    constraint_forbidden = @constraint(model, 0 == 0)
    println("done")
    return model, c_edge, c_commodity, c_node, vars, cpus, artificial_variables, map_vars_paths, constraint_forbidden, edge_contains_paths
end