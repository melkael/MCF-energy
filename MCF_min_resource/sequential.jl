function runSeqOpt(seq, commodities, sn_directed, lim)
    lim_orig = lim
    incumbent = Inf
    tree = nothing
    for i in seq
        i = floor(Int64, i)
        sn = addMetaDestinations(sn_directed, commodities[1:i])
        paths = getInitialPaths(sn, commodities[1:i])

        model, c_edge, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, edge_contains_paths = @suppress makeModel(sn, paths, commodities[1:i])

        #println(paths)
        #stop()
        pb = Problem(sn, commodities[1:i], paths, model, c_edge, edge_contains_paths, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, Dict(k.id=>Set() for k in commodities[1:i]), Dict(i => [] for i in edges(sn)), Dict(i => [] for i in 1:nv(sn)))

        root = NodeBnB(floatmax(Float64)-1, floatmax(Float64)-1, 1, nothing, pb, nothing, false, nothing)
        q = PriorityQueue{Int,Tuple{Float64, Int}}()
        d = Dict{Int, NodeBnB}()
        q[1] = (floatmax(Float64)-1, 1)
        d[1] = root
        tree = BnBTree(sn, commodities[1:i], Inf, nothing, Inf, q, d, 1, [], -Inf, time(), lim, false)
        
        t = @timed BranchAndPrice!(tree)
        lim -= t[2]
        incumbent = tree.incumbent
        println(incumbent)
        if incumbent == Inf
            return Inf, lim_orig -lim
        end
        if tree.incumbent_solution != nothing
            for k in keys(tree.incumbent_solution.problem.vars)
                for v in tree.incumbent_solution.problem.vars[k]
                    if isapprox(value(v), 1, atol=1e-8)
                        retrieve_commodity(commodities, k).dests = [tree.incumbent_solution.problem.map_vars_paths[v][end].src]
                        #println(tree.incumbent_solution.problem.map_vars_paths[v][end].src)
                    end
                    @assert isapprox(value(v), 1, atol=1e-8) || isapprox(value(v), 0, atol=1e-8)
                end
            end
        end
    end

    
    if tree != nothing
        for k in keys(tree.incumbent_solution.problem.vars)
            for (index, p) in enumerate(tree.incumbent_solution.problem.paths[k])
                if isapprox(value(tree.incumbent_solution.problem.vars[k][index]), 1, atol=1e-8)
                    println("Value: ", value(tree.incumbent_solution.problem.vars[k][index]))
                    println("Commodity: ", k)
                    println("Host: ", p[end].src)
                    println(length(tree.incumbent_solution.problem.paths[k]))
                    println(tree.incumbent_solution.problem.cpus[k])
                    println("CPU consumed: ", tree.incumbent_solution.problem.cpus[k][index])
                    println("Path: ", p)
                    println()
                end
            end
        end
    end
    
    return incumbent, lim_orig - lim
end
