function BranchAndPrice!(tree)
    while !terminated(tree)
        if time() - tree.startTime > tree.timeout
            tree.has_timedout = true
            return
        end
        node = getNextNode(tree)
        tree.num_nodes_treated += 1
        if node.parent != nothing && node.parent.problem != nothing
            # we only build the problem at execution time to reduce memory usage (and accelerate execution)
            buildProblem!(node)

            if node.sibling.hasBeenSolved
                node.parent.problem = nothing
            end
        end
        # if the current node is worse than the incumbent before being evaluated (eg parent is worse than incumbent), eliminate
        #println("There are ", length(tree.node_queue), " nodes in queue")
        if node.lb >= tree.incumbent
            #println("closing 1")
            #close_node!(tree, node)
            continue
        end
        
        evaluateNode!(node, tree)
        if node.sibling != nothing && node.sibling.hasBeenSolved
            ub = min(node.ub, node.sibling.ub)
            update_gap!(node.parent, ub, tree)
        end
        #println("evaluated!")
        
        # if node is worse than incumbent, move on to next node
        if node.lb >= tree.incumbent
            #println("closing 2")
            #close_node!(tree, node)
            continue
        end
        # we use big-M relaxation so infeasibilities are indicated by nonzero artificial vars
        if hasNonZeroArtificialVariables(node.problem)
            #println("closing 3")
            #println(node.problem.forbidden_arcs)
            #println(node.problem.forbidden_arcs_constraints)
            #close_node!(tree, node)
            continue
        end

        # if solution is integer, we update incumbent then close node
        if onlyHasIntegerPathVars(node.problem)
            obj_integer = objective_value(node.problem.model)
            if obj_integer < tree.incumbent
                tree.incumbent = obj_integer
                tree.incumbent_solution = node
                println("Found new tree incumbent: ", tree.incumbent)
                println("Gap is ", (tree.incumbent - tree.ub) / tree.ub)
                if (tree.incumbent - tree.ub) / tree.ub <= tree.threshold
                    tree.node_queue = PriorityQueue{Int,Tuple{Float64, Int}}()
                end
            end
            continue
        end
        # otherwise we have to branch on the node
        #@show node.id
        #@show tree.num_nodes
        #println("branching")
        branch!(tree, node)
    end
end

function buildProblem!(node::NodeBnB)
    problem = @suppress copyProblem(node.parent.problem)
    forbidden_com, OA = node.newlyForbiddenArcs
    
    problem.forbidden_arcs[forbidden_com] = union(problem.forbidden_arcs[forbidden_com], Set(OA))

    addForbiddenArcConstraints!(problem, OA, forbidden_com)

    node.problem = problem
end

function terminated(tree)
    return length(tree.node_queue) == 0
end

function getNextNode(tree)
    node_id = dequeue!(tree.node_queue)
    return tree.nodes[node_id]
end

function evaluateNode!(node, tree)
    @suppress optimize!(node.problem.model)

    problem = node.problem
    for k in problem.commodities 
        sn_with_dual_costs = makeBrumGraphEnergy(problem, k)
        for (idx, p) in enumerate(problem.paths[k.id])
            #r = getReducedCost(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
            ##println(value(problem.vars[k.id][idx]))
            ##println(r, " ", reduced_cost(problem.vars[k.id][idx]))
            #r2 = getReducedCostOK(sn_with_dual_costs, Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}(p), k, problem.c_edge, problem.c_commodity, problem.c_node, problem.forbidden_arcs[k.id], problem.knapsackCoverCutsNodes)
            #if !(isapprox(r, r2, atol=1e-8))
            #    println("very good:", reduced_cost(problem.vars[k.id][idx]))
            #    println("good: ", r2)
            #    println("bad: ", r)
            #    stop()
            #end
#
            #if !isapprox(r, reduced_cost(problem.vars[k.id][idx]), atol=1e-8) && r!= Inf
            #    println("very good:", reduced_cost(problem.vars[k.id][idx]))
            #    println("good: ", r2)
            #    println("bad: ", r)
            #    stop()
            #end
        end
    end
    # problem can be infeasible if removing some arcs/setting a(v_i) disconnects the graph.
    if termination_status(node.problem.model) == MOI.INFEASIBLE
        setfield!(node, :lb, Inf)
        return
    end

    setfield!(node, :lb, objective_value(node.problem.model))
    node.ub = node.lb

    added_cuts = true
    
    while added_cuts
        improvingPaths = findImprovingPathsBrum(node.problem)    
        vars, map_vars_paths, paths = columnGeneration(node.problem, improvingPaths, tree)
        
        
        
        if !hasNonZeroArtificialVariables(node.problem)
            #@suppress added_cuts = generateKnapsackCoverCuts(node.problem)
            added_cuts = false
        end

        setfield!(node, :lb, objective_value(node.problem.model))
        setfield!(node.problem, :vars, vars)
        setfield!(node.problem, :map_vars_paths, map_vars_paths)
        setfield!(node.problem, :paths, paths)


        if hasNonZeroArtificialVariables(node.problem)
            break
        end
    end

    node.hasBeenSolved = true
    #if node.id == 1 && !hasNonZeroArtificialVariables(node.problem)
    #    addGomoryCuts(node.problem)
    #end
    #stop()
end

function hasNonZeroArtificialVariables(problem)
    for av in problem.artificial_variables
        if value(av) > 0 && !isapprox(value(av), 0.0, atol=1e-8)
            return true
        end
    end
    return false
end

function onlyHasIntegerPathVars(problem::Problem)
    for k in problem.vars
        for v in k[2]
            if !(isapprox(value(v), 1.0, atol=1e-8) || isapprox(value(v), 0.0, atol=1e-8))
                return false
            end
        end
    end
    return true
end

function branch!(tree, node)
    OA1, OA2, forbidden_com = getOASets(node.problem)

    #problem2 = copyProblem(node.problem)
    #problem2.forbidden_arcs[forbidden_com] = union(node.problem.forbidden_arcs[forbidden_com], Set(OA1))

    #addForbiddenArcConstraints!(problem2, OA1, forbidden_com)

    # now that the model is clean of any forbidden arc, we can put it in the node
    n1 = NodeBnB(node.lb, node.ub, tree.num_nodes + 1, node, nothing, nothing, false, (forbidden_com, OA1))

    addNode!(tree, n1)

    #problem3 = copyProblem(node.problem)
    #problem3.forbidden_arcs[forbidden_com] = union(node.problem.forbidden_arcs[forbidden_com], Set(OA2))

    #addForbiddenArcConstraints!(problem3, OA2, forbidden_com)

    # now that the model is clean of any forbidden arc, we can put it in the node
    n2 = NodeBnB(node.lb, node.ub, tree.num_nodes + 1, node, nothing, n1, false, (forbidden_com, OA2))
    n1.sibling = n2   

    addNode!(tree, n2)
end

function addForbiddenArcConstraints!(problem::Problem, new_cons, commodity_id)
    copy_paths = deepcopy(problem.paths)
    for p in copy_paths[commodity_id]
        for e in p
            if e in new_cons
                for (k, v) in problem.map_vars_paths
                    if v == p
                        set_normalized_coefficient(problem.forbidden_arcs_constraints, k, 1)
                    end
                end
            end
        end
    end
end


function selectSplittedCommodity(problem)
    for K in keys(problem.vars)
        s = 0
        num_coms = 0
        comms = []
        for p in problem.vars[K]
            if !isapprox(value(p), 0.0, atol=1e-8) && !isapprox(value(p), 1.0, atol=1e-8) 
                num_coms += 1
                push!(comms, value(p))
            end
        end
        if num_coms > 1
            return K
        end
    end
    return nothing
end

function getOASets(problem)
    com_id = selectSplittedCommodity(problem)
    @assert com_id != nothing
    p1, p2 = getBranchingPaths(problem, com_id)
    fork_node, n1, n2 = getForkNode(problem.map_vars_paths[p1], problem.map_vars_paths[p2])
    OA1, OA2 = getForbiddenArcs(problem, fork_node, n1, n2, com_id)
    return OA1, OA2, com_id
end

function getForbiddenArcs(problem, fork_node, n1, n2, com_id)
    neigh = deepcopy(outneighbors(problem.sn, fork_node))
    filter!(e->e ≠ n1, neigh)
    filter!(e->e ≠ n2, neigh)
    for e in problem.forbidden_arcs[com_id]
        if e.src == fork_node
            filter!(a->a ≠ e.dst, neigh)
        end
    end

    f1, f2 = [Edge(fork_node, n1)], [Edge(fork_node, n2)]
    for (index, n) in enumerate(neigh)
        if index % 2 == 0
            push!(f1, Edge(fork_node, n))
        else
            push!(f2, Edge(fork_node, n))
        end
    end
    return f1, f2
end

function getBranchingPaths(problem::Problem, com_id::Int64)
    max_1 = 0
    p1 = []
    v1 = nothing
    max_2 = 0
    p2 = []
    v2 = nothing
    for p in problem.vars[com_id]    
        if !isapprox(value(p), 0, atol=1e-8) && !isapprox(value(p), 1.0, atol=1e-8)
            if value(p) > max_2
                if value(p) > max_1
                    max_2 = max_1
                    p2 = p1
                    v2 = v1
                    max_1 = value(p)
                    p1 = problem.map_vars_paths[p]
                    v1 = p
                else
                    max_2 = value(p)
                    p2 = problem.map_vars_paths[p]
                    v2 = p
                end
            end
        end
    end
    return v1, v2
end

function getForkNode(p1, p2)
    p1 = edgePathToNodePath(p1, p1[1].src)
    p2 = edgePathToNodePath(p2, p2[1].src)
    i = 0
    while p1[i+1] == p2[i+1]
        i+=1
    end
    return p1[i], p1[i+1], p2[i+1]
end

function update_gap!(node, ub, tree)
    node.ub = ub
    if node.sibling == nothing
        if tree.ub <= ub    
            tree.ub = ub
            println("Gap is ", (tree.incumbent - tree.ub) / tree.ub)
            if (tree.incumbent - tree.ub) / tree.ub <= tree.threshold
                tree.node_queue = PriorityQueue{Int,Tuple{Float64, Int}}()
            end
        end
    elseif node.sibling.hasBeenSolved
        ub = min(node.sibling.ub, node.ub)
        update_gap!(node.parent, ub, tree)
    end
end