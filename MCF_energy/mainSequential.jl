include("structs.jl")
include("graphLoader.jl")
include("utils.jl")
include("commodityGenerator.jl")
include("lpModel.jl")
include("cpu.jl")
include("columnGeneration.jl")
include("branchAndPrice.jl")
include("cutGeneration.jl")
#using Random
using DataStructures
using NOMAD
using Distributions
using StatsBase

using Suppressor

using HiGHS

OPT = HiGHS
BIG_M = 99999

global Gr
global Fp

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
        pb = Problem(sn, commodities[1:i], paths, model, c_edge, edge_contains_paths, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, Dict(k.id=>Set() for k in commodities[1:i]), Dict(i => [] for i in 1:nv(sn)), Dict(i => [] for i in 1:nv(sn)))

        root = NodeBnB(floatmax(Float64)-1, floatmax(Float64)-1, 1, nothing, pb, nothing, false)
        q = PriorityQueue{Int,Tuple{Float64, Int}}()
        d = Dict{Int, NodeBnB}()
        q[1] = (floatmax(Float64)-1, 1)
        d[1] = root
        tree = BnBTree(sn, commodities[1:i], Inf, nothing, Inf, q, d, 1, [], -Inf, time(), lim, false)
        
        t = @timed BranchAndPrice!(tree)
        lim -= t[2]
        incumbent = tree.incumbent
        if incumbent == Inf
            return Inf
        end
        if tree.incumbent_solution != nothing
            for k in keys(tree.incumbent_solution.problem.vars)
                for v in tree.incumbent_solution.problem.vars[k]
                    if isapprox(value(v), 1, atol=1e-8)
                        retrieve_commodity(commodities, k).dests = [tree.incumbent_solution.problem.map_vars_paths[v][end].src]
                    end
                    @assert isapprox(value(v), 1, atol=1e-8) || isapprox(value(v), 0, atol=1e-8)
                end
            end
        end
    end

    """
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
    """
    
    return incumbent, lim_orig - lim
end

function main()
    #Random.seed!(42)
    #sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/Cogentco.graphml", "/home/maxime/Téléchargements/Cogentco.nodetypes")
    #sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/UsCarrier.graphml", "/home/maxime/Téléchargements/Kdl.nodetypes")
    #sn = loadGraph("/home/maxime/Téléchargements/MCF_final/Arn.graphml", "/home/maxime/Téléchargements/MCF_final/Arn.nodetypes")
    #sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/Ans.graphml", "/home/maxime/Téléchargements/MCF_branch_mieux/Ans.nodetypes")
   
    # sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/Kdl.graphml", "/home/maxime/Téléchargements/Kdl2.nodetypes")
    sn_directed = undirectedGraphToDirected(sn)
    commodities = load("/home/maxime/Téléchargements/MCF_branch_mieux/instances/Arn/Arn1.jld")["commodities"]#generateFlows(sn_directed)

    #save("/home/maxime/Téléchargements/MCF_branch_mieux/instances/AGIS/AGIS1.jld", "commodities", commodities)
    #stop()
    
#    runSeqOpt(, deepcopy(commodities), deepcopy(sn_directed))

    for i in 1:length(commodities)
        s = collect(i:i:length(commodities)) ∪ length(commodities)
        println(s)
        #println(commodities[1])
        @show (runSeqOpt(s, deepcopy(commodities), deepcopy(sn_directed), 600), i)
    end

    """
    current_seq = [i for i in 1:length(commodities)]
    current_seq = [5, 10, 15, 22, 28]
    current_value = runSeqOpt(current_seq, deepcopy(commodities), deepcopy(sn_directed))

    #seen = [seq]
    best_seq = current_seq
    best_value = current_value

    initial_temp = 90
    final_temp = .1
    alpha = 0.01
    
    current_temp = initial_temp

    while current_temp > final_temp
        @show best_value
        @show current_value
        println()
        
        neighbor = StatsBase.sample(collect(getNeighbors(current_seq, length(commodities))))
        res = runSeqOpt(neighbor, deepcopy(commodities), deepcopy(sn_directed))
        if res < best_value
            best_seq = neighbor
            best_value = res
        end

        if res < current_value || rand(Uniform(0, 1)) < exp((res - current_value) / current_temp)
            current_value = res
            current_seq = neighbor
        end       
    end
    """
end
main()
