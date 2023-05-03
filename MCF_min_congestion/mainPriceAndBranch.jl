include("structs.jl")
include("graphLoader.jl")
include("utils.jl")
include("commodityGenerator.jl")
include("lpModel.jl")
include("cpu.jl")
include("columnGeneration.jl")
include("branchAndPrice.jl")
include("cutGeneration.jl")
using Random
using DataStructures

using HiGHS

OPT = HiGHS
BIG_M = 99999

global Gr
global Fp

function main()
    Random.seed!(2)
    sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/Cogentco.graphml", "/home/maxime/Téléchargements/Cogentco.nodetypes")
    sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/UsCarrier.graphml", "/home/maxime/Téléchargements/Kdl.nodetypes")
    sn = loadGraph("/home/maxime/Téléchargements/MCF_final/Arn.graphml", "/home/maxime/Téléchargements/MCF_final/Arn.nodetypes")
    #sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/Ans.graphml", "/home/maxime/Téléchargements/MCF_branch_mieux/Ans.nodetypes")
   
    # sn = loadGraph("/home/maxime/Téléchargements/topology_zoo/Kdl.graphml", "/home/maxime/Téléchargements/Kdl2.nodetypes")
    sn_directed = undirectedGraphToDirected(sn)

    commodities = generateFlows(sn_directed)

    sn = addMetaDestinations(sn_directed, commodities)

    paths = getInitialPaths(sn, commodities)

    model, c_edge, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, edge_contains_paths = makeModel(sn, paths, commodities)

    #println(paths)
    #stop()
    pb = Problem(sn, commodities, paths, model, c_edge, edge_contains_paths, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, Dict(k.id=>Set() for k in commodities), Dict(i => [] for i in 1:nv(sn)), Dict(i => [] for i in 1:nv(sn)))

    @suppress optimize!(pb.model)
    # problem can be infeasible if removing some arcs/setting a(v_i) disconnects the graph.
    if termination_status(pb.model) == MOI.INFEASIBLE
        return Inf
    end
    
    improvingPaths = findImprovingPathsBrum(pb)    

    vars, map_vars_paths, paths = columnGeneration(pb, improvingPaths)
    
    if hasNonZeroArtificialVariables(pb)
        return Inf
    end

    setfield!(pb, :vars, vars)
    setfield!(pb, :map_vars_paths, map_vars_paths)
    setfield!(pb, :paths, paths)

    @show objective_value(pb.model)

    for k in keys(vars)
        for v in vars[k]
            set_binary(v)
        end
    end
    optimize!(pb.model)
    return objective_value(pb.model)
end

main()
