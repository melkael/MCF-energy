include("structs.jl")
include("graphLoader.jl")
include("utils.jl")
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

function main(ARGS)
    path = "../instances/" * ARGS[1] * "/"
    path = joinpath(@__DIR__, path)

    
    sn = loadGraph(path * ARGS[1] * ".graphml", path * ARGS[1] * ".nodetypes", path * ARGS[1] * ".nodecaps", path * ARGS[1] * ".edgecaps")
    sn_directed = undirectedGraphToDirected(sn)

    commodities = load(path * ARGS[1] * ARGS[2] * ".jld")["commodities"]

    sn = addMetaDestinations(sn_directed, commodities)

    paths = getInitialPaths(sn, commodities)

    model, c_edge, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, edge_contains_paths = makeModel(sn, paths, commodities)

    #println(paths)
    #stop()
    pb = Problem(sn, commodities, paths, model, c_edge, edge_contains_paths, c_commodity, c_node, c_energy_a, c_energy_normalized_a, c_energy_l, vars, vars_nodes_a, vars_nodes_l, cpus, artificial_variables, map_vars_paths, c_energy_a_paths_map, constraint_forbidden, Dict(k.id=>Set() for k in commodities), Dict(i => [] for i in 1:nv(sn)), Dict(i => [] for i in 1:nv(sn)))

    root = NodeBnB(floatmax(Float64)-1, floatmax(Float64)-1, 1, nothing, pb, nothing, false, (nothing, nothing))
    q = PriorityQueue{Int,Tuple{Float64, Int}}()
    d = Dict{Int, NodeBnB}()
    q[1] = (floatmax(Float64)-1, 1)
    d[1] = root
    tree = BnBTree(sn, commodities, Inf, nothing, Inf, q, d, 1, [], -Inf, time(), 3600, false, 0.01, 0, false)

    res = @timed BranchAndPrice!(tree)
    println(tree.incumbent)
    if tree.incumbent_solution != nothing
        energy = get_energy_consumption(tree.incumbent_solution.problem)
        resources = get_resource_consumption(tree.incumbent_solution.problem)
        congestion = get_congestion(tree.incumbent_solution.problem)
        return DataFrame(instance=ARGS[2], energy=[energy], resources=[resources], congestion=[congestion], time=[res[2]], num_nodes=[tree.num_nodes_treated], timeout=[tree.has_timedout], cg_broke=[tree.column_generation_has_broken], gap=[(tree.incumbent - tree.ub) / tree.ub])
    end
    return nothing
end

dfMain = DataFrame(instance=[], energy=[], resources=[], congestion=[], time=[], num_nodes=[], timeout=[], cg_broke=[], gap=[])
df = main(split(ARGS[1] * " " * string(1)))
for i in 1:parse(Int64, ARGS[2])
    df = main(split(ARGS[1] * " " * string(i)))
    if df != nothing
        append!(dfMain, df)
        @show dfMain
        CSV.write(ARGS[3] * ARGS[1] * ".csv", dfMain)
    else
        break
    end
end
