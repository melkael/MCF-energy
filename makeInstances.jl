using JLD
using Graphs
using Random
using GraphPlot, Compose
include("./MCF_energy/structs.jl")
include("commodityGenerator.jl")
include("./MCF_energy/graphLoader.jl")
include("./MCF_energy/utils.jl")

# example usage
# julia makeInstances.jl ./instances/Agis/agis.graphml 20 ./instances/Agis/agis.nodetypes ./instances/Agis/agis.nodecaps ./instances/Agis/agis.edgecaps 3 [0.87,0.9,0.95] 600 1200 20 40 60 90 0 2 8 5 15 ./instances/Agis/agis1test.jld


function main(ARGS)
    graph_file = ARGS[1]
    num_commodities = parse(Int64, ARGS[2])
    nodetypes_file = ARGS[3]
    nodecapsFile = ARGS[4]
    edgecapsFile = ARGS[5]
    num_spikes = parse(Int64,   ARGS[6])
    sla        = eval(Meta.parse(ARGS[7]))
    min_delay  = parse(Float64, ARGS[8])
    max_delay  = parse(Float64, ARGS[9])
 
    min_BW     = parse(Float64, ARGS[10])
    max_BW     = parse(Float64, ARGS[11])
    min_cpu    = parse(Float64, ARGS[12])
    max_cpu    = parse(Float64, ARGS[13])


    min_mean_arr = parse(Int64, ARGS[15])
    max_mean_arr = parse(Int64, ARGS[16])
    min_sigma    = parse(Int64, ARGS[17])
    max_sigma    = parse(Int64, ARGS[18])

    Random.seed!(parse(Int64, ARGS[14]))

    sn = loadGraphAndMakeCapsFiles(graph_file, nodetypes_file, nodecapsFile, edgecapsFile)
    sn_directed = undirectedGraphToDirected(sn)
    commodities = generateFlows(sn_directed, num_commodities, num_spikes, sla, union([i for i in min_delay:100:max_delay], max_delay), min_BW, max_BW, min_cpu, max_cpu, min_mean_arr, max_mean_arr, min_sigma, max_sigma)

    println(save(ARGS[19], "commodities", commodities))
end

#main(ARGS)