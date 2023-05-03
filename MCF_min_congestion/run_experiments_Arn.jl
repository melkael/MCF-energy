include("structs.jl")
include("graphLoader.jl")
include("utils.jl")
include("commodityGenerator.jl")
include("lpModel.jl")
include("cpu.jl")
include("columnGeneration.jl")
include("branchAndPrice.jl")
include("cutGeneration.jl")
include("sequential.jl")

using JLD
using DataFrames
using Suppressor
using CSV

using HiGHS

OPT = HiGHS
BIG_M = 99999

# first run sequential solution
function main()
    #sn = loadGraph("/home/maxime/Téléchargements/MCF_branch_mieux/instances/Arn/Arn.graphml", "/home/maxime/Téléchargements/MCF_branch_mieux/instances/Arn/Arn.nodetypes")
    folder = "/home/maxime/Téléchargements/MCF_branch_mieux/instances/Arn/"
    sn = loadGraph(folder * "Arn.graphml", folder * "Arn.nodetypes", folder * "Arn.nodecap", folder * "Arn.edgecap")
    sn_directed = undirectedGraphToDirected(sn)
    df_result = DataFrame(obj=[], time=[], seq=[])
    for exp in ["3"]
        commodities = load("/home/maxime/Téléchargements/MCF_branch_mieux/instances/Arn/Arn" * exp * ".jld")["commodities"]
        commodities = commodities[1:end-3]
        for i in 1:length(commodities)
            s = collect(i:i:length(commodities)) ∪ length(commodities)
            println(s)
            #println(commodities[1])
            (objective, runtime) = runSeqOpt(s, deepcopy(commodities), deepcopy(sn_directed), 3600)
            println(objective, " ", runtime)
            tmp = DataFrame(obj=[objective], time=[runtime], seq=[s])
            append!(df_result, tmp)
            CSV.write("/home/maxime/Téléchargements/MCF_branch_mieux/results/res_testons2_sequential_Arn" * exp * ".csv", df_result)
        end
    end
end

main()
