using MetaGraphs
using Distributions

function generateFlows(sn::MetaDiGraph)
    flows = Vector{Commodity}()
    for id in 1:12
        source = sample(get_prop(sn, :anchors))
        dests = vcat(get_prop(sn, :edges), get_prop(sn, :central))

        arr_coeffs = rand(Dirichlet(2, 1.0))
        arr_mean = [rand(2:4), rand(1:2)]
        arr_sigma = [rand(3:20), rand(5:12)]
        sla = sample([0.85, 0.9, 0.95])
        max_delay = rand(Uniform(1200,2000))
        BW_demanded = rand(18:25)
        cpu_demanded = rand(60:90)

        f = Commodity(id+=1, source, dests, arr_mean, arr_sigma, arr_coeffs, sla, max_delay, cpu_demanded, BW_demanded, 0)
        push!(flows, f)
    end
    return flows
end

function generateFlows(sn::MetaDiGraph, numFlows, numSpikes, slas, delays, minBW, maxBW, minCPU, maxCPU)
    flows = Vector{Commodity}()
    for id in 1:numFlows
        source = sample(get_prop(sn, :anchors))
        dests = vcat(get_prop(sn, :edges), get_prop(sn, :central))

        arr_coeffs = rand(Dirichlet(numSpikes, 1.0))
        arr_mean = [rand(min_mean:max_mean) for i in numSpikes]
        arr_sigma = [rand(min_sigma:max_sigma) for i in numSpikes]
        sla = sample(slas)
        max_delay = sample(delays)
        BW_demanded = rand(minBW:maxBW)
        cpu_demanded = rand(minCPU:maxCPU)

        f = Commodity(id+=1, source, dests, arr_mean, arr_sigma, arr_coeffs, sla, max_delay, cpu_demanded, BW_demanded, 0)
        push!(flows, f)
    end
    return flows
end