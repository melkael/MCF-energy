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

function generateFlows(sn::MetaDiGraph, numFlows, numSpikes, slas, delays, minBW, maxBW, minCPU, maxCPU, min_mean_arr, max_mean_arr, min_sigma, max_sigma)
    flows = Vector{Commodity}()
    for id in 1:numFlows
        source = sample(get_prop(sn, :anchors))
        dests = vcat(get_prop(sn, :edges), get_prop(sn, :central))

        arr_coeffs = rand(Dirichlet(numSpikes, 1.0))
        arr_mean = [rand(min_mean_arr:max_mean_arr) for i in numSpikes]
        arr_mean = []
        arr_sigma = []
        for i in 1:numSpikes
            push!(arr_mean, rand(min_mean_arr:max_mean_arr))
            push!(arr_sigma, rand(min_sigma:max_sigma))
        end
        sla = sample(slas)
        max_delay = sample(delays)
        BW_demanded = rand(minBW:maxBW)
        cpu_demanded = rand(minCPU:maxCPU)

        f = Commodity(id+=1, source, dests, arr_mean, arr_sigma, arr_coeffs, sla, max_delay, cpu_demanded, BW_demanded, 0)
        push!(flows, f)
    end
    flows
end

function addMetaDestinations(G::MetaDiGraph, commodities::Vector{Commodity})
    G2::MetaDiGraph = MetaDiGraph(nv(G) + length(commodities))
    weightfield!(G2, :delay)
    
    for n in 1:nv(G)
        set_props!(G2, n, props(G, n))
    end
    for e in edges(G)
        add_edge!(G2, e)
        set_props!(G2, e, props(G, e))
    end

    for k in 1:length(commodities)
        metaD = k + nv(G)
        commodities[k].meta_destination = metaD
        for d in commodities[k].dests
            add_edge!(G2, d, metaD)
            set_props!(G2, d, metaD, Dict(:BW_max=>2*BIG_M, :BW_used=>0, :delay=>0))
        end
    end

    for n in nv(G)+1:nv(G2)
        set_props!(G2, n, Dict(:cpu_max=>0, :cpu_used=>0, :energy_idle=>0, :energy_max=>0, :type=>"dummy"))
    end
    G2
end