using GraphIO
using EzXML
using Graphs
using MetaGraphs

using DataFrames, CSV

import Cairo, Fontconfig
using GraphPlot, Compose
function loadGraphAndMakeCapsFiles(filename, nodetypesFile, nodecapsFile, edgecapsFile)
    graph::SimpleGraph = loadgraph(filename, GraphIO.GraphML.GraphMLFormat())
    sn::MetaGraph = MetaGraph(graph)
    println(sn)
    nodelabel = 1:nv(sn)
    #draw(PNG("/tmp/Arn.png", 16cm, 16cm), gplot(sn, nodelabel=nodelabel))

    weightfield!(sn, :delay)

    set_prop!(sn, :switches, [])
    set_prop!(sn, :edges, [])
    set_prop!(sn, :anchors, [])
    set_prop!(sn, :central, [])

    i = 1
    for line in eachline(nodetypesFile)
        set_prop!(sn, i, :type, chomp(line))
        set_prop!(sn, i, :cpu_used, 0)
        if chomp(line) == "switch"
            push!(get_prop(sn, :switches), i)
            set_prop!(sn, i, :cpu_max, 0)
            set_prop!(sn, i, :energy_idle, 20)
            set_prop!(sn, i, :energy_max, 20)
        elseif chomp(line) == "anchor"
            push!(get_prop(sn, :anchors), i)
            set_prop!(sn, i, :cpu_max, 30)
            set_prop!(sn, i, :energy_idle, 0)
            set_prop!(sn, i, :energy_max, 0)
        elseif chomp(line) == "edge"
            push!(get_prop(sn, :edges), i)
            set_prop!(sn, i, :cpu_max, 120)
            set_prop!(sn, i, :energy_idle, 200)
            set_prop!(sn, i, :energy_max, 300)
        else
            push!(get_prop(sn, :central), i)
            set_prop!(sn, i, :cpu_max, 400)
            set_prop!(sn, i, :energy_idle, 350)
            set_prop!(sn, i, :energy_max, 550)
        end

        i += 1
    end
    df_nodes = DataFrame(nodes=[], cpu=[], energy_idle=[], energy_max=[])
    for n in 1:nv(sn)
        df_tmp = DataFrame(nodes=[n], cpu=[get_prop(sn, n, :cpu_max)], energy_idle=[get_prop(sn, n, :energy_idle)], energy_max=[get_prop(sn, n, :energy_max)])
        append!(df_nodes, df_tmp)
    end
    CSV.write(nodecapsFile, df_nodes)

    for e in edges(sn)
        capacity = rand(100:200)
        delay = rand(100:200)
        set_prop!(sn, e.src, e.dst, :BW_max, capacity)
        set_prop!(sn, e.src, e.dst, :BW_used, 0)
        set_prop!(sn, e.src, e.dst, :delay, delay)
    end
    df_edges = DataFrame(src=[], dst=[], delay=[], cap=[])
    for e in edges(sn)
        df_tmp = DataFrame(src=[e.src], dst=[e.dst], delay=[get_prop(sn, e, :delay)], cap=[get_prop(sn, e, :BW_max)])
        append!(df_edges, df_tmp)
    end
    CSV.write(edgecapsFile, df_edges)

    sn
end

function loadGraph(filename::String, nodetypesFile::String, nodecapFile::String, edgecapFile::String)
    graph::SimpleGraph = loadgraph(filename, GraphIO.GraphML.GraphMLFormat())
    sn::MetaGraph = MetaGraph(graph)
    #nodelabel = 1:nv(sn)
    #draw(PNG("/tmp/Arn.png", 16cm, 16cm), gplot(sn, nodelabel=nodelabel))

    weightfield!(sn, :delay)

    set_prop!(sn, :switches, [])
    set_prop!(sn, :edges, [])
    set_prop!(sn, :anchors, [])
    set_prop!(sn, :central, [])

    i = 1
    for line in eachline(nodetypesFile)
        set_prop!(sn, i, :type, chomp(line))
        set_prop!(sn, i, :cpu_used, 0)
        if chomp(line) == "switch"
            push!(get_prop(sn, :switches), i)
        elseif chomp(line) == "anchor"
            push!(get_prop(sn, :anchors), i)
        elseif chomp(line) == "edge"
            push!(get_prop(sn, :edges), i)
        else
            push!(get_prop(sn, :central), i)
        end
        i += 1
    end

    df_nodes = DataFrame(CSV.File(nodecapFile))
    df_edges = DataFrame(CSV.File(edgecapFile))
    
    for r in eachrow(df_nodes)
        set_prop!(sn, r[:nodes], :cpu_max, r[:cpu])
        set_prop!(sn, r[:nodes], :energy_max, r[:energy_max])
        set_prop!(sn, r[:nodes], :energy_idle, r[:energy_idle])
    end

    for r in eachrow(df_edges)
        set_prop!(sn, r[:src], r[:dst], :BW_max, r[:cap])
        set_prop!(sn, r[:src], r[:dst], :BW_used, 0)
        set_prop!(sn, r[:src], r[:dst], :delay, r[:delay])
    end
    sn
end

function undirectedGraphToDirected(G::MetaGraph)
    G2::MetaDiGraph = MetaDiGraph(nv(G) + 2 * length(edges(G)))
    weightfield!(G2, :delay)

    set_props!(G2, props(G))
    for n in 1:nv(G)
        set_props!(G2, n, deepcopy(props(G, n)))
    end

    # transformation of page 689 of Ahuja book (exercise 17.21)
    free_nodes = nv(G)+1
    for e in edges(G)
        i_prime = free_nodes
        j_prime = free_nodes+1
        i = e.src
        j = e.dst

        props_prime = Dict(:cpu_max=>0, :energy_idle=>0, :energy_max=>0, :cpu_used=>0, :type=>"dummy")
        set_props!(G2, i_prime, props_prime)
        set_props!(G2, j_prime, deepcopy(props_prime))

        add_edge!(G2, i, i_prime)
        add_edge!(G2, j, i_prime)
        add_edge!(G2, j_prime, i)
        add_edge!(G2, j_prime, j)
        add_edge!(G2, i_prime, j_prime)

        props_edge_prime = Dict(:delay=>0, :BW_max=>typemax(Int64), :BW_used=>0)

        set_props!(G2, i, i_prime, props_edge_prime)
        set_props!(G2, j, i_prime, deepcopy(props_edge_prime))
        set_props!(G2, j_prime, i, deepcopy(props_edge_prime))
        set_props!(G2, j_prime, j, deepcopy(props_edge_prime))

        set_props!(G2, i_prime, j_prime, props(G, i, j))
        set_prop!(G2, i_prime, j_prime, :original_arc, (i, j))
        free_nodes += 2
    end
    G2
end