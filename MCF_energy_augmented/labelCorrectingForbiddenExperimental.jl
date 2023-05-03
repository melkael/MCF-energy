using Graphs, MetaGraphs, DataStructures, PrettyPrint
import Base.parse
using Profile

macro timeout(seconds, expr, fail)
    quote
        tsk = @task $expr
        schedule(tsk)
        Timer($seconds) do timer
            istaskdone(tsk) || Base.throwto(tsk, InterruptException())
        end
        try
            fetch(tsk)
        catch _
            $fail
        end
    end
end

function ⊕(A, x)
    new_A = deepcopy(A)
    for i in 1:length(A)
        for j in 1:length(x)
            new_A[i][j] += x[j]
        end
    end
    return new_A
end

function filterForbidden(A, forbidden_paths, num_labels)
    new_A = deepcopy(A)
    for i in 1:length(A)
        for k in 1:length(forbidden_paths)
            if A[i][k+num_labels] >= length(forbidden_paths[k]) - 1
                filter!(x->x!=A[i], new_A)
                continue
            end
        end
    end
    return new_A
end

function filterConstraints(A, valueConstraints, num_labels)
    new_A = deepcopy(A)
    for i in 1:length(A)
        for k in 1:length(valueConstraints)
            if A[i][k+num_labels] >= valueConstraints[k]
                filter!(x->x!=A[i], new_A)
                continue
            end
        end
    end
    return new_A
end


function all_less(x, y)
    reduce(&, (x .<= y))
end

function all_equal(x, y)
    reduce(&, (x .== y))
end


function naive_merge(X, Y)
    X_cpy = deepcopy(X)
    Y_cpy= deepcopy(Y)
    for x in X
        for y in Y
            if all_less(x, y) && !all_equal(x, y)
                filter!(e->e ≠ y, Y_cpy)
            end
        end
    end
    for y in Y
        for x in X
            if all_less(y, x) && !all_equal(y, x)
                filter!(e->e ≠ x, X_cpy)
            end
        end
    end
    M = vcat(X_cpy, Y_cpy)
    return unique(M)
end

function naive_merge2(X, Y)
    #X = deepcopy(X)
    #Y = deepcopy(Y)
    toFilterY = []
    toFilterX = []
    for x in X
        for y in Y
            if !all_equal(x, y)
                if all_less(x, y)
                    push!(toFilterY, y)
                    #filter!(e->e != y, Y)
                elseif all_less(y, x)
                    push!(toFilterX, x)
                    #filter!(e->e != x, X)
                end
            end
        end
    end
    #println(Y)
    filter!(e->!(e ∈ toFilterY), Y)
    #println(Y)
    filter!(e->!(e ∈ toFilterX), X)
    #M = vcat(X, Y)
    ret = []
    for x in X
        push!(ret, x)
    end
    for y in Y
        if !(y ∈ ret)
            push!(ret, y)
        end
    end
    return ret
end

function dominated(elt, X, numObjs)
    for x in X
        if all_less(x[1:numObjs], elt[1:numObjs]) || all_equal(elt[1:numObjs], x[1:numObjs])
            return true
        end
    end
    return false
end

function addForbiddenLabels!(G, forbidden_paths)
    for e in edges(G)
        label = zeros(Float32, length(forbidden_paths))
        for path_idx in 1:length(forbidden_paths)
            p = forbidden_paths[path_idx]
            for i in 1:length(p)-1
                if e.src == p[i] && e.dst == p[i+1]
                    label[path_idx] = 1
                end
                if e.dst == p[i] && e.src == p[i+1]
                    label[path_idx] = 1
                end
            end
        end
        set_prop!(G, e, :forbiddenWeights, label)
    end
end

function get_weights(G, n1, n2, labelNames)
    map(x->get_prop(G, n1, n2, x), labelNames)
end

function labelCorrectingSPExperimental(G, source, labelsObjective, labelsConstraints, maxValueConstraints, forbidden_paths, targeet=nothing)
    # value constranits merges forbidden paths constraints and normla constraints
    valueConstraints = []
    for p in forbidden_paths
        push!(valueConstraints, length(p)-1)
    end
    valueConstraints = vcat(valueConstraints, maxValueConstraints)

    addForbiddenLabels!(G, forbidden_paths)
    D = Dict()
    DM = []
    Labeled = DataStructures.Queue{Int}()
    for i in 1:nv(G)
        D[i] = []
    end
    D[source] = [zeros(Float32, length(labelsObjective) + length(forbidden_paths) + length(labelsConstraints))]
    numConsideredDomination = length(labelsObjective) + length(forbidden_paths)
    enqueue!(Labeled, source)
    while length(Labeled) > 0
        u = dequeue!(Labeled)
        for j in outneighbors(G, u)
            #println(j, " ", u)
            weightsObj = get_weights(G, u, j, labelsObjective)
            weightsConstraints = get_weights(G, u, j, labelsConstraints)
            forbiddenLabels = get_prop(G, u, j, :forbiddenWeights)
            A = D[j]
            B = D[u] ⊕ vcat(weightsObj, forbiddenLabels, weightsConstraints)
            DM = naive_merge2(A, B)
            DM = filterConstraints(DM, valueConstraints, length(labelsObjective))
            if targeet != j && targeet != nothing && D[targeet] != []
                DM_cpy = deepcopy(DM)
                for elt in DM
                    if dominated(elt, D[targeet], length(labelsObjective))
                        filter!(x->x!=elt, DM_cpy)
                    end
                end
                DM = DM_cpy
            end
            if DM != D[j]
                D[j] = deepcopy(DM)
                if !(j in Labeled)
                    enqueue!(Labeled, j)
                end
            end
        end
    end
    return D
end

function removeDominated!(labels, num_labels)
    to_remove = []
    for x in labels
        for l in labels
            if all_less(x[1:num_labels], l[1:num_labels]) && !all_equal(x[1:num_labels], l[1:num_labels])
                push!(to_remove, x)
                break
            end
        end
    end
    filter!(x->!(x in to_remove), labels)
end

function get_matching_label(labels, G, n, node, l, labelsObjective, labelsConstraints)
    for l2 in labels[n]
        weightsObjective = get_weights(G, n, node, labelsObjective)
        weightsConstraints = get_weights(G, n, node, labelsConstraints)
        forbiddenLabels = get_prop(G, n, node, :forbiddenWeights)
        s = [l2] ⊕ vcat(weightsObjective, forbiddenLabels, weightsConstraints)
        if all_equal(s[1], l)
            return l2
        end
    end
    return nothing
end

function retrievePaths(labels, G, d, labelsObjective, labelsConstraints, forbidden_paths)
    paths = []
    for l in labels[d]
        path = [d]
        node = d
        while sum(labels[node][1]) != 0
            for n in inneighbors(G, node)
                new_label = get_matching_label(labels, G, n, node, l, labelsObjective, labelsConstraints)
                if new_label != nothing
                    node = n
                    l = new_label
                    push!(path, node)
                    break
                end
            end
        end
        reverse!(path)
        push!(paths, path)
    end
    return paths
end

function retrieveFirstPath(labels, G, d, labelsObjective, labelsConstraints, forbidden_paths)
    path = nothing
    for l in labels[d]
        path = [d]
        node = d
        while sum(labels[node][1]) != 0
            for n in inneighbors(G, node)
                new_label = get_matching_label(labels, G, n, node, l, labelsObjective, labelsConstraints)
                if new_label != nothing
                    node = n
                    l = new_label
                    push!(path, node)
                    break
                end
            end
        end
        reverse!(path)
        break
    end
    return path
end


function load_graph(node_number, time_file, distance_file)
    G = MetaGraph(node_number)
    for l in readlines(time_file)
        l = split(l, " ")
        src = parse(Int64, l[2])
        dst = parse(Int64, l[3])
        time = parse(Int64, l[4])
        add_edge!(G, src, dst)
        set_prop!(G, src, dst, :time, time)
    end
    for l in readlines(distance_file)
        l = split(l, " ")
        src = parse(Int64, l[2])
        dst = parse(Int64, l[3])
        dist = parse(Int64, l[4])
        add_edge!(G, src, dst)
        set_prop!(G, src, dst, :dist, dist)
    end
    G
end

function NetMaker(num_nodes, min_arcs, max_arcs, node_interval, num_constraints, numObj)
    labelsConstraints = []
    labelsObjective = []
    G = cycle_graph(num_nodes)
    G = MetaGraph(G)
    for n in 1:nv(G)
        num_arcs = rand(min_arcs:max_arcs)
        already_linked = deepcopy(neighbors(G, n))
        push!(already_linked, n)
        selected_neigh = already_linked[1]
        for i in 1:num_arcs
            while selected_neigh in already_linked
                selected_neigh = mod(rand(n - node_interval:n + node_interval),nv(G))
            end
            add_edge!(G, n, selected_neigh)
        end
    end
    for i in 1:num_constraints
        push!(labelsConstraints, Symbol("c" * string(i)))
    end
    for i in 1:numObj
        push!(labelsObjective, Symbol("o" * string(i)))
    end
    for e in edges(G)
        for o in labelsObjective
            set_prop!(G, e, o, rand(1:33))
        end
        #set_prop!(G, e, :time, rand(1:33))
        #set_prop!(G, e, :dist, rand(67:100))
        for l in labelsConstraints
            set_prop!(G, e, l, rand(5:10))
        end
    end
    G, labelsConstraints, labelsObjective
end

function toDigraph(G)
    dG = MetaDiGraph(nv(G))
    for e in edges(G)
        add_edge!(dG, e.src, e.dst, props(G, e))
        add_edge!(dG, e.dst, e.src, props(G, e))
    end
    dG
end

function addDummyNodes(G, labels::Vector{Symbol})
    G2 = MetaDiGraph(nv(G)*2)
    for e in edges(G)
        properties::Dict{Symbol, Union{Float64, Int64, Tuple{Int64, Int64}}} = props(G::MetaDiGraph, e::Graphs.SimpleGraphs.SimpleEdge{Int64})
        add_edge!(G2, e.src, e.dst, properties)
    end
    for n in 1:nv(G)
        properties = Dict()
        for l in labels
            properties[l] = 0.0
        end
        add_edge!(G2, n, n+nv(G), properties)
    end
    G2
end

function solvePricingExperimental(snWithDual, source, metaDest, forbiddenPaths, dests)
    #G = deepcopy(snWithDual)
    labelsObjective = [:sumWandH, :delay]
    G = addDummyNodes(snWithDual, labelsObjective)

#    forbiddenPaths = deepcopy(forbiddenPaths)

#    println(forbiddenPaths)

    for p in forbiddenPaths
        #push!(forbiddenPaths, vcat(p[1:end-1], [p[end-1]+nv(snWithDual)]))
        p[end] = p[end-1] + nv(snWithDual)
    end
    
    labels = labelCorrectingSPExperimental(G, source, labelsObjective, [], [], forbiddenPaths) #metaDest+nv(snWithDual))

    #paths = retrievePaths(labels, G, metaDest+nv(snWithDual), labelsObjective, [], forbiddenPaths)
    #println(paths)
    paths = []
    for d in dests
        paths = vcat(paths, retrievePaths(labels, G, d+nv(snWithDual), labelsObjective, [], forbiddenPaths))
    end

    m = map(x->vcat(x[1:end-1], [metaDest]), paths)

    delays = []
    for p in paths
        d = 0
        for i in 1:length(p) -1
            d += get_prop(G, p[i], p[i+1], :delay)
        end
        push!(delays, d)
    end
    
    return m, delays
end

function solvePricingSingleObjective(snWithDual, source, metaDest, forbiddenPaths, dests)
    labelsObjective = [:sumWandH, :delay]
    G = addDummyNodes(snWithDual, labelsObjective)
    labelsObjective = [:sumWandH]

    for p in forbiddenPaths
        #push!(forbiddenPaths, vcat(p[1:end-1], [p[end-1]+nv(snWithDual)]))
        p[end] = p[end-1] + nv(snWithDual)
    end
    
    labels = labelCorrectingSPExperimental(G, source, labelsObjective, [], [], forbiddenPaths) #metaDest+nv(snWithDual))

    #paths = retrievePaths(labels, G, metaDest+nv(snWithDual), labelsObjective, [], forbiddenPaths)
    #println(paths)
    paths = []
    for d in dests
        paths = vcat(paths, retrievePaths(labels, G, d+nv(snWithDual), labelsObjective, [], forbiddenPaths))
    end

    m = map(x->vcat(x[1:end-1], [metaDest]), paths)

    delays = []
    for p in paths
        d = 0
        for i in 1:length(p) -1
            d += get_prop(G, p[i], p[i+1], :delay)
        end
        push!(delays, d)
    end
    
    return m, delays
end
