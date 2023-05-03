using HiGHS

function Base.isless(a::VariableRef, b::VariableRef)
    return value(a) < value(b)
end

function partitionCover(cover)
    C1, C2 = [], []
    for v in cover
        if value(v) == 1
            push!(C2, v)
        else
            push!(C1, v)
        end
    end
    C1, C2
end

function cutGeneration(problem)
    for c in problem.c_edge
        if length(problem.edge_contains_paths[c[1]]) > 1
            tot_BW_max = 0
            for (p, BW) in problem.edge_contains_paths[c[1]]
                tot_BW_max += BW 
            end
            # entering this condition means we found an edge constraint which has a chance to be saturated (and hence is interesting to consider for lifting)
            if tot_BW_max > normalized_rhs(c[2])
                nonZ = getNonZeroVariables(problem.edge_contains_paths[c[1]])
                #if 1 in map(x->value(x), nonZ) && !all(y->value(y)==value(nonZ[1]), nonZ)
                if !all(y->value(y)==value(nonZ[1]), nonZ)
                    #println(nonZ)
                    #for v in nonZ
                    #    println(value(v), " ", v)
                    #end
                    #println()
                    #sort!(nonZ, rev=true)
                    #for v in nonZ
                    #    println(value(v))
                    #end

                    sort!(nonZ, rev=true)
                    cover = []
                    bound = normalized_rhs(c[2])
                    for v in nonZ
                        bound -= normalized_coefficient(c[2], v)
                        push!(cover, v)
                        # it has to exceed b from the Gu paper
                        if bound < 0
                            break
                        end
                    end
                    C1, C2 = partitionCover(cover)
                    while !partitionnedCoverIsMinimal(C1, C2, c[2])
                        removeFromC1(C1, C2, c[2])
                    end
                    println(C1)
                    println(C2)
                    stop()
                end
            end
        end
    end
end

function generateKnapsackCoverCutsSCPA(problem)
    init_obj = objective_value(problem.model)
    added_cuts = false
    """
    consToAddEdges = []
    for k in keys(problem.c_edge)
        w = []
        xbar = []
        for com in keys(problem.vars)
            for x in problem.vars[com]
                if normalized_coefficient(problem.c_edge[k], x) > 0
                    push!(w, normalized_coefficient(problem.c_edge[k], x))
                    push!(xbar, x)
                end
            end
        end
        b = normalized_rhs(problem.c_edge[k])
        n = length(w)
        if length(w) > 1 && b <= 10000
            mod = Model(HiGHS.Optimizer)
            F = Dict()

            for i in 1:b
                F[i] = @variable(mod, lower_bound=0)
            end
            for i in 1:b
                for j in 1:b
                    if i+j <= b
                        @constraint(mod, F[i] + F[j] - F[i+j] <= 0)
                    end
                end
            end
            @constraint(mod, F[b] == 1)

            @objective(mod, Max, sum([value(xbar[i]) * F[w[i]] for i in 1:n]) - F[b])
            @suppress optimize!(mod)
            if objective_value(mod) > 0
                println("adding scpa")
                push!(consToAddEdges[k], (F, w, xbar, n, b))
            end
        end
    end
    """
    consToAddNodes = Dict(i => [] for i in 1:nv(problem.sn))
    for k in keys(problem.c_node)
        w = []
        xbar = []
        for com in keys(problem.vars)
            for x in problem.vars[com]
                if normalized_coefficient(problem.c_node[k], x) > 0
                    push!(w, normalized_coefficient(problem.c_node[k], x))
                    push!(xbar, x)
                end
            end
        end
        b = normalized_rhs(problem.c_node[k])
        n = length(w)
        if length(w) > 1 && reduce(&, w .<= b)
            mod = Model(HiGHS.Optimizer)
            F = Dict()

            for i in 1:b
                F[i] = @variable(mod, lower_bound=0)
            end
            for i in 1:b
                for j in 1:b
                    if i+j <= b
                        @constraint(mod, F[i] + F[j] - F[i+j] <= 0)
                    end
                end
            end
            @constraint(mod, F[b] == 1)

            @objective(mod, Max, sum([value(xbar[i]) * F[w[i]] for i in 1:n]) - F[b])
            @suppress optimize!(mod)

            if objective_value(mod) > 0
                push!(consToAddNodes[k], (F, w, xbar, n, b))
                added_cuts = true
            end
        end
    end

    #for (F, w, xbar, n, b) in consToAddEdges
    #    push!(problem.knapsackCoverCutsEdges, @constraint(problem.model, sum([xbar[i] * value(F[w[i]]) for i in 1:n]) <= value(F[b])))
    #end

    for node in keys(consToAddNodes)
        for (F, w, xbar, n, b) in consToAddNodes[node]
            push!(problem.knapsackCoverCutsNodes[node], @constraint(problem.model, sum([xbar[i] * value(F[w[i]]) for i in 1:n]) <= value(F[b])))
        end
    end

    @suppress optimize!(problem.model)
    return added_cuts & !isapprox(objective_value(problem.model), init_obj, atol=1e-8)
end

function generateKnapsackCoverCuts(problem)
    init_obj = objective_value(problem.model)
    added_cuts = false
    """
    consToAddEdges = Dict(e => [] for e in edges(problem.sn))
    for k in keys(problem.c_edge)
        w = []
        xbar = []
        for com in keys(problem.vars)
            for x in problem.vars[com]
                if normalized_coefficient(problem.c_edge[k], x) > 0
                    push!(w, normalized_coefficient(problem.c_edge[k], x))
                    push!(xbar, x)
                end
            end
        end
        b = normalized_rhs(problem.c_edge[k])
        n = length(w)
        if length(w) > 1
            c, zeta =  generateCover(n, w, b, xbar)
            if zeta < 1 && !isapprox(zeta, 1, atol=1e-8)
                alpha = balasLifting(c, n, w, b)
                push!(consToAddEdges[edge], (xbar, c, alpha))
            end
        end
    end
    """
    consToAddNodes = Dict(i => [] for i in 1:nv(problem.sn))
    for k in keys(problem.c_node)
        w = []
        xbar = []
        for com in keys(problem.vars)
            for x in problem.vars[com]
                if normalized_coefficient(problem.c_node[k], x) > 0
                    push!(w, normalized_coefficient(problem.c_node[k], x))
                    push!(xbar, x)
                end
            end
        end
        b = normalized_rhs(problem.c_node[k])
        n = length(w)
        if length(w) > 1
            c, zeta =  generateCover(n, w, b, xbar)
            if zeta < 1 && !isapprox(zeta, 1, atol=1e-8)
                alpha = balasLifting(c, n, w, b)
                push!(consToAddNodes[k], (xbar, c, alpha))
            end
        end
    end

    for k in keys(problem.c_node)
        for (xbar, c, alpha) in consToAddNodes[k]
            push!(problem.knapsackCoverCutsNodes[k], @constraint(problem.model, sum([xbar[i] for i in c]) + sum(xbar[i] * alpha[i] for i in keys(alpha)) <= length(c) - 1))
        end
    end
    """
    for e in keys(consToAddEdges)
        for (xbar, c, alpha) in consToAddEdges[e]
            push!(problem.knapsackCoverCutsEdges[e], @constraint(problem.model, sum([xbar[i] for i in c]) + sum(xbar[i] * alpha[i] for i in keys(alpha)) <= length(c) - 1))
        end
    end
    """
    @suppress optimize!(problem.model)
    return added_cuts & !isapprox(objective_value(problem.model), init_obj, atol=1e-8)
end

function generateCover(n, w, b, xbar)
    m = Model(HiGHS.Optimizer)
    z = [@variable(m, lower_bound=0, binary=true) for i in 1:n]
    @constraint(m, sum([w[j] * z[j] for j in 1:n]) >= b+1)
    @objective(m, Min, sum([(1-value(xbar[j])) * z[j] for j in 1:n]))
    @suppress optimize!(m)
    if primal_status(m) == MOI.NO_SOLUTION
        return nothing, 2
    end
    ret = []
    for j in 1:n
        if value(z[j]) > 0
            push!(ret, j)
        end
    end
    ret, objective_value(m)
end

function sequentialLifting(c, n, w, b)
    alpha = Dict(k => 0.0 for k in setdiff(1:n, c))
    m = Model(HiGHS.Optimizer)
    x_dash = Dict(i => @variable(m, lower_bound=0.0, binary=true) for i in c)
    toBechecked = collect(keys(alpha))
    checked = []

    if length(toBechecked) != 0
        i = popfirst!(toBechecked)
        x_dash[i] = @variable(m, lower_bound = 0.0, binary=true)
        tmpc = sum([x_dash[j] * w[j] for j in c]) + sum([alpha[j] * x_dash[j] for j in checked])
        tmpo = sum([x_dash[j] for j in c]) + sum([alpha[j] * x_dash[j] for j in checked])

        @objective(m, Max, tmpo)
        @constraint(m, tmpc <= b - w[i])
        @suppress optimize!(m)
        if length(c) - 1 - objective_value(m) >= 1
            alpha[i] = length(c) - 1 - objective_value(m)
            push!(checked, i)
        end
    end
    alpha
end

function balasLifting(c, n, w, b)
    orderCoverWeightsList = copy(c)
    lookup = Dict()
    for i in c
        lookup[i] = w[i]
    end
    sort!(orderCoverWeightsList, by=x->lookup[x], rev=true)
    alpha = Dict(k => 0.0 for k in setdiff(1:n, c))
    mu = Dict(0=>0)
    ind = 0
    cum = 0
    for k in orderCoverWeightsList
        if ind == 0
            mu[ind+1] = w[k]
            ind += 1
            cum += w[k]
        else
            mu[ind+1] = w[k] + cum
            ind += 1
            cum += w[k]
        end
    end
    lamb = mu[length(mu)-1] - b
    for j in keys(alpha)
        h_list = []
        for h in keys(mu)
            if h < length(mu) - 1 && mu[h] <= w[j] <= mu[h+1] - lamb
                push!(h_list, h)
            end
        end
        if length(h_list) != 0
            alpha[j] = h_list[1]
        else
            h_list = []
            for h in keys(mu)
                if h < length(mu) - 1 && mu[h] - lamb <= w[j] <= mu[h]
                    push!(h_list, h)
                    alpha[j] = h[end]
                end
            end
        end
    end
    alpha
end
