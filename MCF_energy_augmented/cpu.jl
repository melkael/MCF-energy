function getRo(mu::Float64, arr_mean::Vector{Float64}, arr_coeffs::Vector{Float64})
    arr = 0
    for i in 1:length(arr_coeffs)
        arr += arr_coeffs[i] * arr_mean[i]
    end
    lambd = 1/arr
    return lambd/mu
end

function getNu(mu, arr_mean, arr_sigma, arr_coeffs)
    coeff_square = 0
    coeff_lin = 0
    coeff_const = 0
    # we use the simplified for of the polynom
    for i = 1:length(arr_mean)
        a = arr_coeffs[i]
        b = mu
        c = arr_mean[i]
        d = arr_sigma[i]
        coeff_square += 0.5 * a^2 * b^2 * d^2
        coeff_lin += -(a^2 * b^2 * d^2) + a*b*c
        coeff_const += 0.5 * a^2 * b^2 * d^2 - a*b*c
    end

    # we use the log-form for numerical stability
    function f(x)
        return coeff_square * x^2 + coeff_lin * x + coeff_const - log(x)
    end

    lb = 0.0
    ub = 1-1e-12


    if f(ub) >= 0 
        throw(DomainError("Solution is between 1-1e-12 and 1, try with lower bounds"))
    end

    while true
        c = (lb+ub)/2
        v = f(c)
        if isapprox(v, 0, atol=1e-7)
            return c
        elseif sign(v) != sign(f(lb))    
            ub = c
        elseif sign(v) != sign(f(ub))
            lb = c
        end
    end
end


function getCpu(sn::MetaDiGraph, commodity::Commodity, path::Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}})
    delay_to_meet = commodity.max_delay - 2 * arc_path_delay(sn, path)
    if delay_to_meet < 0
        throw(DomainError(delay_to_meet, "The maximum delay has to be positive! Please select a shorter path"))
    end
    ro = 0
    #P(W <= delay_to_meet) = sla
    UB_START = 100
    lb = 0.01
    ub = UB_START
    mu = -1
    nu = nothing
    mus = []
    while true
        mu = (lb + ub) / 2
        push!(mus, mu)
        ro = getRo(mu, commodity.arr_mean, commodity.arr_coeffs)

        if ro <= 1
            try
                nu = getNu(mu, commodity.arr_mean, commodity.arr_sigma, commodity.arr_coeffs)
            # if we have a domain error it means mu is veeeeeery high
            catch e
                if isa(e, DomainError)
                    ub = mu
                    if length(mus) > 1 && mus[end] == mus[end-1]
                        UB_START = UB_START * .8
                        lb = 0.0
                        ub = UB_START
                        probs = []
                        nu = nothing
                        mu = -1
                    end
                    continue
                end
            end
            prob = 1 - nu * exp((-mu)*(1-nu) * delay_to_meet)
            if prob <= 0
                stop()
            end

            if prob <= commodity.sla + 0.001 && prob >= commodity.sla - 0.001 || (lb==ub)
                break
            end
            if prob > commodity.sla
                ub = mu
            elseif prob < commodity.sla
                lb = mu
            end
        else
            lb = mu
        end

        # if the bounds were too low (because mu does not move e.g. prob << commodity.sla), restart with a higher one
        if length(mus) > 1 && mus[end] == mus[end-1]
            UB_START = UB_START * 1.2
            lb = 0.0
            ub = UB_START
            probs = []
            nu = nothing
            mu = -1
        end
    end

    return ceil(mu * commodity.cpu_request)
end