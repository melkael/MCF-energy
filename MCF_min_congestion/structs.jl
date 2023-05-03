using MetaGraphs
using DataStructures

mutable struct Commodity
    id::Int64
    source::Int64
    dests::Vector{Int64}
    arr_mean::Vector{Float64}
    arr_sigma::Vector{Float64}
    arr_coeffs::Vector{Float64}
    sla::Float64
    max_delay::Float64
    cpu_request::Int64
    BW_demanded::Float64
    meta_destination::Int64
end

mutable struct Problem
    sn::MetaDiGraph   
    commodities::Vector{Commodity}
    paths
    model
    c_edge
    edge_contains_paths
    c_commodity
    c_node
    vars
    cpus
    artificial_variables
    map_vars_paths
    forbidden_arcs_constraints
    forbidden_arcs
    knapsackCoverCutsEdges
    knapsackCoverCutsNodes
    congestion
    c_cpu_cong
    c_BW_cong
end

mutable struct NodeBnB
    lb::Float64
    ub::Float64
    id::Int64
    parent::Union{NodeBnB, Nothing}
    problem::Union{Problem, Nothing}
    sibling::Union{NodeBnB, Nothing}
    hasBeenSolved::Bool
    newlyForbiddenArcs
end

mutable struct BnBTree
    sn
    commodities
    incumbent::Float64
    incumbent_solution
    lb::Float64
    node_queue::PriorityQueue{Int,Tuple{Float64, Int}}
    nodes::Dict{Int, NodeBnB}
    num_nodes::Int
    seen_nodes
    ub
    startTime
    timeout
    has_timedout
    threshold
    num_nodes_treated
    column_generation_has_broken
end