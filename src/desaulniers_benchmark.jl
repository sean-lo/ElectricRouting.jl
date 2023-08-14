include("utils.jl")
using DataStructures
using Printf

Base.@kwdef mutable struct PurePathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    explored::Bool
    served::Vector{Int}
    artificial::Bool = false
end

Base.copy(p::PurePathLabel) = PurePathLabel(
    p.cost,
    copy(p.nodes),
    copy(p.excesses),
    copy(p.slacks),
    p.time_mincharge,
    p.time_maxcharge,
    p.charge_mincharge,
    p.charge_maxcharge,
    p.explored,
    copy(p.served),
    p.artificial,
)

function compute_path_label_cost(
    p::PurePathLabel,
    data::EVRPData,
    graph::EVRPGraph, 
    M::Float64 = 1e10,
    ;
    verbose = false,
)
    if p.artificial
        return M
    elseif length(p.nodes) ≤ 1
        return 0
    end

    arcs = collect(zip(p.nodes[1:end-1], p.nodes[2:end]))
    cost = data.travel_cost_coeff * sum(graph.c[a...] for a in arcs)
    verbose && @printf("Path cost: \t\t%11.3f\n", cost)

    charging_cost = data.charge_cost_coeff * (sum(p.slacks) + sum(p.excesses))
    verbose && @printf("Charging cost: \t\t%11.3f\n", charging_cost)
    cost += charging_cost

    return cost
end


function compute_path_label_modified_cost(
    p::PurePathLabel,
    data::EVRPData,
    graph::EVRPGraph, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    verbose = false,
)
    reduced_cost = compute_path_label_cost(p, data, graph, verbose = verbose)

    service_cost = 0.0
    for (j, c) in enumerate(p.served)
        service_cost += (c * -ν[j])
    end
    verbose && @printf("Service cost: \t\t%11.3f\n", service_cost)
    reduced_cost += service_cost

    verbose && @printf("Starting depot cost: \t%11.3f\n", (- κ[p.nodes[1]]))
    reduced_cost = reduced_cost - κ[p.nodes[1]]

    verbose && @printf("Ending depot cost: \t%11.3f\n", (- μ[p.nodes[end]]))
    reduced_cost = reduced_cost - μ[p.nodes[end]]

    verbose && @printf("Total modified cost: \t%11.3f\n\n", reduced_cost)

    return reduced_cost
end

function add_pure_path_label_to_collection!(
    collection::SortedDict{
        NTuple{N, Int},
        PurePathLabel,
        Base.Order.ForwardOrdering,
    },
    key::NTuple{N, Int},
    path::PurePathLabel,
    ;
    verbose::Bool = false,
) where {N}
    added = true
    for (k, p) in pairs(collection)
        if p.cost ≤ path.cost
            if all(k .≤ key)
                added = false
                if verbose
                    println("$(key), $(path.cost) dominated by $(k), $(p.cost)")
                end
                break
            end
        end
        if path.cost ≤ p.cost
            if all(key .≤ k)
                if verbose
                    println("$(key), $(path.cost) dominates $(k), $(p.cost)")
                end
                pop!(collection, k)
            end
        end
    end
    if added
        if verbose
            println("$(key), $(path.cost) added!")
        end
        insert!(collection, key, path)
    end
    return added
end

function compute_new_pure_path(
    current_path::PurePathLabel,
    current_node::Int,
    next_node::Int,
    data::EVRPData,
    graph::EVRPGraph,
    α::Vector{Int},
    β::Vector{Int},
    modified_costs::Matrix{Float64},
)
    # feasibility checks
    # (1) battery
    excess = max(
        0, 
        graph.q[current_node,next_node] - current_path.charge_mincharge 
    )
    # (2) time windows
    if current_path.time_mincharge + excess + graph.t[current_node,next_node] > β[next_node]
        # println("$(current_path.time_mincharge), $excess, $(t[current_node,next_node]), $(β[next_node])")
        # println("not time windows feasible")
        return (false, nothing)
    end
    if current_path.time_mincharge + excess + graph.t[current_node,next_node] + graph.min_t[next_node] > graph.T
        return (false, nothing)
    end
    # (3) charge interval 
    if (
        (current_node in graph.N_charging_extra && excess > max(graph.B - current_path.charge_mincharge, 0))
        || 
        (!(current_node in graph.N_charging_extra) && excess > max(current_path.charge_maxcharge - current_path.charge_mincharge, 0))
    )
        # if current_node in graph.N_charging_extra
        #     println("$excess, $(B), $(current_path.charge_mincharge)")
        # else
        #     println("$excess, $(current_path.charge_maxcharge), $(current_path.charge_mincharge)")
        # end
        # println("not charge feasible")
        return (false, nothing)
    end

    new_path = copy(current_path)
    push!(new_path.nodes, next_node)
    if next_node in graph.N_customers
        new_path.served[next_node] += 1
    end

    push!(new_path.excesses, excess)
    new_path.time_mincharge = max(
        α[next_node],
        current_path.time_mincharge + graph.t[current_node,next_node] + excess
    )

    if current_node in graph.N_charging_extra
        slack = min(
            new_path.time_mincharge - (current_path.time_mincharge + graph.t[current_node,next_node] + excess),
            graph.B - (current_path.charge_mincharge + excess),
        )
        push!(new_path.slacks, slack)
        new_path.time_maxcharge = min(
            β[next_node],
            max(
                α[next_node],
                current_path.time_mincharge + (graph.B - current_path.charge_mincharge) + graph.t[current_node,next_node],
            )
        )
    else
        slack = min(
            new_path.time_mincharge - (current_path.time_mincharge + graph.t[current_node,next_node] + excess),
            current_path.charge_maxcharge - (current_path.charge_mincharge + excess),
        )
        push!(new_path.slacks, slack)
        new_path.time_maxcharge = min(
            β[next_node],
            max(
                α[next_node],
                current_path.time_maxcharge + graph.t[current_node,next_node],
            )
        )
    end

    new_path.charge_mincharge = (
        current_path.charge_mincharge 
        + excess 
        + slack
        - graph.q[current_node,next_node]
    )
    new_path.charge_maxcharge = (
        new_path.charge_mincharge 
        + new_path.time_maxcharge 
        - new_path.time_mincharge
    )

    new_path.cost += modified_costs[current_node,next_node]
    new_path.cost += data.charge_cost_coeff * (slack + excess)

    return (true, new_path)
end

function find_nondominated_paths(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    single_service::Bool = false,
    time_windows::Bool = true,
    check_customers::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    keylen = check_customers ? graph.n_customers + 3 : 3
    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{keylen, Int}, 
                PurePathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    if check_customers
        # label key here has the following fields:
        # 1) current minimum time T_i(min)
        # 2) negative of current max charge -B_i(max)
        # 3) difference between min time and min charge, T_i(min) - B_i(min)
        # 4) if applicable, whether i-th customer served
        key = (0, -graph.B, -graph.B, zeros(Int, graph.n_customers)...)
    else
        key = (0, -graph.B, -graph.B)
    end

    unexplored_states = SortedSet{NTuple{keylen + 2, Int}}()
    for depot in graph.N_depots
        pure_path_labels[depot][depot][key] = PurePathLabel(
            0.0,
            [depot],
            Int[],
            Int[],
            0,
            0,
            graph.B,
            graph.B,
            false,
            zeros(Int, graph.n_customers),
            false,
        )
        push!(unexplored_states, (key..., depot, depot))
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(pure_path_labels[starting_node][current_node]))
            continue
        end
        current_path = pure_path_labels[starting_node][current_node][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers
                # single-service requirement
                if (
                    single_service 
                    && current_path.served[next_node] > 0
                )
                    # println("already served $next_node")
                    continue
                end
                # Preventing customer 2-cycles (Christofides)
                if christofides 
                    if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == next_node
                        continue
                    end
                end
            end

            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            if !feasible
                continue
            end

            # add new_path to collection
            if check_customers
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge, 
                    new_path.served...,
                )
            else
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge
                )
            end
            added = add_pure_path_label_to_collection!(
                pure_path_labels[starting_node][next_node], 
                new_key, new_path, 
                ;
                verbose = false,
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for depot in graph.N_depots
        for path in values(pure_path_labels[depot][depot])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end
    
    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(pure_path_labels[starting_node][end_node])
                path.cost = path.cost - κ[starting_node] - μ[end_node]
            end
        end
    end

    return pure_path_labels
end

function find_nondominated_paths_ngroute(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    time_windows::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                Tuple{Vararg{Int}}, 
                SortedDict{
                    NTuple{3, Int}, 
                    PurePathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{5, Int}}()
    for depot in graph.N_depots
        key = (0, -graph.B, -graph.B)
        set = (depot,) # NG
        pure_path_labels[depot][depot][set] = SortedDict{
            NTuple{3, Int},
            PurePathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => PurePathLabel(
                0.0,
                [depot],
                Int[],
                Int[],
                0,
                0,
                graph.B,
                graph.B,
                false,
                zeros(Int, graph.n_customers),
                false,
            ),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                depot, # starting_node
                depot, # current_node 
            )
        )
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        for current_set in keys(pure_path_labels[starting_node][current_node])
            if !(current_key in keys(pure_path_labels[starting_node][current_node][current_set]))
                continue
            end
            current_path = pure_path_labels[starting_node][current_node][current_set][current_key]
            for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
                if next_node in current_set
                    # if next_node is a customer not yet visited, proceed
                    # only if one can extend current_subpath along next_node according to ng-route rules
                    continue
                end
                if next_node in graph.N_customers
                    # Preventing customer 2-cycles (Christofides)
                    if christofides 
                        if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == next_node
                            continue
                        end
                    end
                end

                (feasible, new_path) = compute_new_pure_path(
                    current_path, 
                    current_node, next_node, 
                    data, graph,
                    α, β, modified_costs,
                )
                if !feasible
                    continue
                end
                
                new_set = ngroute_create_set(neighborhoods, current_set, next_node)
                if !(new_set in keys(pure_path_labels[starting_node][next_node]))
                    pure_path_labels[starting_node][next_node][new_set] = SortedDict{
                        NTuple{3, Int},
                        PurePathLabel,
                        Base.Order.ForwardOrdering,
                    }(Base.Order.ForwardOrdering())
                end
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge
                )
                added = add_pure_path_label_to_collection!(
                    pure_path_labels[starting_node][next_node][new_set], 
                    new_key, new_path, 
                    ;
                    verbose = false,
                )
                if added && !(next_node in graph.N_depots)
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end
    
    for depot in graph.N_depots
        for set in keys(pure_path_labels[depot][depot])
            for path in values(pure_path_labels[depot][depot][set])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                    path.excesses = [0]
                    path.slacks = [0]
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for set in keys(pure_path_labels[starting_node][end_node])
                for path in values(pure_path_labels[starting_node][end_node][set])
                    path.cost = path.cost - κ[starting_node] - μ[end_node]
                end
            end
        end
    end

    return pure_path_labels
end

function find_nondominated_paths_ngroute_1(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64},
    ;
    time_windows::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => Dict{
                    Tuple{Vararg{Int}}, 
                    SortedDict{
                        NTuple{3, Int}, 
                        PurePathLabel,
                        Base.Order.ForwardOrdering,
                    },
                }()
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{6, Int}}()
    for depot in graph.N_depots
        key = (0, -graph.B, -graph.B)
        set = (depot,) # NG
        pure_path_labels[depot][depot][depot][set] = SortedDict{
            NTuple{3, Int},
            PurePathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => PurePathLabel(
                0.0,
                [depot],
                Int[],
                Int[],
                0,
                0,
                graph.B,
                graph.B,
                false,
                zeros(Int, graph.n_customers),
                false,
            ),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                depot, # starting_node
                depot, # prev_node
                depot, # current_node 
            )
        )
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-3]
        for current_set in keys(pure_path_labels[starting_node][current_node][prev_node])
            if !(current_key in keys(pure_path_labels[starting_node][current_node][prev_node][current_set]))
                continue
            end
            current_path = pure_path_labels[starting_node][current_node][prev_node][current_set][current_key]
            for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
                if next_node in current_set
                    # if next_node is a customer not yet visited, proceed
                    # only if one can extend current_subpath along next_node according to ng-route rules
                    continue
                end
                if next_node in graph.N_customers
                    # Preventing customer 2-cycles (Christofides)
                    if christofides 
                        if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == prev_node == next_node
                            continue
                        end
                    end
                end

                
                (feasible, new_path) = compute_new_pure_path(
                    current_path, 
                    current_node, next_node, 
                    data, graph,
                    α, β, modified_costs,
                )
                if !feasible
                    continue
                end

                # if !(current_node in keys(pure_path_labels[starting_node][next_node]))
                #     pure_path_labels[starting_node][next_node][current_node] = Dict{
                #         Tuple{Vararg{Int}},
                #         SortedDict{
                #             NTuple{3, Int},
                #             PurePathLabel,
                #             Base.Order.ForwardOrdering,
                #         },
                #     }()
                # end
                new_set = ngroute_create_set(neighborhoods, current_set, next_node)
                if !(new_set in keys(pure_path_labels[starting_node][next_node][current_node]))
                    pure_path_labels[starting_node][next_node][current_node][new_set] = SortedDict{
                        NTuple{3, Int},
                        PurePathLabel,
                        Base.Order.ForwardOrdering,
                    }(Base.Order.ForwardOrdering())
                end
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge
                )
                added = add_pure_path_label_to_collection!(
                    pure_path_labels[starting_node][next_node][current_node][new_set], 
                    new_key, new_path, 
                    ;
                    verbose = false,
                )
                if added && !(next_node in graph.N_depots)
                    new_state = (new_key..., starting_node, current_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end
    
    for depot in graph.N_depots
        for prev_node in keys(pure_path_labels[depot][depot])
            for set in keys(pure_path_labels[depot][depot][prev_node])
                for path in values(pure_path_labels[depot][depot][prev_node][set])
                    if length(path.nodes) == 1
                        path.nodes = [depot, depot]
                        path.excesses = [0]
                        path.slacks = [0]
                    end
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for prev_node in keys(pure_path_labels[starting_node][end_node])
                for set in keys(pure_path_labels[starting_node][end_node][prev_node])
                    for path in values(pure_path_labels[starting_node][end_node][prev_node][set])
                        path.cost = path.cost - κ[starting_node] - μ[end_node]
                    end
                end
            end
        end
    end

    return pure_path_labels
end

function find_nondominated_paths_ngroute_sigma(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64},
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    time_windows::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    σ_costs = compute_WSR3_sigma_2costs(σ, graph)  

    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int}, 
                Dict{
                    Tuple{Vararg{Int}}, 
                    SortedDict{
                        NTuple{3, Int}, 
                        PurePathLabel,
                        Base.Order.ForwardOrdering,
                    },
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{7, Int}}()
    for depot in graph.N_depots
        key = (0, -graph.B, -graph.B)
        set = (depot,) # NG
        pure_path_labels[depot][depot][(depot, depot)] = Dict(
            set => SortedDict{
                NTuple{3, Int},
                PurePathLabel,
            }(
                # Base.Order.ForwardOrdering(),
                key => PurePathLabel(
                    0.0,
                    [depot],
                    Int[],
                    Int[],
                    0,
                    0,
                    graph.B,
                    graph.B,
                    false,
                    zeros(Int, graph.n_customers),
                    false,
                ),
            )
        )
        push!(
            unexplored_states, 
            (
                key..., 
                depot, # starting_node
                depot, # prev_prev_node
                depot, # prev_node
                depot, # current_node 
            )
        )
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        prev_prev_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-4]
        for current_set in keys(pure_path_labels[starting_node][current_node][(prev_prev_node, prev_node)])
            if !(current_key in keys(pure_path_labels[starting_node][current_node][(prev_prev_node, prev_node)][current_set]))
                continue
            end
            current_path = pure_path_labels[starting_node][current_node][(prev_prev_node, prev_node)][current_set][current_key]
            for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
                if next_node in current_set
                    # if next_node is a customer not yet visited, proceed
                    # only if one can extend current_subpath along next_node according to ng-route rules
                    continue
                end
                if next_node in graph.N_customers
                    # Preventing customer 2-cycles (Christofides)
                    if christofides 
                        if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == next_node
                            continue
                        end
                    end
                end

                
                (feasible, new_path) = compute_new_pure_path(
                    current_path, 
                    current_node, next_node, 
                    data, graph,
                    α, β, modified_costs,
                )
                if !feasible
                    continue
                end

                new_path.cost += σ_costs[(prev_prev_node, prev_node, current_node, next_node)]

                if !((prev_node, current_node) in keys(pure_path_labels[starting_node][next_node]))
                    pure_path_labels[starting_node][next_node][(prev_node, current_node)] = Dict{
                        Tuple{Vararg{Int}},
                        SortedDict{
                            NTuple{3, Int},
                            PurePathLabel,
                            Base.Order.ForwardOrdering,
                        },
                    }()
                end
                new_set = ngroute_create_set(neighborhoods, current_set, next_node)
                if !(new_set in keys(pure_path_labels[starting_node][next_node][(prev_node, current_node)]))
                    pure_path_labels[starting_node][next_node][(prev_node, current_node)][new_set] = SortedDict{
                        NTuple{3, Int},
                        PurePathLabel,
                        Base.Order.ForwardOrdering,
                    }(Base.Order.ForwardOrdering())
                end
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge
                )
                added = add_pure_path_label_to_collection!(
                    pure_path_labels[starting_node][next_node][(prev_node, current_node)][new_set], 
                    new_key, new_path, 
                    ;
                    verbose = false,
                )
                if added && !(next_node in graph.N_depots)
                    new_state = (new_key..., starting_node, prev_node, current_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end
    
    for depot in graph.N_depots
        for (prev_prev_node, prev_node) in keys(pure_path_labels[depot][depot])
            for set in keys(pure_path_labels[depot][depot][(prev_prev_node, prev_node)])
                for path in values(pure_path_labels[depot][depot][(prev_prev_node, prev_node)][set])
                    if length(path.nodes) == 1
                        path.nodes = [depot, depot]
                        path.excesses = [0]
                        path.slacks = [0]
                    end
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for (prev_prev_node, prev_node) in keys(pure_path_labels[starting_node][end_node])
                for set in keys(pure_path_labels[starting_node][end_node][(prev_prev_node, prev_node)])
                    for path in values(pure_path_labels[starting_node][end_node][(prev_prev_node, prev_node)][set])
                        path.cost = path.cost - κ[starting_node] - μ[end_node]
                    end
                end
            end
        end
    end

    return pure_path_labels
end



function find_nondominated_paths_ngroute_alt(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    time_windows::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{graph.n_nodes_extra + 3, Int}, 
                PurePathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 5, Int}}()
    for depot in graph.N_depots
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[depot] = 1
        key = (0, -graph.B, -graph.B, node_labels...)
        pure_path_labels[depot][depot][key] = PurePathLabel(
            0.0,
            [depot],
            Int[],
            Int[],
            0,
            0,
            graph.B,
            graph.B,
            false,
            zeros(Int, graph.n_customers),
            false,
        )
        push!(unexplored_states, (key..., depot, depot))
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(pure_path_labels[starting_node][current_node]))
            continue
        end
        current_set = state[4:end-2]
        current_path = pure_path_labels[starting_node][current_node][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            if next_node in graph.N_customers
                # Preventing customer 2-cycles (Christofides)
                if christofides 
                    if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == next_node
                        continue
                    end
                end
            end

            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            if !feasible
                continue
            end

            new_set = ngroute_create_set_alt(neighborhoods, collect(current_set), next_node)
            new_key = (
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
                new_set...,
            )
            added = add_pure_path_label_to_collection!(
                pure_path_labels[starting_node][next_node], 
                new_key, new_path, 
                ;
                verbose = false,
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for depot in graph.N_depots
        for path in values(pure_path_labels[depot][depot])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end
    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(pure_path_labels[starting_node][end_node])
                path.cost = path.cost - κ[starting_node] - μ[end_node]
            end
        end
    end

    return pure_path_labels
end



function find_nondominated_paths_ngroute_alt_1(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    time_windows::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => SortedDict{
                    NTuple{graph.n_nodes_extra + 3, Int}, 
                    PurePathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 6, Int}}()
    for depot in graph.N_depots
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[depot] = 1
        key = (0, -graph.B, -graph.B, node_labels...)
        pure_path_labels[depot][depot][depot][key] = PurePathLabel(
            0.0,
            [depot],
            Int[],
            Int[],
            0,
            0,
            graph.B,
            graph.B,
            false,
            zeros(Int, graph.n_customers),
            false,
        )
        push!(
            unexplored_states, 
            (
                key..., 
                depot, # starting_node
                depot, # prev_node
                depot, # current_node
            )
        )
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-3]
        current_set = state[4:end-3]
        current_path = pure_path_labels[starting_node][current_node][prev_node][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            if next_node in graph.N_customers
                # Preventing customer 2-cycles (Christofides)
                if christofides 
                    if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == prev_node == next_node
                        continue
                    end
                end
            end

            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            if !feasible
                continue
            end

            # if !(current_node in keys(pure_path_labels[starting_node][next_node]))
            #     pure_path_labels[starting_node][next_node][current_node] = SortedDict{
            #         NTuple{graph.n_nodes_extra + 3, Int}, 
            #         PurePathLabel,
            #         Base.Order.ForwardOrdering,
            #     }(Base.Order.ForwardOrdering())
            # end
            new_set = ngroute_create_set_alt(neighborhoods, collect(current_set), next_node)
            new_key = (
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
                new_set...,
            )
            added = add_pure_path_label_to_collection!(
                pure_path_labels[starting_node][next_node][current_node], 
                new_key, new_path, 
                ;
                verbose = false,
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for depot in graph.N_depots
        for prev_node in keys(pure_path_labels[depot][depot])
            for path in values(pure_path_labels[depot][depot][prev_node])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                    path.excesses = [0]
                    path.slacks = [0]
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for prev_node in keys(pure_path_labels[starting_node][end_node])
                for path in values(pure_path_labels[starting_node][end_node][prev_node])
                    path.cost = path.cost - κ[starting_node] - μ[end_node]
                end
            end
        end
    end
    
    return pure_path_labels
end


function find_nondominated_paths_ngroute_alt_sigma(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    time_windows::Bool = true,
    christofides::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    σ_costs = compute_WSR3_sigma_2costs(σ, graph)

    pure_path_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int},
                SortedDict{
                    NTuple{graph.n_nodes_extra + 3, Int}, 
                    PurePathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 7, Int}}()
    for depot in graph.N_depots
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[depot] = 1
        key = (0, -graph.B, -graph.B, node_labels...)
        pure_path_labels[depot][depot][(depot, depot)] = SortedDict{
            NTuple{graph.n_nodes_extra + 3, Int}, 
            PurePathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => PurePathLabel(
                0.0,
                [depot],
                Int[],
                Int[],
                0,
                0,
                graph.B,
                graph.B,
                false,
                zeros(Int, graph.n_customers),
                false,
            ),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                depot, # starting_node
                depot, # prev_prev_node
                depot, # prev_node
                depot, # current_node
            )
        )
    end

    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        prev_prev_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-4]
        current_set = state[4:end-4]
        current_path = pure_path_labels[starting_node][current_node][(prev_prev_node, prev_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            if next_node in graph.N_customers
                # Preventing customer 2-cycles (Christofides)
                if christofides 
                    if length(current_path.nodes) ≥ 2 && current_path.nodes[end-1] == next_node
                        continue
                    end
                end
            end

            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            if !feasible
                continue
            end

            new_path.cost += σ_costs[(prev_prev_node, prev_node, current_node, next_node)]

            if !((prev_node, current_node) in keys(pure_path_labels[starting_node][next_node]))
                pure_path_labels[starting_node][next_node][(prev_node, current_node)] = SortedDict{
                    NTuple{graph.n_nodes_extra + 3, Int}, 
                    PurePathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_set = ngroute_create_set_alt(neighborhoods, collect(current_set), next_node)
            new_key = (
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
                new_set...,
            )
            added = add_pure_path_label_to_collection!(
                pure_path_labels[starting_node][next_node][(prev_node, current_node)], 
                new_key, new_path, 
                ;
                verbose = false,
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, prev_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for depot in graph.N_depots
        for (prev_prev_node, prev_node) in keys(pure_path_labels[depot][depot])
            for path in values(pure_path_labels[depot][depot][(prev_prev_node, prev_node)])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                    path.excesses = [0]
                    path.slacks = [0]
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(keys(pure_path_labels[starting_node]), graph.N_depots)
            delete!(pure_path_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for (prev_prev_node, prev_node) in keys(pure_path_labels[starting_node][end_node])
                for path in values(pure_path_labels[starting_node][end_node][(prev_prev_node, prev_node)])
                    path.cost = path.cost - κ[starting_node] - μ[end_node]
                end
            end
        end
    end

    return pure_path_labels
end

unwrap_pure_path_labels(p) = PurePathLabel[p]

function unwrap_pure_path_labels(d::AbstractDict)
    u = PurePathLabel[]
    for v in values(d)
        append!(u, unwrap_pure_path_labels(v))
    end
    return u
end

function get_negative_pure_path_labels_from_pure_path_labels(
    pure_path_labels::Dict{
        Int, 
        Dict{Int, T},
    },
) where {T <: AbstractDict}
    # pure_path_labels_sd_list = []
    
    return PurePathLabel[
        pure_path_label
        for pure_path_label in unwrap_pure_path_labels(pure_path_labels)
            if pure_path_label.cost < -1e-6
    ]
end

# function get_negative_pure_path_labels_from_pure_path_labels(
#     pure_path_labels::Dict{
#         Int, 
#         Dict{Int, T},
#     },
# ) where {T <: AbstractDict}
#     pure_path_labels_sd_list = []

#     return PurePathLabel[
#         pure_path_label
#         for pure_path_label in unwrap_pure_path_labels(pure_path_labels)
#             if pure_path_label.cost < -1e-6
#     ]
# end

function subproblem_iteration_benchmark(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    time_windows::Bool = false,
    path_single_service::Bool = true,
    path_check_customers::Bool = true,
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)
    start_time = time()
    if ngroute && !ngroute_alt
        if length(σ) == 0
            pure_path_labels_result = @timed find_nondominated_paths_ngroute_1(
                data, graph, neighborhoods, κ, μ, ν,
                ;
                time_windows = time_windows,
                christofides = christofides,
                time_limit = time_limit - (time() - start_time),
            )
        else
            pure_path_labels_result = @timed find_nondominated_paths_ngroute_sigma(
                data, graph, neighborhoods, κ, μ, ν, σ,
                ;
                time_windows = time_windows,
                christofides = christofides,
                time_limit = time_limit - (time() - start_time),
            )
        end
    elseif ngroute && ngroute_alt
        if length(σ) == 0
            pure_path_labels_result = @timed find_nondominated_paths_ngroute_alt_1(
                data, graph, neighborhoods, κ, μ, ν,
                ;
                time_windows = time_windows,
                christofides = christofides,
                time_limit = time_limit - (time() - start_time),
            )
        else
            pure_path_labels_result = @timed find_nondominated_paths_ngroute_alt_sigma(
                data, graph, neighborhoods, κ, μ, ν, σ,
                ;
                time_windows = time_windows,
                christofides = christofides,
                time_limit = time_limit - (time() - start_time),
            )
        end
    else
        pure_path_labels_result = @timed find_nondominated_paths(
            data, graph, κ, μ, ν,
            ;
            time_windows = time_windows, 
            single_service = path_single_service, 
            check_customers = path_check_customers,
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    end
    pure_path_labels_time = pure_path_labels_result.time
    negative_pure_path_labels = get_negative_pure_path_labels_from_pure_path_labels( pure_path_labels_result.value)
    negative_pure_path_labels_count = length(negative_pure_path_labels)
    return (negative_pure_path_labels, negative_pure_path_labels_count, pure_path_labels_time)
end

