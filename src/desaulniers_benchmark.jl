include("utils.jl")

Base.@kwdef mutable struct PurePathLabel <: Label
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
        return (false, current_path)
    end
    if current_path.time_mincharge + excess + graph.t[current_node,next_node] + graph.min_t[next_node] > graph.T
        return (false, current_path)
    end
    # (3) charge interval 
    if (
        (current_node in graph.N_charging && excess > max(graph.B - current_path.charge_mincharge, 0))
        || 
        (!(current_node in graph.N_charging) && excess > max(current_path.charge_maxcharge - current_path.charge_mincharge, 0))
    )
        # if current_node in graph.N_charging
        #     println("$excess, $(B), $(current_path.charge_mincharge)")
        # else
        #     println("$excess, $(current_path.charge_maxcharge), $(current_path.charge_mincharge)")
        # end
        # println("not charge feasible")
        return (false, current_path)
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

    if current_node in graph.N_charging
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

    return (true, new_path)
end


function find_nondominated_paths(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    α::Vector{Int},
    β::Vector{Int},
    ;
    elementary::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    if elementary
        # label key here has the following fields:
        # 0) reduced cost
        # 1) current minimum time T_i(min)
        # 2) negative of current max charge -B_i(max)
        # 3) difference between min time and min charge, T_i(min) - B_i(min)
        # 4) if applicable, whether i-th customer served
        key = (0.0, 0, -graph.B, -graph.B, falses(graph.n_customers),)
        pure_path_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int, Int, BitVector}, 
                PurePathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_nodes
        )
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Int, BitVector, Int, Int}}()
    else
        key = (0.0, 0, -graph.B, -graph.B)
        pure_path_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int, Int}, 
                PurePathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_nodes
        )
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Int, Int, Int}}()
    end

    for depot in graph.N_depots
        pure_path_labels[(depot, depot)][key] = PurePathLabel(
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

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(pure_path_labels[(starting_node, current_node)]))
            continue
        end
        current_path = pure_path_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if (
                next_node in graph.N_customers 
                && elementary 
                && current_path.served[next_node] > 0
            )
                # single-service requirement
                # println("already served $next_node")
                continue
            end

            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            !feasible && continue

            # add new_path to collection
            if elementary
                new_key = (
                    new_path.cost,
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge, 
                    BitVector(new_path.served),
                )
            else
                new_key = (
                    new_path.cost,
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge
                )
            end
            added = add_label_to_collection!(
                pure_path_labels[(starting_node, next_node)], 
                new_key, new_path, 
                ;
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for depot in graph.N_depots
        for path in values(pure_path_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end
    
    for starting_node in graph.N_depots
        for end_node in setdiff(graph.N_nodes, graph.N_depots)
            delete!(pure_path_labels, (starting_node, end_node))
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(pure_path_labels[(starting_node, end_node)])
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
    α::Vector{Int},
    β::Vector{Int},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    pure_path_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, Int, Int, BitVector}, 
            PurePathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) current minimum time T_i(min)
        # 2) negative of current max charge -B_i(max)
        # 3) difference between min time and min charge, T_i(min) - B_i(min)
        # 4) whether i-th node is in forward ng-set
        key = (0.0, 0, -graph.B, -graph.B, ng_resources)
        pure_path_labels[(depot, depot)][key] = PurePathLabel(
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
                depot, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(pure_path_labels[(starting_node, current_node)]))
            continue
        end
        (_, _, _, _, current_ng_resources) = current_key
        current_path = pure_path_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_ng_resources) = ngroute_check_create_fset(
                neighborhoods, current_ng_resources, next_node,
            )
            !feasible && continue
            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            !feasible && continue

            new_key = (
                new_path.cost,
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
                new_ng_resources,
            )
            added = add_label_to_collection!(
                pure_path_labels[(starting_node, next_node)], 
                new_key, new_path, 
                ;
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for depot in graph.N_depots
        for path in values(pure_path_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(graph.N_nodes, graph.N_depots)
            delete!(pure_path_labels, (starting_node, end_node))
        end
    end
    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(pure_path_labels[(starting_node, end_node)])
                path.cost = path.cost - κ[starting_node] - μ[end_node]
            end
        end
    end

    return pure_path_labels
end



function find_nondominated_paths_ngroute_lambda(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{NTuple{3, Int}, Float64},
    α::Vector{Int},
    β::Vector{Int},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    λvals, λcust = prepare_lambda(λ, graph.n_nodes)

    pure_path_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, BitVector, Int, Int, Int, BitVector}, 
            PurePathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        λ_labels = falses(length(λ))
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) current minimum time T_i(min)
        # 3) negative of current max charge -B_i(max)
        # 4) difference between min time and min charge, T_i(min) - B_i(min)
        # 5) whether i-th node is in forward ng-set
        key = (0.0, λ_labels, 0, -graph.B, -graph.B, ng_resources,)
        pure_path_labels[(depot, depot)][key] = PurePathLabel(
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
                depot, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(pure_path_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_λ_labels, 
            _, _, _, current_ng_resources,
        ) = current_key
        current_path = pure_path_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_ng_resources) = ngroute_check_create_fset(
                neighborhoods, current_ng_resources, next_node,
            )
            !feasible && continue
            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            !feasible && continue

            (new_λ_labels, λ_cost) = compute_new_lambda_labels_cost(
                next_node, current_λ_labels, λvals, λcust,
            )
            new_path.cost += λ_cost

            new_key = (
                new_path.cost,
                new_λ_labels,
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
                new_ng_resources,
            )
            added = add_label_to_collection_cuts!(
                pure_path_labels[(starting_node, next_node)], 
                new_key, new_path, λvals,
                ;
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for depot in graph.N_depots
        for path in values(pure_path_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(graph.N_nodes, graph.N_depots)
            delete!(pure_path_labels, (starting_node, end_node))
        end
    end
    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(pure_path_labels[(starting_node, end_node)])
                path.cost = path.cost - κ[starting_node] - μ[end_node]
            end
        end
    end

    return pure_path_labels
end


function find_nondominated_paths_ngroute_lambda_lmSR3(
    data::EVRPData, 
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}, Float64},
    α::Vector{Int},
    β::Vector{Int},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)

    pure_path_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, BitVector, Int, Int, Int, BitVector,}, 
            PurePathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        λ_labels = falses(length(λ))
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) current minimum time T_i(min)
        # 3) negative of current max charge -B_i(max)
        # 4) difference between min time and min charge, T_i(min) - B_i(min)
        # 5) whether i-th node is in forward ng-set
        key = (0.0,  λ_labels, 0, -graph.B, -graph.B, ng_resources,)
        pure_path_labels[(depot, depot)][key] = PurePathLabel(
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
                depot, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(pure_path_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_λ_labels, 
            _, _, _, current_ng_resources,
        ) = current_key
        current_path = pure_path_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_ng_resources) = ngroute_check_create_fset(
                neighborhoods, current_ng_resources, next_node,
            )
            !feasible && continue
            (feasible, new_path) = compute_new_pure_path(
                current_path, 
                current_node, next_node, 
                data, graph,
                α, β, modified_costs,
            )
            !feasible && continue

            (new_λ_labels, λ_cost) = compute_lambda_flabels_cost_lmSR3(
                next_node, current_λ_labels, λvals, λcust, λmemory,
            )
            new_path.cost += λ_cost

            new_key = (
                new_path.cost,
                new_λ_labels,
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
                new_ng_resources,
            )
            added = add_label_to_collection_cuts!(
                pure_path_labels[(starting_node, next_node)], 
                new_key, new_path, λvals,
                ;
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for depot in graph.N_depots
        for path in values(pure_path_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(graph.N_nodes, graph.N_depots)
            delete!(pure_path_labels, (starting_node, end_node))
        end
    end
    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(pure_path_labels[(starting_node, end_node)])
                path.cost = path.cost - κ[starting_node] - μ[end_node]
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
        NTuple{2, Int}, 
        T,
    },
) where {T <: AbstractDict}
    return PurePathLabel[
        pure_path_label
        for pure_path_label in unwrap_pure_path_labels(pure_path_labels)
            if pure_path_label.cost < -1e-6
    ]
end

function subproblem_iteration_benchmark(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{T, Float64},
    ;
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    time_windows::Bool = false,
    elementary::Bool = true,
    time_limit::Float64 = Inf,
) where {T}
    start_time = time()
    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes)
        β = fill(graph.T, graph.n_nodes)
    end

    if ngroute
        if length(λ) == 0
            pure_path_labels_result = @timed find_nondominated_paths_ngroute(
                data, graph, neighborhoods, κ, μ, ν, α, β,
                ;
                time_limit = time_limit - (time() - start_time),
            )
        else
            if keytype(λ) == NTuple{3, Int}
                pure_path_labels_result = @timed find_nondominated_paths_ngroute_lambda(
                    data, graph, neighborhoods, κ, μ, ν, λ, α, β,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            elseif keytype(λ) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
                pure_path_labels_result = @timed find_nondominated_paths_ngroute_lambda_lmSR3(
                    data, graph, neighborhoods, κ, μ, ν, λ, α, β,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            else
                error("Unrecognized key type for λ: $(keytype(λ))")
            end
        end
    else
        pure_path_labels_result = @timed find_nondominated_paths(
            data, graph, κ, μ, ν, α, β,
            ;
            elementary = elementary,
            time_limit = time_limit - (time() - start_time),
        )
    end
    pure_path_labels_time = pure_path_labels_result.time
    negative_pure_path_labels = get_negative_pure_path_labels_from_pure_path_labels(pure_path_labels_result.value)
    negative_pure_path_labels_count = length(negative_pure_path_labels)
    return (negative_pure_path_labels, negative_pure_path_labels_count, pure_path_labels_time)
end

