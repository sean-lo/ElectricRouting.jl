using JuMP
using Gurobi
using Suppressor
using DataStructures
using Printf

include("subpath_stitching.jl")
include("path_formulation.jl")
include("utils.jl")

function compute_new_subpath_nocharge(
    current_subpath::BaseSubpathLabel,
    graph::EVRPGraph,
    current_node::Int,
    next_node::Int,
    modified_costs::Matrix{Float64},
    T_heuristic::Int,
)

    # time feasibility
    new_time_taken = current_subpath.time_taken + graph.t[current_node, next_node]
    if new_time_taken + graph.min_t[next_node] > T_heuristic
        return (false, current_subpath)
    end

    new_subpath = copy(current_subpath)
    new_subpath.time_taken = new_time_taken
    new_subpath.cost += modified_costs[current_node, next_node]
    push!(new_subpath.nodes, next_node)
    if next_node in graph.N_customers
        new_subpath.served[next_node] += 1
    end
    return (true, new_subpath)
end

function find_nondominated_paths_nocharge(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    elementary::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    keylen = elementary ? graph.n_customers + 1 : 1
    base_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Vararg{Int, keylen}}, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in union(graph.N_depots, graph.N_customers)
    )

    if elementary
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time taken
        # 2) if applicable, number of times i-th customer served
        key = (0.0, 0, zeros(Int, graph.n_customers))
    else
        key = (0.0, 0,)
    end

    unexplored_states = SortedSet{Tuple{Float64, Vararg{Int, keylen + 2}}}()
    for node in graph.N_depots
        base_labels[(node, node)][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(unexplored_states, (key..., node, node))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(base_labels[(starting_node, current_node)]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging,
        )
            if (
                next_node in graph.N_customers
                && elementary
                && current_path.served[next_node] > 0
            )
                # single-service requirement
                # println("already served $next_node")
                continue
            end

            (feasible, new_subpath) = compute_new_subpath_nocharge(
                current_subpath, graph,
                current_node, next_node, 
                modified_costs, T_heuristic,
            )
            !feasible && continue

            if elementary
                new_key = (
                    new_subpath.cost,
                    new_subpath.time_taken,
                    new_subpath.served...,
                )            
                added = add_subpath_longlabel_to_collection_nodelabels!(
                    base_labels[(starting_node, next_node)], 
                    new_key, new_subpath,
                    ;
                )
            else
                new_key = (
                    new_subpath.cost,
                    new_subpath.time_taken,
                )
                added = add_subpath_longlabel_to_collection!(
                    base_labels[(starting_node, next_node)], 
                    new_key, new_subpath,
                    ;
                )
            end

            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for depot in graph.N_depots
        for path in values(base_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    base_labels = Dict(
        (starting_node, current_node) => Dict{
            BitVector, 
            SortedDict{
                Tuple{Float64, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        }()
        for starting_node in graph.N_depots,
            current_node in union(graph.N_depots, graph.N_customers)
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, Int, Int}}()
    for node in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        key = (0.0, 0,)
        base_labels[(node, node)][node_labels] = SortedDict{
            Tuple{Float64, Int}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            )
        )
        push!(
            unexplored_states, 
            (
                key...,
                node_labels,
                node, # starting_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[end-2]
        current_key = state[1:end-3]
        if !(current_key in keys(base_labels[(starting_node, current_node)][current_set]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_set][current_key]
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging,
        )
            (feasible, new_set) = ngroute_check_create_fset(
                neighborhoods, current_set, next_node,
            )
            !feasible && continue

            (feasible, new_subpath) = compute_new_subpath_nocharge(
                current_subpath, graph,
                current_node, next_node, 
                modified_costs, T_heuristic,
            )
            !feasible && continue

            if !(new_set in keys(base_labels[(starting_node, next_node)]))
                base_labels[(starting_node, next_node)][new_set] = SortedDict{
                    Tuple{Float64, Int},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (
                new_subpath.cost,
                new_subpath.time_taken,
            )
            added = add_subpath_longlabel_to_collection!(
                base_labels[(starting_node, next_node)][new_set],
                new_key, new_subpath,
                ;
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set, starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for depot in graph.N_depots
        for set in keys(base_labels[(depot, depot)])
            for path in values(base_labels[(depot, depot)][set])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute_alt(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    
    base_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, BitVector}, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in union(graph.N_depots, graph.N_customers)
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, Int, Int}}()
    for node in graph.N_depots
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        key = (0.0, 0, node_labels,)
        base_labels[(node, node)][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node, # starting_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[end-2]
        current_key = state[1:end-2]
        if !(current_key in keys(base_labels[(starting_node, current_node)]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging,
        )
            (feasible, new_set) = ngroute_check_create_fset(
                neighborhoods, current_set, next_node,
            )
            !feasible && continue

            (feasible, new_subpath) = compute_new_subpath_nocharge(
                current_subpath, graph,
                current_node, next_node, 
                modified_costs, T_heuristic,
            )
            !feasible && continue

            new_key = (
                new_subpath.cost,
                new_subpath.time_taken,
                new_set,
            )
            added = add_subpath_longlabel_to_collection_ngroute_alt!(
                base_labels[(starting_node, next_node)],
                new_key, new_subpath,
                ;
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for depot in graph.N_depots
        for path in values(base_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{NTuple{3, Int}, Float64},
    T_heuristic::Int,
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    λvals, λcust = prepare_lambda(λ, graph.n_nodes)

    base_labels = Dict(
        (starting_node, current_node) => Dict{
            BitVector, 
            SortedDict{
                Tuple{Float64, Int, BitVector}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        }()
        for starting_node in graph.N_depots,
            current_node in union(graph.N_depots, graph.N_customers)
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, BitVector, Int, Int}}()
    for node in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        λ_labels = falses(length(λ))
        key = (0.0, 0, λ_labels)
        base_labels[(node, node)][node_labels] = SortedDict{
            Tuple{Float64, Int, BitVector}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            )
        )
        push!(
            unexplored_states, 
            (
                key...,
                node_labels,
                node, # starting_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[end-2]
        current_key = state[1:end-3]
        current_λ_labels = state[end-3]
        if !(current_key in keys(base_labels[(starting_node, current_node)][current_set]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_set][current_key]
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging,
        )
            (feasible, new_set) = ngroute_check_create_fset(
                neighborhoods, current_set, next_node,
            )
            !feasible && continue

            (feasible, new_subpath) = compute_new_subpath_nocharge(
                current_subpath, graph,
                current_node, next_node, 
                modified_costs, T_heuristic,
            )
            !feasible && continue

            new_λ_labels = compute_new_lambda_labels!(
                new_subpath, current_λ_labels, λvals, λcust,
            )

            if !(new_set in keys(base_labels[(starting_node, next_node)]))
                base_labels[(starting_node, next_node)][new_set] = SortedDict{
                    Tuple{Float64, Int, BitVector},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (
                new_subpath.cost,
                new_subpath.time_taken,
                new_λ_labels,
            )
            added = add_subpath_longlabel_to_collection_ngroute_lambda!(
                base_labels[(starting_node, next_node)][new_set],
                new_key, new_subpath, λvals,
                ;
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set, starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for depot in graph.N_depots
        for set in keys(base_labels[(depot, depot)])
            for path in values(base_labels[(depot, depot)][set])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute_alt_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{NTuple{3, Int}, Float64},
    T_heuristic::Int,
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    λvals, λcust = prepare_lambda(λ, graph.n_nodes)

    base_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, BitVector, BitVector}, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in union(graph.N_depots, graph.N_customers)
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, BitVector, Int, Int}}()
    for node in graph.N_depots
        λ_labels = falses(length(λ))
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        key = (0.0, 0, λ_labels, node_labels,)
        base_labels[(node, node)][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node, # starting_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(base_labels[(starting_node, current_node)]))
            continue
        end
        current_λ_labels = state[3]
        current_set = state[4][1:graph.n_nodes]
        current_subpath = base_labels[(starting_node, current_node)][current_key]
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging,
        )
            (feasible, new_set) = ngroute_check_create_fset(
                neighborhoods, current_set, next_node,
            )
            !feasible && continue

            (feasible, new_subpath) = compute_new_subpath_nocharge(
                current_subpath, graph,
                current_node, next_node, 
                modified_costs, T_heuristic,
            )
            !feasible && continue

            new_λ_labels = compute_new_lambda_labels!(
                new_subpath, current_λ_labels, λvals, λcust,
            )

            new_key = (
                new_subpath.cost,
                new_subpath.time_taken,
                new_λ_labels,
                new_set,
            )
            added = add_subpath_longlabel_to_collection_ngroute_alt_lambda!(
                base_labels[(starting_node, next_node)],
                new_key, new_subpath, λvals,
                ;
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for depot in graph.N_depots
        for path in values(base_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end

unwrap_base_labels(b::BaseSubpathLabel) = BaseSubpathLabel[b]

function unwrap_base_labels(d::AbstractDict)
    u = BaseSubpathLabel[]
    for v in values(d)
        append!(u, unwrap_base_labels(v))
    end
    return u
end

function get_negative_base_labels_from_base_labels(
    base_labels::Dict{
        NTuple{2, Int}, 
        T,
    },
) where {T <: AbstractDict}
    return BaseSubpathLabel[
        base_label
        for base_label in unwrap_base_labels(base_labels)
            if base_label.cost < -1e-6
    ]
end

function subproblem_iteration_nocharge(
    data,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{
        T,
        Float64,
    },
    T_heuristic::Int,
    ;
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    elementary::Bool = true,
    time_limit::Float64 = Inf,
) where {T}
    start_time = time()
    if ngroute && !ngroute_alt
        if length(λ) == 0
            base_labels_result = @timed find_nondominated_paths_nocharge_ngroute(
                data, graph, neighborhoods, 
                κ, μ, ν, T_heuristic,
                ;
                time_limit = time_limit - (time() - start_time),
            )
        else
            base_labels_result = @timed find_nondominated_paths_nocharge_ngroute_lambda(
                data, graph, neighborhoods, 
                κ, μ, ν, λ, T_heuristic,
                ;
                time_limit = time_limit - (time() - start_time),
            )
        end
    elseif ngroute && ngroute_alt
        if length(λ) == 0
            base_labels_result = @timed find_nondominated_paths_nocharge_ngroute_alt(
                data, graph, neighborhoods, 
                κ, μ, ν, T_heuristic,
                ;
                time_limit = time_limit - (time() - start_time),
            )
        else
            base_labels_result = @timed find_nondominated_paths_nocharge_ngroute_alt_lambda(
                data, graph, neighborhoods, 
                κ, μ, ν, λ, T_heuristic,
                ;
                time_limit = time_limit - (time() - start_time),
            )
        end
    else
        base_labels_result = @timed find_nondominated_paths_nocharge(
            data, graph, 
            κ, μ, ν, T_heuristic,
            ;
            elementary = elementary,
            time_limit = time_limit - (time() - start_time),
        )
    end

    negative_base_labels = get_negative_base_labels_from_base_labels(base_labels_result.value)
    negative_base_labels_count = length(negative_base_labels)
    return (negative_base_labels, negative_base_labels_count, base_labels_result.time)
end

function get_paths_from_negative_base_labels(
    graph::EVRPGraph,
    base_labels::Vector{BaseSubpathLabel},
)
    generated_paths = Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    }()
    for base_label in base_labels
        p = Path(
            subpaths = [
                Subpath(
                    n_customers = graph.n_customers,
                    starting_node = base_label.nodes[1],
                    starting_time = 0,
                    starting_charge = graph.B,
                    current_node = base_label.nodes[end],
                    arcs = collect(zip(base_label.nodes[1:end-1], base_label.nodes[2:end])),
                    current_time = base_label.time_taken,
                    current_charge = graph.B,
                    served = base_label.served,
                )
            ],
            charging_arcs = ChargingArc[],
            served = base_label.served,
        )
        add_path_to_generated_paths!(generated_paths, p)
    end
    return generated_paths
end

function path_formulation_column_generation_nocharge!(
    model::Model,
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    },
    SR3_constraints::Dict{
        T,
        ConstraintRef,
    },
    data::EVRPData, 
    graph::EVRPGraph,
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    path_costs::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Int},
    },
    path_service::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        Vector{Int},
    },
    printlist::Vector{String},
    ;
    time_heuristic_slack::Float64 = 0.9,
    elementary::Bool = true,
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    verbose::Bool = true,
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
) where {T}

    start_time = time()
    counter = 0
    converged = false
    local CGLP_results = Dict{String, Any}()

    T_heuristic = Int(round(time_heuristic_slack * (graph.T + graph.B) * graph.μ / (1 + graph.μ)))

    CG_params = Dict{String, Any}()
    CG_params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
    CG_params["objective"] = Float64[]
    CG_params["κ"] = Dict{Int, Float64}[]
    CG_params["μ"] = Dict{Int, Float64}[]
    CG_params["ν"] = Vector{Float64}[]
    CG_params["λ"] = Dict{keytype(SR3_constraints), Float64}[]
    CG_params["lp_relaxation_solution_time_taken"] = Float64[]
    CG_params["sp_base_time_taken"] = Float64[]
    CG_params["sp_full_time_taken"] = Float64[]
    CG_params["sp_total_time_taken"] = Float64[]
    CG_params["lp_relaxation_constraint_time_taken"] = Float64[]
    CG_params["number_of_new_paths"] = Int[]
    CG_params["converged"] = false

    while (
        !converged
        && time_limit ≥ (time() - start_time)
        && max_iters > counter
    )
        counter += 1
        mp_solution_start_time = time()
        @suppress optimize!(model)
        mp_solution_end_time = time()
        CGLP_results = Dict(
            "objective" => objective_value(model),
            "z" => Dict(
                (key, p) => value.(z[(key, p)])
                for (key, p) in keys(z)
            ),
            "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
            "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
            "ν" => dual.(model[:ν]).data,
            "λ" => Dict{keytype(SR3_constraints), Float64}(
                S => dual(SR3_constraints[S])
                for S in keys(SR3_constraints)
            ),
        )
        push!(CG_params["objective"], CGLP_results["objective"])
        push!(CG_params["κ"], CGLP_results["κ"])
        push!(CG_params["μ"], CGLP_results["μ"])
        push!(CG_params["ν"], CGLP_results["ν"])
        push!(CG_params["λ"], CGLP_results["λ"])
        push!(CG_params["lp_relaxation_solution_time_taken"], round(mp_solution_end_time - mp_solution_start_time, digits = 3))

        local negative_base_labels
        local base_labels_time
        try 
            (negative_base_labels, _, base_labels_time) = subproblem_iteration_nocharge(
                data, graph, 
                CGLP_results["κ"], 
                CGLP_results["μ"], 
                CGLP_results["ν"], 
                CGLP_results["λ"], 
                T_heuristic,
                ;
                neighborhoods = neighborhoods,
                ngroute = ngroute,
                ngroute_alt = ngroute_alt,
                elementary = elementary,
                time_limit = (time_limit - (time() - start_time)),
            )
        catch e
            if isa(e, TimeLimitException)
                break
            else
                throw(e)
            end
        end
        (generated_paths) = get_paths_from_negative_base_labels(
            graph, negative_base_labels,
        )
        push!(
            CG_params["sp_base_time_taken"],
            round(base_labels_time, digits=3)
        )
        push!(
            CG_params["sp_full_time_taken"],
            0.0
        )
        push!(
            CG_params["sp_total_time_taken"],
            round(base_labels_time, digits=3)
        )
        
        if length(generated_paths) == 0
            push!(CG_params["number_of_new_paths"], 0)
            converged = true
        else
            push!(
                CG_params["number_of_new_paths"],
                sum(length(v) for v in values(generated_paths))
            )
        end

        mp_constraint_time = add_paths_to_path_model!(
            model,
            z,
            some_paths, 
            path_costs,
            path_service,
            generated_paths,
            Dict{Tuple{Int, 3}, ConstraintRef}(),
            data, graph,
        )

        push!(
            CG_params["number_of_paths"], 
            sum(length(v) for v in values(some_paths))
        )
        push!(
            CG_params["lp_relaxation_constraint_time_taken"],
            mp_constraint_time,
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %9.3f | %9.3f | %9.3f | %9.3f | %10d \n", 
                counter,
                CG_params["objective"][counter],
                CG_params["number_of_paths"][counter],
                CG_params["lp_relaxation_solution_time_taken"][counter],
                CG_params["sp_base_time_taken"][counter],
                CG_params["sp_full_time_taken"][counter],                
                CG_params["lp_relaxation_constraint_time_taken"][counter],
                CG_params["number_of_new_paths"][counter],
            ),
            verbose,
        )
    end

    CG_params["converged"] = converged
    CG_params["counter"] = counter
    end_time = time() 
    time_taken = round(end_time - start_time, digits = 3)
    CG_params["time_taken"] = time_taken
    CG_params["time_limit_reached"] = (time_taken > time_limit)
    CG_params["lp_relaxation_time_taken"] = sum.(zip(CG_params["lp_relaxation_constraint_time_taken"], CG_params["lp_relaxation_solution_time_taken"]))
    CG_params["lp_relaxation_time_taken_total"] = sum(CG_params["lp_relaxation_time_taken"])
    CG_params["sp_base_time_taken_total"] = sum(CG_params["sp_base_time_taken"])
    CG_params["sp_full_time_taken_total"] = sum(CG_params["sp_full_time_taken"])
    CG_params["sp_time_taken_total"] = CG_params["sp_base_time_taken_total"] + CG_params["sp_full_time_taken_total"]
    CG_params["lp_relaxation_time_taken_mean"] = CG_params["lp_relaxation_time_taken_total"] / length(CG_params["lp_relaxation_time_taken"])
    CG_params["sp_base_time_taken_mean"] = CG_params["sp_base_time_taken_total"] / length(CG_params["sp_base_time_taken"])
    CG_params["sp_full_time_taken_mean"] = CG_params["sp_full_time_taken_total"] / length(CG_params["sp_full_time_taken"])
    CG_params["sp_time_taken_mean"] = CG_params["sp_base_time_taken_mean"] + CG_params["sp_full_time_taken_mean"]

    add_message!(
        printlist, 
        @sprintf(
            "Total         |            | %10d | %9.3f | %9.3f | %9.3f | %9.3f | \n", 
            CG_params["number_of_paths"][end],
            sum(CG_params["lp_relaxation_solution_time_taken"]),
            sum(CG_params["sp_base_time_taken"]),
            sum(CG_params["sp_full_time_taken"]),
            sum(CG_params["lp_relaxation_constraint_time_taken"]),
        ),
        verbose,
    )

    return CGLP_results, CG_params
end

function get_postcharge_shortest_pure_path_label(
    data::EVRPData, 
    graph::EVRPGraph,
    nodelist::Vector{Int},
    ;
    time_windows::Bool = false,
)
    modified_costs = compute_arc_modified_costs(graph, data, zeros(Float64, graph.n_customers))
    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes)
        β = repeat([graph.T], graph.n_nodes)
    end

    pure_path_labels = Dict(
        (nodelist[1:j]...,) => Dict(
            current_node => SortedDict{
                Tuple{Float64, Vararg{Int, 3}}, 
                # Tuple{Vararg{Int, 3}},
                PurePathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes
        )
        for j in eachindex(nodelist)
    )
    key = (0.0, 0, -graph.B, graph.B)
    # key = (0, -graph.B, graph.B)
    starting_node = nodelist[1]
    nodeseq = (starting_node,)
    pure_path_labels[nodeseq][starting_node][key] = PurePathLabel(
        0.0, 
        [starting_node],
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
    # unexplored_states = SortedSet{Tuple{Float64, Vararg{Int}}}()
    unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
    push!(unexplored_states, (key..., starting_node, nodeseq...))

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        current_node = state[5]
        current_nodeseq = state[6:end]
        current_key = state[1:4]
        # current_node = state[4]
        # current_nodeseq = state[5:end]
        # current_key = state[1:3]
        if !(current_key in keys(pure_path_labels[current_nodeseq][current_node]))
            continue
        end
        current_path = pure_path_labels[current_nodeseq][current_node][current_key]
        # println("current_path: $(current_path.nodes)")
        eventual_next_node = nodelist[length(current_nodeseq) + 1]
        for next_node in setdiff(
            vcat(eventual_next_node, graph.N_charging), 
            current_node,
        )
            if !(next_node in outneighbors(graph.G, current_node))
                continue
            end
            
            (feasible, new_path) = compute_new_pure_path(
                current_path,
                current_node,
                next_node, 
                data, graph,
                α, β, modified_costs,
            )
            !feasible && continue

            # add new_path to collection
            if next_node in union(graph.N_customers, graph.N_depots)
                new_nodeseq = (current_nodeseq..., next_node)
            else
                new_nodeseq = current_nodeseq
            end
            new_key = (
                # new_path.cost,
                0.0,
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
            )
            # println("$next_node, $new_nodeseq, $new_key")
            added = add_pure_path_label_to_collection!(
                pure_path_labels[new_nodeseq][next_node], 
                new_key, new_path, 
                ;
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., next_node, new_nodeseq...,)
                push!(unexplored_states, new_state)
            end
        end
    end
    if length(pure_path_labels[(nodelist...)][nodelist[end]]) == 0
        return (false, nothing)
    else
        return (true, first(values(pure_path_labels[(nodelist...)][nodelist[end]])))
    end
end

function delete_paths_by_time_length_from_model!(
    model::Model,
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
            Int,
        }, 
        VariableRef,
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    path_costs::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Int},
    },
    path_service::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        Vector{Int},
    },
    T_cutoff::Int,
    graph::EVRPGraph,
)
    for state_pair in keys(some_paths)
        delete_inds = Int[]
        for (count, path) in enumerate(some_paths[state_pair])
            # check if path needs to be removed
            removed = path.subpaths[end].current_time > T_cutoff
            if removed
                push!(delete_inds, count)
            end
        end
        # delete variables from model
        orig_inds = 1:length(some_paths[state_pair])
        for count in delete_inds
            delete(model, z[(state_pair, count)])
            pop!(z, (state_pair, count))
        end
        # rename variable bindings in dictionary `z`
        for (i, ind) in enumerate(sort(setdiff(orig_inds, delete_inds)))
            z[(state_pair, i)] = pop!(z, (state_pair, ind))
        end
        # delete paths, their costs, and service info via `delete_inds`
        deleteat!(some_paths[state_pair], delete_inds)
        deleteat!(path_costs[state_pair], delete_inds)
        for node in graph.N_customers
            deleteat!(path_service[(state_pair, node)], delete_inds)
        end
    end
end

function path_formulation_decomposition_heuristic(
    data::EVRPData,
    graph::EVRPGraph,
    ;
    Env = nothing,
    time_heuristic_slack::Float64 = 0.9, 
    elementary::Bool = true,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size::String = "small",
    ngroute_neighborhood_charging_size::String = "small",
    verbose::Bool = true,
    use_adaptive_ngroute::Bool = true,
    use_SR3_cuts::Bool = true,
    use_lmSR3_cuts::Bool = true,
    max_SR3_cuts::Int = 10, 
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
)
    # if use_SR3_cuts || use_lmSR3_cuts
    #     error("Not yet implemented!")
    # end

    start_time = time()

    if ngroute && isnothing(neighborhoods)
        neighborhoods = compute_ngroute_neighborhoods(
            graph,
            ngroute_neighborhood_size; 
            depots_size = ngroute_neighborhood_depots_size,
            charging_size = ngroute_neighborhood_charging_size,
        )
    end

    some_paths = generate_artificial_paths(data, graph)
    path_costs = compute_path_costs(
        data, graph, 
        some_paths,
    )
    path_service = compute_path_service(
        graph,
        some_paths,
    )

    printlist = String[]

    start_printlist = String[]
    push!(start_printlist, 
        @sprintf(
            """
            Starting column generation on the path formulation.
            # customers:                    %3d
            # depots:                       %3d
            # charging stations:            %3d
            # vehicles:                     %3d
            """,
            graph.n_customers,
            graph.n_depots,
            graph.n_charging,
            data.n_vehicles,
        )
    )
    if ngroute
        push!(start_printlist, 
            @sprintf(
                """
                ngroute:                        %s
                ngroute_alt:                    %s
                ngroute neighborhood size:
                    customers                   %3d
                    depots                      %s
                    charging                    %s
                """,
                ngroute,
                ngroute_alt,
                ngroute_neighborhood_size,
                ngroute_neighborhood_depots_size,
                ngroute_neighborhood_charging_size,
            )
        )
        push!(start_printlist,
            @sprintf(
                "use_adaptive_ngroute:           %s\n",
                use_adaptive_ngroute
            )
        )
        if use_SR3_cuts
            push!(start_printlist,
                @sprintf(
                    """
                    use_SR3_cuts:                   %s
                        use_lmSR3_cuts:             %s
                        max_SR3_cuts:               %3d
                    """,
                    use_SR3_cuts,
                    use_lmSR3_cuts,
                    max_SR3_cuts,
                )
            )
        else
            push!(start_printlist,
                @sprintf(
                    """
                    use_SR3_cuts:                   %s
                    """,
                    use_SR3_cuts,
                )
            )
        end
    else
        push!(start_printlist,
            @sprintf(
                """
                ngroute:                        %s
                elementary:                     %s
                """,
                ngroute,
                elementary,
            )
        )
    end

    for message in start_printlist
        add_message!(printlist, message, verbose)
    end

    add_message!(
        printlist,
        @sprintf("              |  Objective | #    paths | Time (LP) | Time (SP) | Time (SP) | Time (LP) | #    paths \n"),
        verbose,
    )
    add_message!(
        printlist,
        @sprintf("              |            |            |           |    (base) |    (full) |    (cons) |      (new) \n"),
        verbose,
    )

    model, z = path_formulation_build_model(
        data, graph, some_paths, path_costs, path_service,
        ; 
        Env = Env,
    )

    local CGLP_results
    local CG_params
    local CGIP_results
    local converged = false

    if use_lmSR3_cuts
        SR3_constraints = Dict{
            Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
            ConstraintRef,
        }()
        SR3_list = Tuple{Float64, NTuple{3, Int}, Tuple{Vararg{Int}}}[]
    else
        SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
        SR3_list = Tuple{Float64, NTuple{3, Int}}[]
    end
    CGLP_all_results = []
    CGIP_all_results = []
    all_params = []
    CG_all_params = []
    if ngroute
        CG_all_neighborhoods = BitMatrix[]
    else
        CG_all_neighborhoods = nothing
    end

    local results_paths_withcharge = Tuple{Float64, Path}[]
    # outer loop for heuristic
    while true
        # loop for column generation w/ cuts
        while true
            continue_flag = false
            iteration_params = Dict{String, Any}()

            CGLP_results, CG_params = path_formulation_column_generation_nocharge!(
                model, z, SR3_constraints, 
                data, graph,
                some_paths, path_costs, path_service,
                printlist,
                ;
                time_heuristic_slack = time_heuristic_slack,
                elementary = elementary,
                neighborhoods = neighborhoods,
                ngroute = ngroute,
                ngroute_alt = ngroute_alt,
                verbose = verbose,
                time_limit = time_limit - (time() - start_time),
                max_iters = max_iters,
            )
            CGIP_results = path_formulation_solve_integer_model!(
                model,
                z,
            )
            push!(CGLP_all_results, CGLP_results)
            push!(CGIP_all_results, CGIP_results)
            push!(CG_all_params, CG_params)
            if ngroute
                push!(CG_all_neighborhoods, copy(neighborhoods))
            end

            # Termination criteria
            CG_params["CGLP_objective"] = CGLP_results["objective"]
            CG_params["CGIP_objective"] = CGIP_results["objective"]
            CG_params["LP_IP_gap"] = 1.0 - CGLP_results["objective"] / CGIP_results["objective"]
            for message in [
                @sprintf("\n"),
                @sprintf("Time taken (s):       %9.3f s\n", CG_params["time_taken"]),
                @sprintf("(CGLP) Objective:         %.4e\n", CGLP_results["objective"]),
                @sprintf("(CGIP) Objective:         %.4e\n", CGIP_results["objective"]),
                @sprintf("%% gap:                %9.3f %%\n", CG_params["LP_IP_gap"] * 100.0),
            ]
                add_message!(printlist, message, verbose)
            end

            if CG_params["converged"]
                CGLP_results["paths"] = collect_path_solution_support(
                    CGLP_results, some_paths, data, graph
                )
                CGIP_results["paths"] = collect_path_solution_support(
                    CGIP_results, some_paths, data, graph
                )
            end

            iteration_params["CGLP_objective"] = CG_params["CGLP_objective"]
            iteration_params["CGIP_objective"] = CG_params["CGIP_objective"]
            iteration_params["CG_LP_IP_gap"] = CG_params["LP_IP_gap"]
            iteration_params["CG_time_taken"] = CG_params["time_taken"]
            iteration_params["CG_sp_time_taken_mean"] = CG_params["sp_time_taken_mean"]
            iteration_params["method"] = "none"
            iteration_params["ngroute_neighborhood_size"] = ngroute ? (graph.n_customers + graph.n_charging) * mean(neighborhoods[graph.N_customers, vcat(graph.N_customers, graph.N_charging)]) : 0
            iteration_params["cycles_lookup_length"] = 0
            iteration_params["implemented_SR3_cuts_count"] = 0
            iteration_params["converged"] = false
            iteration_params["time_limit_reached"] = false
            # check if converged
            if CGIP_results["objective"] ≈ CGLP_results["objective"]
                converged = true
                iteration_params["converged"] = true
            end

            if CG_params["converged"] && !continue_flag && !converged && ngroute && use_adaptive_ngroute
                # see if any path in solution was non-elementary
                cycles_lookup = detect_cycles_in_path_solution([p for (val, p) in CGLP_results["paths"]], graph)
                if length(cycles_lookup) > 0 
                    delete_paths_with_found_cycles_from_model!(model, z, some_paths, path_costs, path_service, cycles_lookup, graph)
                    modify_neighborhoods_with_found_cycles!(neighborhoods, cycles_lookup)
                    continue_flag = true
                    iteration_params["method"] = "use_adaptive_ngroute"
                    iteration_params["cycles_lookup_length"] = length(cycles_lookup)
                    add_message!(printlist, "Expanded ng-route neighborhoods by $(length(cycles_lookup))\n", verbose)
                    add_message!(printlist, "\n", verbose)
                end
            end

            if CG_params["converged"] && !continue_flag && !converged && ngroute && use_SR3_cuts
                # if no path in solution was non-elementary, 
                # generate violated WSR3 inequalities
                generated_WSR3_list = enumerate_violated_path_WSR3_inequalities(
                    CGLP_results["paths"], 
                    graph,
                )
                generated_SR3_list = [(t[1], t[2:4]) for t in generated_WSR3_list]
                if length(generated_SR3_list) == 0
                    # backup: generate violated SR3 inequalities
                    generated_SR3_list = enumerate_violated_path_SR3_inequalities(
                        CGLP_results["paths"],
                        graph,
                    )
                    add_message!(printlist, "Found SR3 cuts:      $(length(generated_SR3_list))\n", verbose)
                    # sample cuts if too many
                    if length(generated_SR3_list) ≤ max_SR3_cuts
                        implemented_SR3_list = generated_SR3_list
                    else
                        implemented_SR3_list = sample(
                            generated_SR3_list, 
                            Weights([val for (val, _) in generated_SR3_list]),
                            max_SR3_cuts, 
                            replace = false,
                        )
                        add_message!(printlist, "Sampled SR3 cuts:    $(length(implemented_SR3_list))\n", verbose)
                    end
                else
                    add_message!(printlist, "Found WSR3 cuts:     $(length(generated_SR3_list))\n", verbose)
                    implemented_SR3_list = generated_SR3_list
                end
                if length(implemented_SR3_list) != 0
                    if use_lmSR3_cuts
                        implemented_SR3_list = [
                            (val, S, compute_memory_set_of_lmSRnk_inequality(CGLP_results["paths"], S, 2))
                            for (val, S) in implemented_SR3_list
                        ]
                        append!(SR3_list, implemented_SR3_list)
                        add_lmSR3_constraints_to_path_model!(
                            model, z, some_paths, 
                            SR3_constraints, implemented_SR3_list,
                        )
                        iteration_params["method"] = "use_lmSR3_cuts"
                        iteration_params["implemented_SR3_cuts_count"] = length(implemented_SR3_list)
                        continue_flag = true
                        add_message!(printlist, "Imposed lm-SR3 cuts: $(length(implemented_SR3_list))\n", verbose)
                        for (val, S, M) in implemented_SR3_list
                            add_message!(printlist, "$S, $M: $val\n", verbose)
                        end
                        add_message!(printlist, "\n", verbose)
                    else
                        append!(SR3_list, implemented_SR3_list)
                        # Add violated inequalities to master problem
                        add_SR3_constraints_to_path_model!(
                            model, z, some_paths, 
                            SR3_constraints, implemented_SR3_list, 
                        )
                        iteration_params["method"] = "use_SR3_cuts"
                        iteration_params["implemented_SR3_cuts_count"] = length(implemented_SR3_list)
                        continue_flag = true
                        add_message!(printlist, "Imposed SR3 cuts:    $(length(implemented_SR3_list))\n", verbose)
                        for (val, S) in implemented_SR3_list
                            add_message!(printlist, "$S: $val\n", verbose)
                        end
                        add_message!(printlist, "\n", verbose)
                    end
                end
            end

            push!(all_params, iteration_params)

            if !(time_limit > time() - start_time)
                iteration_params["time_limit_reached"] = true
                break
            end
            
            if !continue_flag
                break
            end
        end
        
        feasible = true

        results_paths_withcharge = Tuple{Float64, Path}[]
        for (val, p) in CGIP_all_results[end]["paths"]
            if p.subpaths[1].artificial
                p_feasible = false
                feasible = false
                break
            end
            nodelist = vcat(p.arcs[1][1], [x[2] for x in p.arcs])
            if length(nodelist) == 2 && p.subpaths[1].current_time == 0
                push!(results_paths_withcharge, (val, p))
                continue
            end
            (p_feasible, pure_path_label) = get_postcharge_shortest_pure_path_label(
                data, graph, nodelist,
                ;
                time_windows = false,
            )
            feasible = feasible && p_feasible
            if !p_feasible
                break
            end
            p_ = convert_pure_path_label_to_path(pure_path_label, graph)
            push!(results_paths_withcharge, (val, p_))
        end

        if feasible
            break
        else
            time_heuristic_slack -= 0.05
            if time_heuristic_slack ≤ 0.0
                heuristic_results = Dict(
                    "objective" => Inf,
                    "paths" => Tuple{Float64, Path}[],
                )
                return (
                    CGIP_all_results[end],
                    heuristic_results,
                    time_heuristic_slack,
                )
            end
        end

        # filter all paths that are shorter than time_heuristic
        T_heuristic = Int(round(time_heuristic_slack * (graph.T + graph.B) * graph.μ / (1 + graph.μ)))
        delete_paths_by_time_length_from_model!(
            model, z, 
            some_paths, path_costs, path_service,
            T_heuristic, graph,
        )
    end

    heuristic_results = Dict(
        "objective" => sum(
            val * compute_path_cost(data, graph, p)
            for (val, p) in results_paths_withcharge
        ),
        "paths" => results_paths_withcharge,
    )

    return (
        CGIP_all_results[end],
        heuristic_results,
        time_heuristic_slack,
    )

end