include("utils.jl")
using DataStructures
using Printf

mutable struct BaseSubpathLabel
    time_taken::Int
    charge_taken::Int
    cost::Float64
    nodes::Vector{Int}
    served::Vector{Int}
end

Base.copy(s::BaseSubpathLabel) = BaseSubpathLabel(
    s.time_taken,
    s.charge_taken,
    s.cost, 
    copy(s.nodes),
    copy(s.served),
)

mutable struct PathLabel
    cost::Float64
    subpath_labels::Vector{BaseSubpathLabel}
    charging_actions::Vector{Int}
    served::Vector{Int}
end

Base.copy(p::PathLabel) = PathLabel(
    p.cost,
    [s for s in p.subpath_labels],
    [t for t in p.charging_actions],
    copy(p.served),
)

function add_subpath_longlabel_to_collection!(
    collection::SortedDict{
        NTuple{N, Int},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::NTuple{N, Int},
    v1::BaseSubpathLabel,
    ;
    verbose::Bool = false,
) where {N}
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if v2.cost ≤ v1.cost
            if all(k2 .≤ k1)
                added = false
                if verbose
                    println("$(k1), $(v1.cost) dominated by $(k2), $(v2.cost)")
                end
                break
            end
        end
        # check if v1 dominates v2
        if v1.cost ≤ v2.cost
            if all(k1 .≤ k2)
                if verbose
                    println("$(k1), $(v1.cost) dominates $(k2), $(v2.cost)")
                end
                pop!(collection, k2)
            end
        end
    end
    if added
        if verbose
            println("$(k1), $(v1.cost) added!")
        end
        insert!(collection, k1, v1)
    end
    return added
end

function add_path_label_to_collection!(
    collection::SortedDict{
        NTuple{N, Int},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::NTuple{N, Int},
    v1::PathLabel,
    ;
    verbose::Bool = false,
) where {N}
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if v2.cost ≤ v1.cost
            if all(k2 .≤ k1)
                added = false
                if verbose
                    println("$(k1), $(v1.cost) dominated by $(k2), $(v2.cost)")
                end
                break
            end
        end
        # check if v1 dominates v2
        if v1.cost ≤ v2.cost
            if all(k1 .≤ k2)
                if verbose
                    println("$(k1), $(v1.cost) dominates $(k2), $(v2.cost)")
                end
                pop!(collection, k2)
            end
        end
    end
    if added
        if verbose
            println("$(k1), $(v1.cost) added!")
        end
        insert!(collection, k1, v1)
    end
    return added
end


function ngroute_extend_partial_path_check(
    neighborhoods::BitMatrix,
    set::Tuple{Vararg{Int}},
    s::BaseSubpathLabel,
)
    new_set = collect(set)
    for next_node in s.nodes[2:end]
        if next_node in new_set
            return (nothing, false)
        end
        new_set = [
            node for node in new_set
                if neighborhoods[next_node, node]
        ]
        push!(new_set, next_node)
        # println("$next_node, $new_set")
    end
    return (Tuple(sort(unique(new_set))), true)
end

function ngroute_extend_partial_path_check_alt(
    neighborhoods::BitMatrix,
    set::Vector{Int},
    s::BaseSubpathLabel,
    graph::EVRPGraph,
)
    new_set = copy(set)
    for next_node in s.nodes[2:end]
        if new_set[next_node] == 1
            return (nothing, false)
        end
        for node in graph.N_nodes_extra
            if new_set[node] == 1 && !(neighborhoods[next_node, node])
                new_set[node] = 0
            end
        end
        new_set[next_node] = 1
        # println("$next_node, $new_set")
    end
    return (new_set, true)
end

function compute_next_subpath(
    current_subpath::BaseSubpathLabel,
    graph::EVRPGraph,
    current_node::Int,
    next_node::Int,
    modified_costs::Matrix{Float64},
)

    # time and charge feasibility
    # if current_subpath.time_taken + graph.t[current_node, next_node] + graph.min_t[next_node] > graph.T
    #     return (false, nothing)
    # end 

    if current_subpath.charge_taken + graph.q[current_node, next_node] + graph.min_q[next_node] > graph.B
        return (false, nothing)
    end

    new_subpath = copy(current_subpath)
    new_subpath.time_taken += graph.t[current_node, next_node]
    new_subpath.charge_taken += graph.q[current_node, next_node]
    new_subpath.cost += modified_costs[current_node, next_node]
    push!(new_subpath.nodes, next_node)
    if next_node in graph.N_customers
        new_subpath.served[next_node] += 1
    end
    return (true, new_subpath)
end

function generate_base_labels_nonsingleservice(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{1, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{3, Int}}()
    for node in graph.N_depots_charging_extra
        base_labels[node][node][(0,)] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                0, 
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
        current_key = state[1:end-2]
        if !(current_key in keys(base_labels[starting_node][current_node]))
            continue
        end
        current_subpath = base_labels[starting_node][current_node][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            # time and charge feasibility
            # if current_subpath.time_taken + graph.t[current_node, next_node] + graph.min_t[next_node] > graph.T
            #     continue
            # end 
            if current_subpath.charge_taken + graph.q[current_node, next_node] + graph.min_q[next_node] > graph.B
                continue
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken += graph.t[current_node, next_node]
            new_subpath.charge_taken += graph.q[current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in graph.N_customers
                new_subpath.served[next_node] += 1
            end

            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node], 
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], (0,))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end    
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end

function generate_base_labels_nonsingleservice_christofides(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int},
                SortedDict{
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{5, Int}}()
    for node in graph.N_depots_charging_extra
        base_labels[node][node][(node, node)] = SortedDict{
            NTuple{1, Int}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            (0,) => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            ),
        )
        push!(
            unexplored_states, 
            (
                0, 
                node, # starting_node 
                node, # first_node
                node, # prev_node
                node, # current_node
            )            
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        first_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-4]
        if !(current_key in keys(base_labels[starting_node][current_node][(first_node, prev_node)]))
            continue
        end
        current_subpath = base_labels[starting_node][current_node][(first_node, prev_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            # Preventing customer 2-cycles (Christofides)
            if next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == prev_node == next_node
                    continue
                end
            end
            # time and charge feasibility
            # if current_subpath.time_taken + graph.t[current_node, next_node] + graph.min_t[next_node] > graph.T
            #     continue
            # end 
            if current_subpath.charge_taken + graph.q[current_node, next_node] + graph.min_q[next_node] > graph.B
                continue
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken += graph.t[current_node, next_node]
            new_subpath.charge_taken += graph.q[current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in graph.N_customers
                new_subpath.served[next_node] += 1
            end

            if first_node == starting_node
                first_node = current_node
            end
            if !((first_node, current_node) in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][(first_node, current_node)] = SortedDict{
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][(first_node, current_node)], 
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, first_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], (node, node))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][(first_node, prev_node)])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][(first_node, prev_node)])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels

end

function generate_base_labels_singleservice(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    check_customers::Bool = false,
    christofides::Bool = false,
    time_limit::Float64 = Inf,
) # UNTESTED
    function add_subpath_label_to_collection!(
        collection::SortedDict{
            Int, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
        k1::Int,
        v1::BaseSubpathLabel,
        ;
        last::Int = 0,
    )
        last_assigned = false
        for (k2, v2) in collection
            if k2 < last
                continue
            end
            if k2 < k1
                if v2.cost ≤ v1.cost
                    return (false, k2)
                end
            else
                if !last_assigned
                    last_assigned = true
                    last = k2
                end
                if k1 < k2
                    if v1.cost ≤ v2.cost
                        pop!(collection, k2)
                    end
                end
            end
        end
        insert!(collection, k1, v1)
        last = k1
        return (true, k1)
    end

    function direct_sum_of_collections(
        labels1::SortedDict{
            Int, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
        labels2::SortedDict{
            Int, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
        ;
    )
        keys1 = collect(keys(labels1))
        keys2 = collect(keys(labels2))
    
        new = []
        for (t, cost, i, j) in sort([
            (k1 + k2, s1.cost + s2.cost, i, j)
            for (i, (k1, s1)) in enumerate(pairs(labels1)),
                (j, (k2, s2)) in enumerate(pairs(labels2))
                if s1.charge_taken + s2.charge_taken ≤ graph.B
                    && all(s1.served .+ s2.served .≤ 1)
                    # TODO: if I relax this, think of how to impose Christofides 2-cycle condition
        ])
            if length(new) == 0
                push!(new, (t, cost, i, j))
            elseif new[end][1] == t
                if new[end][2] > cost
                    new = new[1:end-1]
                    push!(new, (t, cost, i, j))
                end
            else
                if new[end][2] > cost
                    push!(new, (t, cost, i, j))
                end
            end
        end       
        new_labels = SortedDict{
            Int, 
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(
            Base.Order.ForwardOrdering(),
            t => BaseSubpathLabel(
                t,
                labels1[keys1[i]].charge_taken + labels2[keys2[j]].charge_taken,
                cost,
                vcat(labels1[keys1[i]].nodes, labels2[keys2[j]].nodes[2:end]),
                labels1[keys1[i]].served .+ labels2[keys2[j]].served
            )
            for (t, cost, i, j) in new
        )
        return new_labels
    end

    function merge_collections!(
        labels1::SortedDict{
            Int, 
            BaseSubpathLabel,            
            Base.Order.ForwardOrdering,
        },
        labels2::SortedDict{
            Int, 
            BaseSubpathLabel,            
            Base.Order.ForwardOrdering,
        },
    )
        last = 0
        while length(labels2) > 0
            (k, v) = first(labels2)
            pop!(labels2, k)
            (added, last) = add_subpath_label_to_collection!(
                labels1,
                k, v;
                last = last
            )
        end
        return 
    end

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    keylen = check_customers ? graph.n_customers + 1 : 1
    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{keylen, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_nodes_extra
    )

    for edge in edges(graph.G)
        starting_node = edge.src
        current_node = edge.dst
        time_taken = graph.t[starting_node, current_node]
        served = zeros(Int, graph.n_customers)
        if current_node in graph.N_customers
            served[current_node] = 1
        end
        if check_customers
            key = (time_taken, served...)
        else
            key = (time_taken,)
        end
        base_labels[starting_node][current_node][key] = BaseSubpathLabel(
            time_taken,
            graph.q[starting_node, current_node],
            modified_costs[starting_node, current_node],
            [starting_node, current_node],
            served,
        )
    end

    for new_node in graph.N_customers
        for starting_node in setdiff(graph.N_nodes_extra, new_node)
            if length(base_labels[starting_node][new_node]) == 0
                continue
            end
            for end_node in setdiff(graph.N_nodes_extra, new_node)
                if length(base_labels[new_node][end_node]) == 0
                    continue
                end
                if true
                    for (k1, s1) in pairs(base_labels[starting_node][new_node])
                        for (k2, s2) in pairs(base_labels[new_node][end_node])
                            if time_limit < time() - start_time
                                throw(TimeLimitException())
                            end
                            # Preventing customer 2-cycles (Christofides)
                            if christofides && s1.nodes[end-1] in graph.N_customers && s1.nodes[end-1] == s2.nodes[2]
                                continue
                            end
                            k = k1 .+ k2
                            if !all(s1.served .+ s2.served .≤ 1)
                                continue
                            end
                            if s1.charge_taken + s2.charge_taken > graph.B
                                continue
                            end
                            s = BaseSubpathLabel(
                                s1.time_taken + s2.time_taken,
                                s1.charge_taken + s2.charge_taken,
                                s1.cost + s2.cost,
                                vcat(s1.nodes, s2.nodes[2:end]),
                                s1.served .+ s2.served,
                            )
                            add_subpath_longlabel_to_collection!(
                                base_labels[starting_node][end_node],
                                k, s,
                            )
                        end
                    end
                else
                    merge_collections!(
                        base_labels[starting_node][end_node],
                        direct_sum_of_collections(
                            base_labels[starting_node][new_node], 
                            base_labels[new_node][end_node],
                        )
                    )
                end
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν)
    end

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{graph.n_nodes_extra, Int},
                SortedDict{
                    NTuple{1, Int},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 3, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        base_labels[node][node][(node_labels...)] = SortedDict{
            NTuple{1, Int}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            (0,) => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            )
        )
        push!(
            unexplored_states, 
            (
                0, 
                node_labels...,
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
        current_set = state[2:end-2]
        current_key = state[1:1]
        current_subpath = get(base_labels[starting_node][current_node][current_set], current_key, nothing)
        if isnothing(current_subpath)
            continue
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
                
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            if !(new_set in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][new_set] = SortedDict{ 
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][new_set],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        delete!(base_labels[node][node][(node_labels...)], (0,))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_christofides(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64},
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν)
    end

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int}, 
                Dict{
                    NTuple{graph.n_nodes_extra, Int},
                    SortedDict{
                        NTuple{1, Int},
                        BaseSubpathLabel,
                        Base.Order.ForwardOrdering,
                    },
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 5, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        base_labels[node][node][(node, node)] = Dict{
            NTuple{graph.n_nodes_extra, Int},
            SortedDict{
                NTuple{1, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }
        }(
            (node_labels...) => SortedDict{
                NTuple{1, Int}, 
                BaseSubpathLabel,
            }(
                Base.Order.ForwardOrdering(),
                (0,) => BaseSubpathLabel(
                    0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
                )
            )
        )
        push!(
            unexplored_states, 
            (
                0, 
                node_labels...,
                node, # starting_node
                node, # first_node
                node, # prev_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        first_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_set = state[2:end-4]
        current_key = state[1:1]
        current_subpath = get(base_labels[starting_node][current_node][(first_node, prev_node)][current_set], current_key, nothing)
        if isnothing(current_subpath)
            continue
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == prev_node == next_node
                    continue
                end
            end
                                
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            if first_node == starting_node
                first_node = current_node 
            end

            if !((first_node, current_node) in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][(first_node, current_node)] = Dict{
                    NTuple{graph.n_nodes_extra, Int},
                    SortedDict{
                        NTuple{1, Int},
                        BaseSubpathLabel,
                        Base.Order.ForwardOrdering,
                    },
                }()
            end

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            if !(new_set in keys(base_labels[starting_node][next_node][(first_node, current_node)]))
                base_labels[starting_node][next_node][(first_node, current_node)][new_set] = SortedDict{ 
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][(first_node, current_node)][new_set],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set..., starting_node, first_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], (node, node))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for set in keys(base_labels[starting_node][end_node][(first_node, prev_node)])
                    for v in values(base_labels[starting_node][end_node][(first_node, prev_node)][set])
                        v.cost = v.cost - κ[starting_node]
                    end
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for set in keys(base_labels[starting_node][end_node][(first_node, prev_node)])
                    for v in values(base_labels[starting_node][end_node][(first_node, prev_node)][set])
                        v.cost = v.cost - μ[end_node]
                    end
                end
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_sigma(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν; σ = σ)
    end
    σ_costs = compute_WSR3_sigma_costs(σ, graph)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => Dict{
                    NTuple{graph.n_nodes_extra, Int},
                    SortedDict{
                        NTuple{1, Int},
                        BaseSubpathLabel,
                        Base.Order.ForwardOrdering,
                    },
                }(Base.Order.ForwardOrdering())
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 4, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        base_labels[node][node][node][(node_labels...)] = SortedDict{
            NTuple{1, Int}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            (0,) => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            )
        )
        push!(
            unexplored_states, 
            (
                0, 
                node_labels...,
                node, # starting_node
                node, # prev_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_set = state[2:end-3]
        current_key = state[1:1]
        current_subpath = get(base_labels[starting_node][current_node][prev_node][current_set], current_key, nothing)
        if isnothing(current_subpath)
            continue
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
                                
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_subpath.cost += σ_costs[(prev_node, current_node, next_node)]

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            if !(new_set in keys(base_labels[starting_node][next_node][current_node]))
                base_labels[starting_node][next_node][current_node][new_set] = SortedDict{ 
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][current_node][new_set],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set..., starting_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], node)
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for prev_node in keys(base_labels[starting_node][end_node])
                for set in keys(base_labels[starting_node][end_node][prev_node])
                    for v in values(base_labels[starting_node][end_node][prev_node][set])
                        v.cost = v.cost - κ[starting_node]
                    end
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for prev_node in keys(base_labels[starting_node][end_node])
                for set in keys(base_labels[starting_node][end_node][prev_node])
                    for v in values(base_labels[starting_node][end_node][prev_node][set])
                        v.cost = v.cost - μ[end_node]
                    end
                end
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_sigma_christofides(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν; σ = σ)
    end
    σ_costs = compute_WSR3_sigma_costs(σ, graph)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int}, 
                Dict{
                    NTuple{graph.n_nodes_extra, Int},
                    SortedDict{
                        NTuple{1, Int},
                        BaseSubpathLabel,
                        Base.Order.ForwardOrdering,
                    },
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 5, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        base_labels[node][node][(node, node)] = Dict{
            NTuple{graph.n_nodes_extra, Int},
            SortedDict{
                NTuple{1, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }
        }(
            (node_labels...) => SortedDict{
                NTuple{1, Int}, 
                BaseSubpathLabel,
            }(
                Base.Order.ForwardOrdering(),
                (0,) => BaseSubpathLabel(
                    0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
                )
            )
        )
        push!(
            unexplored_states, 
            (
                0, 
                node_labels...,
                node, # starting_node
                node, # first_node
                node, # prev_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        first_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_set = state[2:end-4]
        current_key = state[1:1]
        current_subpath = get(base_labels[starting_node][current_node][(first_node, prev_node)][current_set], current_key, nothing)
        if isnothing(current_subpath)
            continue
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == prev_node == next_node
                    continue
                end
            end
                                
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_subpath.cost += σ_costs[(prev_node, current_node, next_node)]

            if first_node == starting_node
                first_node = current_node 
            end

            if !((first_node, current_node) in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][(first_node, current_node)] = Dict{
                    NTuple{graph.n_nodes_extra, Int},
                    SortedDict{
                        NTuple{1, Int},
                        BaseSubpathLabel,
                        Base.Order.ForwardOrdering,
                    },
                }()
            end

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            if !(new_set in keys(base_labels[starting_node][next_node][(first_node, current_node)]))
                base_labels[starting_node][next_node][(first_node, current_node)][new_set] = SortedDict{ 
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][(first_node, current_node)][new_set],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set..., starting_node, first_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], (node, node))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for set in keys(base_labels[starting_node][end_node][(first_node, prev_node)])
                    for v in values(base_labels[starting_node][end_node][(first_node, prev_node)][set])
                        v.cost = v.cost - κ[starting_node]
                    end
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for set in keys(base_labels[starting_node][end_node][(first_node, prev_node)])
                    for v in values(base_labels[starting_node][end_node][(first_node, prev_node)][set])
                        v.cost = v.cost - μ[end_node]
                    end
                end
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_alt(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()    
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν; σ = σ)
    end
    σ_costs = compute_WSR3_sigma_costs(σ, graph)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{graph.n_nodes_extra + 1, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 3, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        key = (0, node_labels...)
        base_labels[node][node][key] = BaseSubpathLabel(
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
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(base_labels[starting_node][current_node]))
            continue
        end
        current_set = state[2:end-2]
        current_subpath = base_labels[starting_node][current_node][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
                            
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_subpath.cost += σ_costs[(prev_node, current_node, next_node)]
            
            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            new_key = (new_subpath.time_taken, new_set...)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        key = (0, node_labels...)
        delete!(base_labels[node][node], key)
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_alt_christofides(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()    
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν)
    end

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int},
                SortedDict{
                    NTuple{graph.n_nodes_extra + 1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 5, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        key = (0, node_labels...)
        base_labels[node][node][(node, node)] = SortedDict{
            NTuple{graph.n_nodes_extra + 1, Int},
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            ),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node, # starting_node
                node, # first_node 
                node, # prev_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        first_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-4]
        if !(current_key in keys(base_labels[starting_node][current_node][(first_node, prev_node)]))
            continue
        end
        current_set = state[2:end-4]
        current_subpath = base_labels[starting_node][current_node][(first_node, prev_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == prev_node == next_node
                    continue
                end
            end
                            
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_subpath.cost += σ_costs[(prev_node, current_node, next_node)]

            if first_node == starting_node
                first_node = current_node
            end
            if !((first_node, current_node) in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][(first_node, current_node)] = SortedDict{
                    NTuple{graph.n_nodes_extra + 1, Int},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            new_key = (new_subpath.time_taken, new_set...)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][(first_node, current_node)],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, first_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], (node, node))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][(first_node, prev_node)])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][(first_node, prev_node)])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels
end


function generate_base_labels_ngroute_alt_sigma(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()    
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν; σ = σ)
    end
    σ_costs = compute_WSR3_sigma_costs(σ, graph)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => SortedDict{
                    NTuple{graph.n_nodes_extra + 1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 4, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        key = (0, node_labels...)
        base_labels[node][node][node][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node, # starting_node
                node, # prev_node
                node, # current_node
            )
        )
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
        if !(current_key in keys(base_labels[starting_node][current_node][prev_node]))
            continue
        end
        current_set = state[2:end-4]
        current_subpath = base_labels[starting_node][current_node][prev_node][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
                           
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_subpath.cost += σ_costs[(prev_node, current_node, next_node)]

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            new_key = (new_subpath.time_taken, new_set...)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][current_node],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], node)
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for prev_node in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][prev_node])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for prev_node in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][prev_node])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_alt_sigma_christofides(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    σ::Dict{Tuple{Vararg{Int}}, Float64},
    ;
    modified_costs::Union{Nothing, Matrix{Float64}} = nothing,
    time_limit::Float64 = Inf,
)

    start_time = time()    
    if isnothing(modified_costs)
        modified_costs = compute_arc_modified_costs(graph, data, ν; σ = σ)
    end
    σ_costs = compute_WSR3_sigma_costs(σ, graph)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{2, Int},
                SortedDict{
                    NTuple{graph.n_nodes_extra + 1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in graph.N_nodes_extra
        )
        for starting_node in graph.N_depots_charging_extra
    )
    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 5, Int}}()
    for node in graph.N_depots_charging_extra
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[node] = 1
        key = (0, node_labels...)
        base_labels[node][node][(node, node)] = SortedDict{
            NTuple{graph.n_nodes_extra + 1, Int},
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            ),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node, # starting_node
                node, # first_node 
                node, # prev_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-3]
        first_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-4]
        if !(current_key in keys(base_labels[starting_node][current_node][(first_node, prev_node)]))
            continue
        end
        current_set = state[2:end-4]
        current_subpath = base_labels[starting_node][current_node][(first_node, prev_node)][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if next_node in graph.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == prev_node == next_node
                    continue
                end
            end
                            
            (feasible, new_subpath) = compute_next_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            if !feasible 
                continue 
            end

            new_subpath.cost += σ_costs[(prev_node, current_node, next_node)]

            if first_node == starting_node
                first_node = current_node
            end
            if !((first_node, current_node) in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][(first_node, current_node)] = SortedDict{
                    NTuple{graph.n_nodes_extra + 1, Int},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end

            new_set = ngroute_create_set_alt(neighborhoods, current_set, next_node)
            new_key = (new_subpath.time_taken, new_set...)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][(first_node, current_node)],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, first_node, current_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging_extra
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for node in graph.N_depots_charging_extra
        delete!(base_labels[node][node], (node, node))
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][(first_node, prev_node)])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging_extra
            for (first_node, prev_node) in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][(first_node, prev_node)])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels
end

function compute_new_path(
    current_path::PathLabel,
    s::BaseSubpathLabel,
    state::NTuple{2, Int},
    next_node::Int,
    data::EVRPData,
    graph::EVRPGraph,
    christofides::Bool,
)
    # don't count initial subpath again
    if (
        next_node in graph.N_depots
        && s.time_taken == 0
    )
        return (false, nothing, nothing, nothing)
    end
    # Preventing customer 2-cycles (Christofides)
    if christofides
        if length(current_path.subpath_labels) ≥ 1
            prev_subpath = current_path.subpath_labels[end]
            if (
                prev_subpath.nodes[end-1] in graph.N_customers 
                && prev_subpath.nodes[end-1] == s.nodes[2]
            )
                return (false, nothing, nothing, nothing)
            end
        end
    end

    # time horizon and charge feasibility
    (delta, end_time, end_charge) = charge_to_specified_level(
        - state[2], # current charge
        s.charge_taken, 
        state[1], # current time
    )
    if end_time + s.time_taken + graph.min_t[next_node] > graph.T
        return (false, nothing, nothing, nothing)
    end

    new_path = copy(current_path)
    new_path.cost += s.cost
    push!(new_path.subpath_labels, s)
    new_path.served += s.served
    if length(current_path.subpath_labels) > 0
        push!(new_path.charging_actions, delta)
        new_path.cost += data.charge_cost_coeff * delta
    end

    return (true, new_path, end_time, end_charge)
end

function find_nondominated_paths_notimewindows(
    data::EVRPData,
    graph::EVRPGraph,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                NTuple{2, Int},
                SortedDict{
                    T,
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            },
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ;
    single_service::Bool = false,
    check_customers::Bool = false,
    christofides::Bool = false,
    time_limit::Float64 = Inf,
) where {T <: Tuple{Vararg{Int}}}

    start_time = time()

    keylen = check_customers ? graph.n_customers + 2 : 2
    full_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => SortedDict{
                    NTuple{keylen, Int}, 
                    PathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_depots_charging_extra
        )
        for starting_node in graph.N_depots
    )

    if check_customers
        # label key here has the following fields:
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0, -graph.B, zeros(Int, graph.n_customers)...)
    else
        key = (0, -graph.B)
    end
    unexplored_states = SortedSet{NTuple{keylen + 3, Int}}()
    for depot in graph.N_depots
        full_labels[depot][depot][depot][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                depot, # starting_node
                depot, # prev_node 
                depot, # current_node
            ),
        )
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
        if !(current_key in keys(full_labels[starting_node][current_node][prev_node]))
            continue
        end
        current_path = full_labels[starting_node][current_node][prev_node][current_key]
        for next_node in graph.N_depots_charging_extra
            for (fn, pn) in keys(base_labels[current_node][next_node])
                for s in values(base_labels[current_node][next_node][(fn, pn)])
                    # single-service requirement
                    if (
                        single_service
                        && any(s.served + current_path.served .> 1)
                    )
                        continue
                    end
                    (feasible, new_path, end_time, end_charge) = compute_new_path(
                        current_path, s, state[1:2], next_node, 
                        data, graph, 
                        christofides,
                    )
                    if !feasible
                        continue
                    end

                    new_prev_node = s.nodes[end-1]

                    if check_customers
                        new_key = (
                            end_time + s.time_taken, 
                            - (end_charge - s.charge_taken),
                            new_path.served...
                        )
                    else
                        new_key = (
                            end_time + s.time_taken, 
                            - (end_charge - s.charge_taken),
                        )
                    end
                    added = add_path_label_to_collection!(
                        full_labels[starting_node][next_node][new_prev_node],
                        new_key, new_path,
                        ;
                        verbose = false,
                    )
                    if added && next_node in graph.N_charging_extra
                        new_state = (new_key..., starting_node, new_prev_node, next_node)
                        push!(unexplored_states, new_state)
                    end
                end
            end
        end
    end

    if check_customers
        # label key here has the following fields:
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0, -graph.B, zeros(Int, graph.n_customers)...)
    else
        key = (0, -graph.B)
    end
    for depot in graph.N_depots
        push!(
            full_labels[depot][depot][depot],
            key => PathLabel(
                - κ[depot] - μ[depot],
                [
                    BaseSubpathLabel(
                        0,
                        0,
                        - κ[depot] - μ[depot],
                        [depot, depot],
                        zeros(Int, graph.n_customers),
                    )
                ],
                Int[],
                zeros(Int, graph.n_customers),
            )
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging_extra
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

function find_nondominated_paths_notimewindows_ngroute_sigma(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                NTuple{2, Int},
                Dict{
                    Tuple{Vararg{Int}},
                    SortedDict{
                        NTuple{1, Int},
                        BaseSubpathLabel,
                        Base.Order.ForwardOrdering,
                    },
                },
            },
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ;
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)

    start_time = time()
    full_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => Dict{
                    Tuple{Vararg{Int}}, 
                    SortedDict{
                        NTuple{2, Int}, 
                        PathLabel,
                        Base.Order.ForwardOrdering,
                    },
                }()
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_depots_charging_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{5, Int}}()
    for depot in graph.N_depots
        key = (0, -graph.B)
        set = (depot,)
        full_labels[depot][depot][depot][set] = SortedDict{
            NTuple{2, Int},
            PathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => PathLabel(
                0.0,
                BaseSubpathLabel[],
                NTuple{2, Int}[],
                zeros(Int, graph.n_customers),
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

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-3]
        for current_set in keys(full_labels[starting_node][current_node][prev_node])
            if !(current_key in keys(full_labels[starting_node][current_node][prev_node][current_set]))
                continue
            end
            current_path = full_labels[starting_node][current_node][prev_node][current_set][current_key]
            for next_node in graph.N_depots_charging_extra
                for (fn, pn) in keys(base_labels[current_node][next_node])
                    for set in keys(base_labels[current_node][next_node][(fn, pn)])
                        for s in values(base_labels[current_node][next_node][(fn, pn)][set])
                            # ngroute stitching subpaths check
                            (new_set, check) = ngroute_extend_partial_path_check(neighborhoods, current_set, s)
                            if !check
                                continue
                            end
                            (feasible, new_path, end_time, end_charge) = compute_new_path(
                                current_path, s, state[1:2], next_node, 
                                data, graph, christofides,
                            )
                            if !feasible
                                continue
                            end

                            new_key = (
                                end_time + s.time_taken, 
                                - (end_charge - s.charge_taken),
                            )
                            new_prev_node = s.nodes[end-1]

                            if !(new_set in keys(full_labels[starting_node][next_node][new_prev_node]))
                                full_labels[starting_node][next_node][new_prev_node][new_set] = SortedDict{
                                    NTuple{2, Int},
                                    PathLabel,
                                    Base.Order.ForwardOrdering,
                                }(Base.Order.ForwardOrdering())
                            end
                            added = add_path_label_to_collection!(
                                full_labels[starting_node][next_node][new_prev_node][new_set],
                                new_key, new_path,
                                ;
                                verbose = false,
                            )
                            if added && next_node in graph.N_charging_extra
                                new_state = (new_key..., starting_node, new_prev_node,next_node)
                                push!(unexplored_states, new_state)
                            end
                        end
                    end
                end
            end
        end
    end

    # label key here has the following fields:
    # 1) current time
    # 2) negative of current charge
    key = (0, -graph.B)
    for depot in graph.N_depots
        push!(
            full_labels[depot][depot][depot][(depot,)],
            key => PathLabel(
                - κ[depot] - μ[depot],
                [
                    BaseSubpathLabel(
                        0,
                        0,
                        - κ[depot] - μ[depot],
                        [depot, depot],
                        zeros(Int, graph.n_customers),
                    )
                ],
                Int[],
                zeros(Int, graph.n_customers),
            )
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging_extra
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

function find_nondominated_paths_notimewindows_ngroute_alt_sigma(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                NTuple{2, Int},
                SortedDict{
                    T,
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            },
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ;
    christofides::Bool = false,
    time_limit::Float64 = Inf,
) where {T <: Tuple{Vararg{Int}}}

    start_time = time()
    full_labels = Dict(
        starting_node => Dict(
            current_node => Dict(
                prev_node => SortedDict{
                    NTuple{graph.n_nodes_extra + 2, Int}, 
                    PathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
                for prev_node in graph.N_nodes_extra
            )
            for current_node in graph.N_depots_charging_extra
        )
        for starting_node in graph.N_depots
    )

    unexplored_states = SortedSet{NTuple{graph.n_nodes_extra + 5, Int}}()
    for depot in graph.N_depots
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[depot] = 1
        key = (0, -graph.B, node_labels...)
        full_labels[depot][depot][depot][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, graph.n_customers),
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

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-2]
        prev_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-3]
        if !(current_key in keys(full_labels[starting_node][current_node][prev_node]))
            continue
        end
        current_set = state[3:end-3]
        current_path = full_labels[starting_node][current_node][prev_node][current_key]
        for next_node in graph.N_depots_charging_extra
            for (fn, pn) in keys(base_labels[current_node][next_node])
                for s in values(base_labels[current_node][next_node][(fn, pn)])
                    # ngroute stitching subpaths check
                    (new_set, check) = ngroute_extend_partial_path_check_alt(neighborhoods, collect(current_set), s, graph)
                    if !check
                        continue
                    end
                    (feasible, new_path, end_time, end_charge) = compute_new_path(
                        current_path, s, state[1:2], next_node, 
                        data, graph, christofides,
                    )
                    if !feasible
                        continue
                    end

                    new_key = (
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                        new_set...,
                    )
                    new_prev_node = s.nodes[end-1]

                    added = add_path_label_to_collection!(
                        full_labels[starting_node][next_node][new_prev_node],
                        new_key, new_path,
                        ;
                        verbose = false,
                    )
                    if added && next_node in graph.N_charging_extra
                        new_state = (new_key..., starting_node, new_prev_node, next_node)
                        push!(unexplored_states, new_state)
                    end
                end
            end
        end
    end

    # label key here has the following fields:
    # 1) current time
    # 2) negative of current charge
    key = (0, -graph.B)
    for depot in graph.N_depots
        node_labels = zeros(Int, graph.n_nodes_extra)
        node_labels[depot] = 1
        key = (0, -graph.B, node_labels...,)
        full_labels[depot][depot][depot][key] = PathLabel(
            - κ[depot] - μ[depot],
            [
                BaseSubpathLabel(
                    0,
                    0,
                    - κ[depot] - μ[depot],
                    [depot, depot],
                    zeros(Int, graph.n_customers),
                )
            ],
            Int[],
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging_extra
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

unwrap_path_labels(p::PathLabel) = PathLabel[p]

function unwrap_path_labels(d::AbstractDict)
    u = PathLabel[]
    for v in values(d)
        append!(u, unwrap_path_labels(v))
    end
    return u
end

function get_negative_path_labels_from_path_labels(
    path_labels::Dict{
        Int, 
        Dict{Int, T},
    },
) where {T <: AbstractDict}
    return PathLabel[
        path_label
        for path_label in unwrap_path_labels(path_labels)
            if path_label.cost < -1e-6
    ]
end

function normalize_path_labels!(
    path_labels::Vector{PathLabel},
    graph::EVRPGraph,
)
    for path_label in path_labels
        for subpath_label in path_label.subpath_labels
            subpath_label.nodes = [
                graph.nodes_extra_to_nodes_map[node]
                for node in subpath_label.nodes
            ]
        end
    end
    return path_labels
end

function subproblem_iteration_ours(
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
    subpath_single_service::Bool = true,        
    subpath_check_customers::Bool = true,
    path_single_service::Bool = true,
    path_check_customers::Bool = true,
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)
    start_time = time()
    if ngroute && !ngroute_alt
        if length(σ) == 0
            if christofides
                base_labels_result = @timed generate_base_labels_ngroute_christofides(
                    data, graph, neighborhoods, 
                    κ, μ, ν,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            else
                base_labels_result = @timed generate_base_labels_ngroute(
                    data, graph, neighborhoods, 
                    κ, μ, ν,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            end
        else
            if christofides
                base_labels_result = @timed generate_base_labels_ngroute_sigma_christofides(
                    data, graph, neighborhoods, 
                    κ, μ, ν, σ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            else
                base_labels_result = @timed generate_base_labels_ngroute_sigma(
                    data, graph, neighborhoods, 
                    κ, μ, ν, σ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            end
        end
        base_labels_time = base_labels_result.time
        full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_sigma(
            data, graph, neighborhoods, 
            base_labels_result.value, κ, μ,
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
        full_labels_time = full_labels_result.time
    elseif ngroute && ngroute_alt
        if length(σ) == 0
            if christofides
                base_labels_result = @timed generate_base_labels_ngroute_alt_christofides(
                    data, graph, neighborhoods, 
                    κ, μ, ν,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            else
                base_labels_result = @timed generate_base_labels_ngroute_alt(
                    data, graph, neighborhoods, 
                    κ, μ, ν,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            end
        else
            if christofides
                base_labels_result = @timed generate_base_labels_ngroute_alt_sigma_christofides(
                    data, graph, neighborhoods, 
                    κ, μ, ν, σ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            else
                base_labels_result = @timed generate_base_labels_ngroute_alt_sigma(
                    data, graph, neighborhoods, 
                    κ, μ, ν, σ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
            end
        end
        base_labels_time = base_labels_result.time
        full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_alt_sigma(
            data, graph, neighborhoods, 
            base_labels_result.value, κ, μ,
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
        full_labels_time = full_labels_result.time
    elseif subpath_single_service
        base_labels_result = @timed generate_base_labels_singleservice(
            data, graph, κ, μ, ν,
            ;
            check_customers = subpath_check_customers,
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
        base_labels_time = base_labels_result.time
        full_labels_result = @timed find_nondominated_paths_notimewindows(
            data, graph, base_labels_result.value, κ, μ,
            ;
            single_service = path_single_service,
            check_customers = path_check_customers,
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
        full_labels_time = full_labels_result.time
    else 
        base_labels_result = @timed generate_base_labels_nonsingleservice(
            data, graph, κ, μ, ν,
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
        base_labels_time = base_labels_result.time
        full_labels_result = @timed find_nondominated_paths_notimewindows(
            data, graph, base_labels_result.value, κ, μ,
            ;
            single_service = path_single_service,
            check_customers = path_check_customers,
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
        full_labels_time = full_labels_result.time
    end

    negative_full_labels = get_negative_path_labels_from_path_labels(full_labels_result.value)
    normalize_path_labels!(negative_full_labels, graph)
    negative_full_labels_count = length(negative_full_labels)
    return (negative_full_labels, negative_full_labels_count, base_labels_time, full_labels_time)
end