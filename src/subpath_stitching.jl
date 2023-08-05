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
        Base.ForwardOrdering,
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
    neighborhoods::NGRouteNeighborhood,
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
                if node in neighborhoods.x[next_node]
        ]
        push!(new_set, next_node)
        # println("$next_node, $new_set")
    end
    return (Tuple(sort(unique(new_set))), true)
end

function ngroute_extend_partial_path_check_alt(
    neighborhoods::NGRouteNeighborhood,
    set::Vector{Int},
    s::BaseSubpathLabel,
)
    new_set = copy(set)
    for next_node in s.nodes[2:end]
        if new_set[next_node] == 1
            return (nothing, false)
        end
        for node in data.N_nodes
            if new_set[node] == 1 && !(node in neighborhoods.x[next_node])
                new_set[node] = 0
            end
        end
        new_set[next_node] = 1
        # println("$next_node, $new_set")
    end
    return (new_set, true)
end

function generate_base_labels_nonsingleservice(
    data::EVRPData, 
    G::SimpleDiGraph{Int},
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{1, Int}, 
                BaseSubpathLabel,
            }()
            for current_node in data.N_nodes
        )
        for starting_node in data.N_depots_charging
    )
    unexplored_states = SortedSet{NTuple{3, Int}}()
    for node in data.N_depots_charging
        base_labels[node][node][(0,)] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, data.n_customers),
        )
        push!(unexplored_states, (0, node, node))
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
        for next_node in setdiff(outneighbors(G, current_node), current_node)
            # Preventing customer 2-cycles (Christofides)
            if christofides && next_node in data.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            # time and charge feasibility
            # if current_subpath.time_taken + data.t[current_node, next_node] + data.min_t[next_node] > data.T
            #     continue
            # end 
            if current_subpath.charge_taken + data.q[current_node, next_node] + data.min_q[next_node] > data.B
                continue
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken += data.t[current_node, next_node]
            new_subpath.charge_taken += data.q[current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in data.N_customers
                new_subpath.served[next_node] += 1
            end

            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node], 
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in data.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in data.N_depots_charging
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    # remove self-loops with nonnegative cost
    for node in data.N_depots_charging
        for (k, v) in pairs(base_labels[node][node])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node], k)
            end
        end
    end

    return base_labels

end

function generate_base_labels_singleservice(
    data::EVRPData, 
    G::SimpleDiGraph{Int},
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    check_customers::Bool = false,
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)
    function add_subpath_label_to_collection!(
        collection::SortedDict{Int, BaseSubpathLabel},
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
        labels1::SortedDict{Int, BaseSubpathLabel},
        labels2::SortedDict{Int, BaseSubpathLabel},
        ;
    )
        keys1 = collect(keys(labels1))
        keys2 = collect(keys(labels2))
    
        new = []
        for (t, cost, i, j) in sort([
            (k1 + k2, s1.cost + s2.cost, i, j)
            for (i, (k1, s1)) in enumerate(pairs(labels1)),
                (j, (k2, s2)) in enumerate(pairs(labels2))
                if s1.charge_taken + s2.charge_taken ≤ data.B
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
        new_labels = SortedDict{Int, BaseSubpathLabel}(
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
        labels1::SortedDict{Int, BaseSubpathLabel},
        labels2::SortedDict{Int, BaseSubpathLabel},
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
    modified_costs = compute_arc_modified_costs(data, ν)

    keylen = check_customers ? data.n_customers + 1 : 1
    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{keylen, Int}, 
                BaseSubpathLabel,
            }()
            for current_node in data.N_nodes
        )
        for starting_node in data.N_nodes
    )

    for edge in edges(G)
        starting_node = edge.src
        current_node = edge.dst
        time_taken = data.t[starting_node, current_node]
        served = zeros(Int, data.n_customers)
        if current_node in data.N_customers
            served[current_node] = 1
        end
        if check_customers
            key = (time_taken, served...)
        else
            key = (time_taken,)
        end
        base_labels[starting_node][current_node][key] = BaseSubpathLabel(
            time_taken,
            data.q[starting_node, current_node],
            modified_costs[starting_node, current_node],
            [starting_node, current_node],
            served,
        )
    end

    for new_node in data.N_customers
        for starting_node in setdiff(data.N_nodes, new_node)
            if length(base_labels[starting_node][new_node]) == 0
                continue
            end
            for end_node in setdiff(data.N_nodes, new_node)
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
                            if christofides && s1.nodes[end-1] in data.N_customers && s1.nodes[end-1] == s2.nodes[2]
                                continue
                            end
                            k = k1 .+ k2
                            if !all(s1.served .+ s2.served .≤ 1)
                                continue
                            end
                            if s1.charge_taken + s2.charge_taken > data.B
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

    for starting_node in data.N_depots_charging
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    # remove self-loops with nonnegative cost
    for node in data.N_depots_charging
        for (k, v) in pairs(base_labels[node][node])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node], k)
            end
        end
    end

    return base_labels
end


function generate_base_labels_ngroute(
    data::EVRPData, 
    G::SimpleDiGraph{Int},
    neighborhoods::NGRouteNeighborhood,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                Tuple{Vararg{Int}}, 
                SortedDict{
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                },
            }()
            for current_node in data.N_nodes
        )
        for starting_node in data.N_depots_charging
    )
    unexplored_states = SortedSet{NTuple{3, Int}}()
    for node in data.N_depots_charging
        base_labels[node][node][(node,)] = SortedDict{
            NTuple{1, Int},
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            (0,) => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, data.n_customers),
            )
        )
        push!(unexplored_states, (0, node, node))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        for current_set in keys(base_labels[starting_node][current_node])
            if !(current_key in keys(base_labels[starting_node][current_node][current_set]))
                continue
            end
            current_subpath = base_labels[starting_node][current_node][current_set][current_key]
            for next_node in setdiff(outneighbors(G, current_node), current_node)
                if next_node in current_set
                    # if next_node is a customer not yet visited, proceed
                    # only if one can extend current_subpath along next_node according to ng-route rules
                    continue
                end
                # Preventing customer 2-cycles (Christofides)
                if christofides && next_node in data.N_customers
                    if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                        continue
                    end
                end
                # time and charge feasibility
                # if current_subpath.time_taken + data.t[current_node, next_node] + data.min_t[next_node] > data.T
                #     continue
                # end 
                if current_subpath.charge_taken + data.q[current_node, next_node] + data.min_q[next_node] > data.B
                    continue
                end

                new_subpath = copy(current_subpath)
                new_subpath.time_taken += data.t[current_node, next_node]
                new_subpath.charge_taken += data.q[current_node, next_node]
                new_subpath.cost += modified_costs[current_node, next_node]
                push!(new_subpath.nodes, next_node)
                if next_node in data.N_customers
                    new_subpath.served[next_node] += 1
                end

                new_set = ngroute_create_set(neighborhoods, current_set, next_node)
                if !(new_set in keys(base_labels[starting_node][next_node]))
                    base_labels[starting_node][next_node][new_set] = SortedDict{
                        Tuple{Vararg{Int}},
                        BaseSubpathLabel,
                    }()
                end
                new_key = (new_subpath.time_taken,)
                added = add_subpath_longlabel_to_collection!(
                    base_labels[starting_node][next_node][new_set],
                    new_key, new_subpath,
                    ;
                    verbose = false,
                )
                if added && next_node in data.N_customers
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    for starting_node in data.N_depots_charging
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots_charging
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots_charging
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    # remove self-loops with nonnegative cost
    for node in data.N_depots_charging
        for set in keys(base_labels[node][node])
            for (k, v) in pairs(base_labels[node][node][set])
                if v.cost ≥ 0.0
                    pop!(base_labels[node][node][set], k)
                end
            end
        end
    end

    return base_labels
end

function generate_base_labels_ngroute_alt(
    data::EVRPData, 
    neighborhoods::NGRouteNeighborhood, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    christofides::Bool = false,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{data.n_nodes + 1, Int}, 
                BaseSubpathLabel,
            }()
            for current_node in data.N_nodes
        )
        for starting_node in data.N_depots_charging
    )
    unexplored_states = SortedSet{NTuple{data.n_nodes + 3, Int}}()
    for node in data.N_depots_charging
        node_labels = zeros(Int, data.n_nodes)
        node_labels[node] = 1
        key = (0, node_labels...)
        base_labels[node][node][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, data.n_customers),
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
        if !(current_key in keys(base_labels[starting_node][current_node]))
            continue
        end
        current_set = state[2:end-2]
        current_subpath = base_labels[starting_node][current_node][current_key]
        for next_node in setdiff(outneighbors(G, current_node), current_node)
            if next_node in data.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if christofides && next_node in data.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            # time and charge feasibility
            # if current_subpath.time_taken + data.t[current_node, next_node] + data.min_t[next_node] > data.T
            #     continue
            # end 
            if current_subpath.charge_taken + data.q[current_node, next_node] + data.min_q[next_node] > data.B
                continue
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken += data.t[current_node, next_node]
            new_subpath.charge_taken += data.q[current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in data.N_customers
                new_subpath.served[next_node] += 1
            end

            new_set = ngroute_create_set_alt(neighborhoods, collect(current_set), next_node)
            new_key = (new_subpath.time_taken, new_set...)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in data.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in data.N_depots_charging
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    # remove self-loops with nonnegative cost
    for node in data.N_depots_charging
        for (k, v) in pairs(base_labels[node][node])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node], k)
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
    christofides::Bool,
)
    # don't count initial subpath again
    if (
        next_node in data.N_depots
        && s.time_taken == 0
    )
        return (false, nothing, nothing, nothing)
    end
    # Preventing customer 2-cycles (Christofides)
    if christofides
        if length(current_path.subpath_labels) ≥ 1
            prev_subpath = current_path.subpath_labels[end]
            if (
                prev_subpath.nodes[end-1] in data.N_customers 
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
    if end_time + s.time_taken + data.min_t[next_node] > data.T
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
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            SortedDict{
                T,
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
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

    keylen = check_customers ? data.n_customers + 2 : 2
    full_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{keylen, Int}, 
                PathLabel,
            }()
            for current_node in data.N_depots_charging
        )
        for starting_node in data.N_depots
    )

    if check_customers
        # label key here has the following fields:
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0, -data.B, zeros(Int, data.n_customers)...)
    else
        key = (0, -data.B)
    end
    unexplored_states = SortedSet{NTuple{keylen + 2, Int}}()
    for depot in data.N_depots
        full_labels[depot][depot][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, data.n_customers),
        )
        push!(unexplored_states, (key..., depot, depot))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(full_labels[starting_node][current_node]))
            continue
        end
        current_path = full_labels[starting_node][current_node][current_key]
        for next_node in data.N_depots_charging
            for s in values(base_labels[current_node][next_node])
                # single-service requirement
                if (
                    single_service
                    && any(s.served + current_path.served .> 1)
                )
                    continue
                end
                (feasible, new_path, end_time, end_charge) = compute_new_path(
                    current_path, s, state[1:2], next_node, 
                    data, christofides,
                )
                if !feasible
                    continue
                end

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
                    full_labels[starting_node][next_node],
                    new_key, new_path,
                    ;
                    verbose = false,
                )
                if added && next_node in data.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    if check_customers
        # label key here has the following fields:
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0, -data.B, zeros(Int, data.n_customers)...)
    else
        key = (0, -data.B)
    end
    for depot in data.N_depots
        push!(
            full_labels[depot][depot],
            key => PathLabel(
                - κ[depot] - μ[depot],
                [
                    BaseSubpathLabel(
                        0,
                        0,
                        - κ[depot] - μ[depot],
                        [depot, depot],
                        zeros(Int, data.n_customers),
                    )
                ],
                Int[],
                zeros(Int, data.n_customers),
            )
        )
    end

    for starting_node in data.N_depots
        for end_node in data.N_charging
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

function find_nondominated_paths_notimewindows_ngroute(
    data::EVRPData,
    neighborhoods::NGRouteNeighborhood,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                Tuple{Vararg{Int}},
                SortedDict{
                    NTuple{1, Int},
                    BaseSubpathLabel,
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
            current_node => Dict{
                Tuple{Vararg{Int}}, 
                SortedDict{
                    NTuple{2, Int}, 
                    PathLabel,
                },
            }()
            for current_node in data.N_depots_charging
        )
        for starting_node in data.N_depots
    )

    unexplored_states = SortedSet{NTuple{4, Int}}()
    for depot in data.N_depots
        key = (0, -data.B)
        set = (depot,)
        full_labels[depot][depot][set] = SortedDict{
            NTuple{2, Int},
            PathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => PathLabel(
                0.0,
                BaseSubpathLabel[],
                NTuple{2, Int}[],
                zeros(Int, data.n_customers),
            ),
        )
        push!(unexplored_states, (key..., depot, depot))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        for current_set in keys(full_labels[starting_node][current_node])
            if !(current_key in keys(full_labels[starting_node][current_node][current_set]))
                continue
            end
            current_path = full_labels[starting_node][current_node][current_set][current_key]
            for next_node in data.N_depots_charging
                for set in keys(base_labels[current_node][next_node])
                    for s in values(base_labels[current_node][next_node][set])
                        # ngroute stitching subpaths check
                        (new_set, check) = ngroute_extend_partial_path_check(neighborhoods, current_set, s)
                        if !check
                            continue
                        end
                        (feasible, new_path, end_time, end_charge) = compute_new_path(
                            current_path, s, state[1:2], next_node, 
                            data, christofides,
                        )
                        if !feasible
                            continue
                        end

                        new_key = (
                            end_time + s.time_taken, 
                            - (end_charge - s.charge_taken),
                        )
                        if !(new_set in keys(full_labels[starting_node][next_node]))
                            full_labels[starting_node][next_node][new_set] = SortedDict{
                                Tuple{Vararg{Int}},
                                PathLabel,
                            }()
                        end
                        added = add_path_label_to_collection!(
                            full_labels[starting_node][next_node][new_set],
                            new_key, new_path,
                            ;
                            verbose = false,
                        )
                        if added && next_node in data.N_charging
                            new_state = (new_key..., starting_node, next_node)
                            push!(unexplored_states, new_state)
                        end
                    end
                end
            end
        end
    end

    # label key here has the following fields:
    # 1) current time
    # 2) negative of current charge
    key = (0, -data.B)
    for depot in data.N_depots
        push!(
            full_labels[depot][depot][(depot,)],
            key => PathLabel(
                - κ[depot] - μ[depot],
                [
                    BaseSubpathLabel(
                        0,
                        0,
                        - κ[depot] - μ[depot],
                        [depot, depot],
                        zeros(Int, data.n_customers),
                    )
                ],
                Int[],
                zeros(Int, data.n_customers),
            )
        )
    end

    for starting_node in data.N_depots
        for end_node in data.N_charging
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

function find_nondominated_paths_notimewindows_ngroute_alt(
    data::EVRPData,
    neighborhoods::NGRouteNeighborhood,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            SortedDict{
                T,
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
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
            current_node => SortedDict{
                NTuple{data.n_nodes + 2, Int}, 
                PathLabel,
            }()
            for current_node in data.N_depots_charging
        )
        for starting_node in data.N_depots
    )

    unexplored_states = SortedSet{NTuple{data.n_nodes + 4, Int}}()
    for depot in data.N_depots
        node_labels = zeros(Int, data.n_nodes)
        node_labels[depot] = 1
        key = (0, -data.B, node_labels...)
        full_labels[depot][depot][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, data.n_customers),
        )
        push!(unexplored_states, (key..., depot, depot))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(full_labels[starting_node][current_node]))
            continue
        end
        current_set = state[3:end-2]
        current_path = full_labels[starting_node][current_node][current_key]
        for next_node in data.N_depots_charging
            for s in values(base_labels[current_node][next_node])
                # ngroute stitching subpaths check
                (new_set, check) = ngroute_extend_partial_path_check_alt(neighborhoods, collect(current_set), s)
                if !check
                    continue
                end
                (feasible, new_path, end_time, end_charge) = compute_new_path(
                    current_path, s, state[1:2], next_node, 
                    data, christofides,
                )
                if !feasible
                    continue
                end

                new_key = (
                    end_time + s.time_taken, 
                    - (end_charge - s.charge_taken),
                    new_set...,
                )
                added = add_path_label_to_collection!(
                    full_labels[starting_node][next_node],
                    new_key, new_path,
                    ;
                    verbose = false,
                )
                if added && next_node in data.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 1) current time
    # 2) negative of current charge
    key = (0, -data.B)
    for depot in data.N_depots
        node_labels = zeros(Int, data.n_nodes)
        node_labels[depot] = 1
        key = (0, -data.B, node_labels...,)
        full_labels[depot][depot][key] = PathLabel(
            - κ[depot] - μ[depot],
            [
                BaseSubpathLabel(
                    0,
                    0,
                    - κ[depot] - μ[depot],
                    [depot, depot],
                    zeros(Int, data.n_customers),
                )
            ],
            Int[],
            zeros(Int, data.n_customers),
        )
    end

    for starting_node in data.N_depots
        for end_node in data.N_charging
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

function get_negative_path_labels_from_path_labels(
    data::EVRPData, 
    path_labels::Dict{
        Int, 
        Dict{
            Int, 
            SortedDict{
                T,
                PathLabel,
                Base.Order.ForwardOrdering,
            },
        },
    },
) where {T <: Tuple{Vararg{Int}}}
    return PathLabel[
        path_label
        for starting_node in data.N_depots
            for end_node in data.N_depots
                for (key, path_label) in path_labels[starting_node][end_node]
                    if path_label.cost < -1e-6
    ]
end

function get_negative_path_labels_from_path_labels_ngroute(
    data::EVRPData, 
    path_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                Tuple{Vararg{Int}}, 
                SortedDict{
                    NTuple{2, Int},
                    PathLabel,
                },
            },
        },
    },
)
    return PathLabel[
        path_label
        for starting_node in data.N_depots
            for end_node in data.N_depots
                for set in keys(path_labels[starting_node][end_node])
                    for (key, path_label) in path_labels[starting_node][end_node][set]
                        if path_label.cost < -1e-6
    ]
end

function subproblem_iteration_ours(
    data::EVRPData, 
    G::SimpleDiGraph{Int},
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    neighborhoods::Union{Nothing, NGRouteNeighborhood} = nothing,
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
        base_labels_result = @timed generate_base_labels_ngroute(
            data, G, neighborhoods, κ, μ, ν, 
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    elseif ngroute && ngroute_alt
        base_labels_result = @timed generate_base_labels_ngroute_alt(
            data, G, neighborhoods, κ, μ, ν, 
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    elseif subpath_single_service
        base_labels_result = @timed generate_base_labels_singleservice(
            data, G, κ, μ, ν,
            ;
            check_customers = subpath_check_customers,
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    else 
        base_labels_result = @timed generate_base_labels_nonsingleservice(
            data, G, κ, μ, ν,
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    end
    base_labels_time = base_labels_result.time
    if ngroute && !ngroute_alt
        full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute(
            data, neighborhoods, 
            base_labels_result.value, κ, μ,
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    elseif ngroute && ngroute_alt
        full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_alt(
            data, neighborhoods, 
            base_labels_result.value, κ, μ,
            ;
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    else
        full_labels_result = @timed find_nondominated_paths_notimewindows(
            data, base_labels_result.value, κ, μ,
            ;
            single_service = path_single_service,
            check_customers = path_check_customers,
            christofides = christofides,
            time_limit = time_limit - (time() - start_time),
        )
    end
    full_labels_time = full_labels_result.time
    if ngroute && !ngroute_alt
        negative_full_labels = get_negative_path_labels_from_path_labels_ngroute(data, full_labels_result.value)
    elseif (ngroute && ngroute_alt) || !ngroute
        negative_full_labels = get_negative_path_labels_from_path_labels(data, full_labels_result.value)
    end
    negative_full_labels_count = length(negative_full_labels)
    return (negative_full_labels, negative_full_labels_count, base_labels_time, full_labels_time)
end