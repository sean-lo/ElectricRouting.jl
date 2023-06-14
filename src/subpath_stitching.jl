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


function generate_base_labels_nonsingleservice(
    G, 
    data, 
    κ,
    μ,
    ν,
    ;
)
    function add_subpath_longlabel_to_collection!(
        collection::SortedDict{
            Int,
            BaseSubpathLabel,
        },
        k1::Int,
        v1::BaseSubpathLabel,
        ;
        verbose::Bool = false,
    )
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

    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for j in data["N_customers"]
        for i in data["N_nodes"]
            modified_costs[i,j] -= ν[j]
        end
    end

    base_labels = Dict(
        start_node => Dict(
            current_node => SortedDict{Int, BaseSubpathLabel}()
            for current_node in data["N_nodes"]
        )
        for start_node in union(data["N_depots"], data["N_charging"])
    )
    for node in union(data["N_depots"], data["N_charging"])
        base_labels[node][node][0] = BaseSubpathLabel(
            0,
            0,
            0.0,
            [node],
            zeros(Int, data["n_customers"]),
        )
    end

    unexplored_states = SortedSet(
        [
            (0.0, node, node)
            for node in union(data["N_depots"], data["N_charging"])
        ]
    )

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        if !(state[1] in keys(base_labels[starting_node][current_node]))
            continue
        end
        current_subpath = base_labels[starting_node][current_node][state[1]]
        for next_node in setdiff(outneighbors(G, current_node), current_node)
            if current_subpath.charge_taken + data["q"][current_node, next_node] + data["min_q"][next_node] > data["B"]
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if next_node in data["N_customers"]
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            # if current_subpath.time_taken + data["t"][current_node, next_node] + data["min_t"][next_node] > data["T"]
            #     continue
            # end 
            new_subpath = copy(current_subpath)
            new_subpath.time_taken += data["t"][current_node, next_node]
            new_subpath.charge_taken += data["q"][current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in data["N_customers"]
                new_subpath.served[next_node] += 1
            end
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node], 
                new_subpath.time_taken, new_subpath,
                ;
                verbose = false
            )
            if added && next_node in data["N_customers"]
                next_state = (new_subpath.time_taken, starting_node, next_node)
                push!(unexplored_states, next_state)
            end
        end
    end

    for start_node in vcat(data["N_depots"], data["N_charging"])
        for end_node in data["N_customers"]
            delete!(base_labels[start_node], end_node)
        end
    end

    for start_node in data["N_depots"]
        for end_node in vcat(data["N_depots"], data["N_charging"])
            for v in values(base_labels[start_node][end_node])
                v.cost = v.cost - κ[start_node]
            end
        end
    end
    for end_node in data["N_depots"]
        for start_node in vcat(data["N_depots"], data["N_charging"])
            for v in values(base_labels[start_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    # remove self-loops with nonnegative cost
    for node in union(data["N_depots"], data["N_charging"])
        for (k, v) in pairs(base_labels[node][node])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node], k)
            end
        end
    end

    return base_labels

end

function generate_base_labels_singleservice(
    G, 
    data, 
    κ,
    μ,
    ν,
    ;
    check_customers::Bool = false,
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

    function add_subpath_longlabel_to_collection!(
        collection::SortedDict{
            Tuple{Vararg{Int}},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering
        },
        k1::Tuple{Vararg{Int}},
        v1::BaseSubpathLabel,
        ;
        verbose::Bool = false,
    )
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
                if s1.charge_taken + s2.charge_taken ≤ data["B"]
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

    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for j in data["N_customers"]
        for i in data["N_nodes"]
            modified_costs[i,j] -= ν[j]
        end
    end

    base_labels = Dict(
        start_node => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                BaseSubpathLabel,
            }()
            for current_node in data["N_nodes"]
        )
        for start_node in data["N_nodes"]
    )

    for edge in edges(G)
        start_node = edge.src
        current_node = edge.dst
        time_taken = data["t"][start_node, current_node]
        served = zeros(Int, data["n_customers"])
        if current_node in data["N_customers"]
            served[current_node] = 1
        end
        if check_customers
            key = (time_taken, served...)
        else
            key = (time_taken,)
        end
        base_labels[start_node][current_node][key] = BaseSubpathLabel(
            time_taken,
            data["q"][start_node, current_node],
            modified_costs[start_node, current_node],
            [start_node, current_node],
            served,
        )
    end

    for new_node in data["N_customers"]
        for start_node in setdiff(data["N_nodes"], new_node)
            if length(base_labels[start_node][new_node]) == 0
                continue
            end
            for end_node in setdiff(data["N_nodes"], new_node)
                if length(base_labels[new_node][end_node]) == 0
                    continue
                end
                if true
                    for (k1, s1) in pairs(base_labels[start_node][new_node])
                        for (k2, s2) in pairs(base_labels[new_node][end_node])
                            # Preventing customer 2-cycles (Christofides)
                            if s1.nodes[end-1] in data["N_customers"] && s1.nodes[end-1] == s2.nodes[2]
                                continue
                            end
                            k = k1 .+ k2
                            if !all(s1.served .+ s2.served .≤ 1)
                                continue
                            end
                            if s1.charge_taken + s2.charge_taken > data["B"]
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
                                base_labels[start_node][end_node],
                                k, s,
                            )
                        end
                    end
                else
                    merge_collections!(
                        base_labels[start_node][end_node],
                        direct_sum_of_collections(
                            base_labels[start_node][new_node], 
                            base_labels[new_node][end_node],
                        )
                    )
                end
            end
        end
    end

    for start_node in vcat(data["N_depots"], data["N_charging"])
        for end_node in data["N_customers"]
            delete!(base_labels[start_node], end_node)
        end
    end

    for start_node in data["N_depots"]
        for end_node in vcat(data["N_depots"], data["N_charging"])
            for v in values(base_labels[start_node][end_node])
                v.cost = v.cost - κ[start_node]
            end
        end
    end
    for end_node in data["N_depots"]
        for start_node in vcat(data["N_depots"], data["N_charging"])
            for v in values(base_labels[start_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    # remove self-loops with nonnegative cost
    for node in union(data["N_depots"], data["N_charging"])
        for (k, v) in pairs(base_labels[node][node])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node], k)
            end
        end
    end

    return base_labels
end


function find_nondominated_paths_notimewindows(
    data,
    base_labels,
    κ,
    μ,
    ;
    single_service::Bool = false,
    check_customers::Bool = false,
)
    function charge_to_specified_level(
        start_charge::Int, 
        desired_end_charge::Int, 
        start_time::Int, 
    )
        if desired_end_charge ≤ start_charge
            return (0, start_time, start_charge)
        end
        delta = desired_end_charge - start_charge
        end_time = start_time + delta
        return (delta, end_time, desired_end_charge)
    end

    function add_path_label_to_collection!(
        collection::SortedDict{
            Tuple{Vararg{Int}},
            PathLabel,
            Base.Order.ForwardOrdering
        },
        k1::Tuple{Vararg{Int}},
        v1::PathLabel,
        ;
        verbose::Bool = false,
    )
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

    full_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                PathLabel,
            }()
            for current_node in union(data["N_charging"], data["N_depots"])
        )
        for starting_node in data["N_depots"]
    )

    if check_customers
        # label key here has the following fields:
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0, -data["B"], zeros(Int, data["n_customers"])...)
    else
        key = (0, -data["B"])
    end
    for depot in data["N_depots"]
        full_labels[depot][depot][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Tuple{Int, Int}[],
            zeros(Int, data["n_customers"]),
        )
    end
    unexplored_states = SortedSet(
        [
            (key..., depot, depot)
            for depot in data["N_depots"]
        ]
    )

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        if !(state[1:end-2] in keys(full_labels[starting_node][current_node]))
            continue
        end
        current_path = full_labels[starting_node][current_node][state[1:end-2]]
        for next_node in union(data["N_depots"], data["N_charging"])
            for s in values(base_labels[current_node][next_node])
                if (
                    next_node in data["N_depots"]
                    && s.time_taken == 0
                )
                    continue
                end
                if (
                    single_service
                    && any(s.served + current_path.served .> 1)
                )
                    continue
                end
                # Preventing customer 2-cycles (Christofides)
                if length(current_path.subpath_labels) ≥ 1
                    prev_subpath = current_path.subpath_labels[end]
                    if (
                        prev_subpath.nodes[end-1] in data["N_customers"] 
                        && prev_subpath.nodes[end-1] == s.nodes[2]
                    )
                        continue
                    end
                end
                (delta, end_time, end_charge) = charge_to_specified_level(
                    - state[2], # current charge
                    s.charge_taken, 
                    state[1], # current time
                )
                if end_time + s.time_taken + data["min_t"][next_node] > data["T"]
                    continue
                end

                new_path = copy(current_path)
                new_path.cost += s.cost
                push!(new_path.subpath_labels, s)
                new_path.served += s.served
                if length(current_path.subpath_labels) > 0
                    push!(new_path.charging_actions, delta)
                    new_path.cost += data["charge_cost_coeff"] * delta
                end

                if check_customers
                    key = (
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                        new_path.served...
                    )
                else
                    key = (
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                    )
                end

                # println(key)
                added = add_path_label_to_collection!(
                    full_labels[starting_node][next_node],
                    key, new_path,
                )
                if added && next_node in data["N_charging"]
                    next_state = (key..., starting_node, next_node)
                    push!(unexplored_states, next_state)
                end
            end
        end
    end

    if check_customers
        # label key here has the following fields:
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0, -data["B"], zeros(Int, data["n_customers"])...)
    else
        key = (0, -data["B"])
    end
    for depot in data["N_depots"]
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
                        zeros(Int, data["n_customers"]),
                    )
                ],
                Int[],
                zeros(Int, data["n_customers"]),
            )
        )
    end

    for starting_node in data["N_depots"]
        for end_node in data["N_charging"]
            delete!(full_labels[starting_node], end_node)
        end
    end

    return full_labels

end

function get_negative_path_labels_from_path_labels(
    data, 
    path_labels::Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PathLabel,
        Base.Order.ForwardOrdering
    }}}
)
    return Dict(
        starting_node => Dict(
            end_node => SortedDict{
                Tuple{Vararg{Int}},
                PathLabel,
            }(
                key => path_label
                for (key, path_label) in path_labels[starting_node][end_node] 
                    if path_label.cost < -1e-6
            )
            for end_node in keys(path_labels[starting_node])
        )
        for starting_node in data["N_depots"]
    )
end

function subproblem_iteration_ours(
    G, data, κ, μ, ν,
    ;
    subpath_single_service::Bool = true,        
    subpath_check_customers::Bool = true,
    path_single_service::Bool = true,
    path_check_customers::Bool = true,
)
    if subpath_single_service
        base_labels_result = @timed generate_base_labels_singleservice(
            G, data, κ, μ, ν,
            ;
            check_customers = subpath_check_customers,
        )
    else 
        base_labels_result = @timed generate_base_labels_nonsingleservice(
            G, data, κ, μ, ν,
            ;
        )
    end
    base_labels_time = base_labels_result.time
    full_labels_result = @timed find_nondominated_paths_notimewindows(
        data, base_labels_result.value, κ, μ,
        ;
        single_service = path_single_service,
        check_customers = path_check_customers,
    )
    full_labels_time = full_labels_result.time
    negative_full_labels = get_negative_path_labels_from_path_labels(data, full_labels_result.value)
    negative_full_labels_count = sum(
        length(negative_full_labels[starting_node][end_node]) 
        for starting_node in data["N_depots"]
            for end_node in keys(negative_full_labels[starting_node]) 
    )
    return (negative_full_labels, negative_full_labels_count, base_labels_time, full_labels_time)
end