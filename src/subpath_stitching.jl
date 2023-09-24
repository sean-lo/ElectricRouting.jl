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
        Tuple{Float64, Int},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int},
    v1::BaseSubpathLabel,
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        if all(k2 .≤ k1)
            added = false
            break
        end
        if all(k1 .≤ k2)
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_subpath_longlabel_to_collection_nodelabels!(
    collection::SortedDict{
        Tuple{Float64, Int, Vararg{Int, N}},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Vararg{Int, N}},
    v1::BaseSubpathLabel,
    ;
) where {N}
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if all(k2 .≤ k1)
            added = false
            break
        end
        # check if v1 dominates v2
        if all(k1 .≤ k2)
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_subpath_longlabel_to_collection_ngroute!(
    collection::SortedDict{
        Tuple{Float64, Int, BitVector},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, BitVector},
    v1::BaseSubpathLabel,
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[1] ≤ k1[1]
            && k2[2] ≤ k1[2]
            && all(k2[3] .≤ k1[3])
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[1] ≤ k2[1]
            && k1[2] ≤ k2[2]
            && all(k1[3] .≤ k2[3])
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_subpath_longlabel_to_collection_ngroute_alt!(
    collection::SortedDict{
        Tuple{Float64, Int, BitVector},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, BitVector},
    v1::BaseSubpathLabel,
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[1] ≤ k1[1]
            && k2[2] ≤ k1[2] 
            && all(k2[3] .≤ k1[3])
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[1] ≤ k2[1]
            && k1[2] ≤ k2[2] 
            && all(k1[3] .≤ k2[3])
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_subpath_longlabel_to_collection_ngroute_lambda!(
    collection::SortedDict{
        Tuple{Float64, Int, BitVector},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, BitVector},
    v1::BaseSubpathLabel,
    λvals::Vector{Float64},
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[2] ≤ k1[2]
            && k2[1] - sum(λvals[k2[3] .& .~k1[3]]) ≤ k1[1]
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[2] ≤ k2[2]
            && k1[1] - sum(λvals[k1[3] .& .~k2[3]]) ≤ k2[1]
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end


function add_subpath_longlabel_to_collection_ngroute_alt_lambda!(
    collection::SortedDict{
        Tuple{Float64, Int, BitVector, BitVector},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, BitVector, BitVector},
    v1::BaseSubpathLabel,
    λvals::Vector{Float64},
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[2] ≤ k1[2]
            && all(k2[4] .≤ k1[4])
            && k2[1] - sum(λvals[k2[3] .& .~k1[3]]) ≤ k1[1]
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[2] ≤ k2[2]
            && all(k1[4] .≤ k2[4])
            && k1[1] - sum(λvals[k1[3] .& .~k2[3]]) ≤ k2[1]
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end



function add_subpath_longlabel_to_collection_ngroute_lambda_lmSR3!(
    collection::SortedDict{
        Tuple{Float64, Int, BitVector},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, BitVector},
    v1::BaseSubpathLabel,
    λvals::Vector{Float64},
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[2] ≤ k1[2]
            && all(k2[3][length(λvals)+1:end] .≤ k1[3][length(λvals)+1:end])
            && k2[1] - sum(λvals[k2[3][1:length(λvals)] .& .~k1[3][1:length(λvals)]]) ≤ k1[1]
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[2] ≤ k2[2]
            && all(k1[3][length(λvals)+1:end] .≤ k2[3][length(λvals)+1:end])
            && k1[1] - sum(λvals[k1[3][1:length(λvals)] .& .~k2[3][1:length(λvals)]]) ≤ k2[1]
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end


function add_subpath_longlabel_to_collection_ngroute_alt_lambda_lmSR3!(
    collection::SortedDict{
        Tuple{Float64, Int, BitVector, BitVector},
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, BitVector, BitVector},
    v1::BaseSubpathLabel,
    λvals::Vector{Float64},
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[2] ≤ k1[2]
            && all(k2[4] .≤ k1[4])
            && all(k2[3][length(λvals)+1:end] .≤ k1[3][length(λvals)+1:end])
            && k2[1] - sum(λvals[k2[3][1:length(λvals)] .& .~k1[3][1:length(λvals)]]) ≤ k1[1]
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[2] ≤ k2[2]
            && all(k1[4] .≤ k2[4])
            && all(k1[3][length(λvals)+1:end] .≤ k2[3][length(λvals)+1:end])
            && k1[1] - sum(λvals[k1[3][1:length(λvals)] .& .~k2[3][1:length(λvals)]]) ≤ k2[1]
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end


function add_path_label_to_collection!(
    collection::SortedDict{
        Tuple{Float64, Int, Int},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Int},
    v1::PathLabel,
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if all(k2 .≤ k1)
            added = false
            break
        end
        # check if v1 dominates v2
        if all(k1 .≤ k2)
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end


function add_path_label_to_collection_nodelabels!(
    collection::SortedDict{
        Tuple{Float64, Int, Int, Vararg{Int, N}},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Int, Vararg{Int, N}},
    v1::PathLabel,
    ;
) where {N}
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if all(k2 .≤ k1)
            added = false
            break
        end
        # check if v1 dominates v2
        if all(k1 .≤ k2)
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_path_label_to_collection_ngroute!(
    collection::SortedDict{
        Tuple{Float64, Int, Int},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Int},
    v1::PathLabel,
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if all(k2 .≤ k1)
            added = false
            break
        end
        # check if v1 dominates v2
        if all(k1 .≤ k2)
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_path_label_to_collection_ngroute_alt!(
    collection::SortedDict{
        Tuple{Float64, Int, Int, BitVector},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Int, BitVector},
    v1::PathLabel,
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            all(k2[1:3] .≤ k1[1:3])
            && all(k2[4] .≤ k1[4])
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            all(k1[1:3] .≤ k2[1:3])
            && all(k1[4] .≤ k2[4])
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end


function add_path_label_to_collection_ngroute_lambda!(
    collection::SortedDict{
        Tuple{Float64, Int, Int, BitVector},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Int, BitVector},
    v1::PathLabel,
    λvals::Vector{Float64},
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            k2[2] .≤ k1[2]
            && k2[3] .≤ k1[3]
            && k2[1] - sum(λvals[k2[4] .& .~k1[4]]) ≤ k1[1]
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            k1[2] .≤ k2[2]
            && k1[3] .≤ k2[3]
            && k1[1] - sum(λvals[k1[4] .& .~k2[4]]) ≤ k2[1]
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function add_path_label_to_collection_ngroute_alt_lambda!(
    collection::SortedDict{
        Tuple{Float64, Int, Int, BitVector, BitVector},
        PathLabel,
        Base.Order.ForwardOrdering,
    },
    k1::Tuple{Float64, Int, Int, BitVector, BitVector},
    v1::PathLabel,
    λvals::Vector{Float64},
    ;
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if (
            all(k2[2:3] .≤ k1[2:3])
            && all(k2[5] .≤ k1[5])
            && k2[1] - sum(λvals[k2[4] .& .~k1[4]]) ≤ k1[1]
        )
            added = false
            break
        end
        # check if v1 dominates v2
        if (
            all(k1[2:3] .≤ k2[2:3])
            && all(k1[5] .≤ k2[5])
            && k1[1] - sum(λvals[k1[4] .& .~k2[4]]) ≤ k2[1]
        )
            pop!(collection, k2)
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

function ngroute_extend_partial_path_check(
    neighborhoods::BitMatrix,
    set::BitVector,
    nodes::Vector{Int},
)
    new_set = copy(set)
    for next_node in nodes[2:end]
        if new_set[next_node]
            return (false, set)
        end
        for node in eachindex(new_set)
            if new_set[node] && !(neighborhoods[node, next_node])
                new_set[node] = 0
            end
        end
        # new_set .&= neighborhoods[:, next_node]
        new_set[next_node] = 1
        # println("$next_node, $new_set")
    end
    return (true, new_set)
end

function compute_new_subpath(
    current_subpath::BaseSubpathLabel,
    graph::EVRPGraph,
    current_node::Int,
    next_node::Int,
    modified_costs::Matrix{Float64},
)

    # time and charge feasibility
    # if current_subpath.time_taken + graph.t[current_node, next_node] + graph.min_t[next_node] > graph.T
    #     return (false, current_subpath)
    # end 

    if current_subpath.charge_taken + graph.q[current_node, next_node] + graph.min_q[next_node] > graph.B
        return (false, current_subpath)
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


function compute_new_lambda_labels!(
    new_subpath::BaseSubpathLabel,
    current_λ_labels::BitVector,
    λvals::Vector{Float64},
    λcust::BitMatrix,
)
    next_node = new_subpath.nodes[end]
    # 1: create new λ_labels 
    new_λ_labels = current_λ_labels .⊻ λcust[:, next_node]
    # 2: modify cost of new_subpath
    new_subpath.cost -= sum(λvals[current_λ_labels .& λcust[:, next_node]])
    return new_λ_labels
end

function compute_flabels_cost_lmSR3(
    next_node::Int,
    current_λ_flabels::BitVector,
    λvals::Vector{Float64},
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    # 1: create new λ_flabels 
    λ_flabels = current_λ_flabels .& λmemory[:, next_node]
    new_λ_flabels = λ_flabels .⊻ λcust[:, next_node]
    # 2: modify cost of new_subpath
    new_cost = - sum(λvals[λ_flabels .& λcust[:, next_node]])
    return (new_λ_flabels, new_cost)
end

function compute_new_subpath_lambda_flabels_lmSR3!(
    new_subpath::BaseSubpathLabel,
    current_λ_flabels::BitVector,
    λvals::Vector{Float64},
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    (new_λ_flabels, new_cost) = compute_flabels_cost_lmSR3(
        new_subpath.nodes[end],
        current_λ_flabels,
        λvals,
        λcust,
        λmemory,
    )
    new_subpath.cost += new_cost
    return new_λ_flabels
end

function compute_new_subpath_lambda_blabels_lmSR3(
    next_node::Int,
    current_λ_blabels::BitVector,
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    return current_λ_blabels .⊻ (λmemory[:, next_node] .& λcust[:, next_node])
end

function generate_base_labels(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    elementary::Bool = true,
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)

    if elementary
        base_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Vararg{Int, graph.n_customers}},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots_charging,
                current_node in graph.N_nodes
        )        
        # label key here has the following fields:
        # 0) reduced cost
        # 1) current minimum time T_i(min)
        # 2) if applicable, whether i-th customer served
        key = (0.0, 0, zeros(Int, graph.n_customers)...,)
        unexplored_states = SortedSet{Tuple{Float64, Int, Vararg{Int, graph.n_customers + 2}}}()
    else
        base_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots_charging,
                current_node in graph.N_nodes
        )
        key = (0.0, 0,)
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Int}}()
    end 
    
    for node in graph.N_depots_charging
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
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            if (
                next_node in graph.N_customers 
                && elementary 
                && current_subpath.served[next_node] > 0
            )
                # single-service requirement
                # println("already served $next_node")
                continue
            end
            
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph, current_node, next_node, modified_costs,
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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots_charging
        if elementary
            delete!(base_labels[(node, node)], (0.0, 0, zeros(Int, graph.n_customers)...,))
        else
            delete!(base_labels[(node, node)], (0.0, 0,))
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end    
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
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
                Tuple{Float64, Vararg{Int, keylen}}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes
        )
        for starting_node in graph.N_nodes
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
            key = (modified_costs[starting_node, current_node], time_taken, served...)
        else
            key = (modified_costs[starting_node, current_node], time_taken,)
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
        for starting_node in setdiff(graph.N_nodes, new_node)
            if length(base_labels[starting_node][new_node]) == 0
                continue
            end
            for end_node in setdiff(graph.N_nodes, new_node)
                if length(base_labels[new_node][end_node]) == 0
                    continue
                end
                if true
                    for (k1, s1) in pairs(base_labels[starting_node][new_node])
                        for (k2, s2) in pairs(base_labels[new_node][end_node])
                            if time_limit < time() - start_time
                                throw(TimeLimitException())
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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
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
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, Int, Int}}()
    for node in graph.N_depots
        key = (0.0, 0,)
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
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
    for node in graph.N_charging
        key = (0.0, 0,)
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
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
        current_key = state[1:2]
        current_node_labels = state[3]
        current_fset = current_node_labels[1:graph.n_nodes]
        if starting_node in graph.N_charging
            current_bset = current_node_labels[graph.n_nodes+1:end]
        end
        if !(current_key in keys(base_labels[(starting_node, current_node)][current_node_labels]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_node_labels][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_fset) = ngroute_check_create_fset(
                neighborhoods, current_fset, next_node,
            )
            !feasible && continue
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            !feasible && continue
            if starting_node in graph.N_depots
                new_set = new_fset
            else
                new_bset = ngroute_create_bset(
                    neighborhoods, new_subpath.nodes, current_bset,
                )
                new_set = [new_fset; new_bset]
            end

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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        delete!(base_labels[(node, node)][node_labels], (0.0, 0,))
    end
    for node in graph.N_charging
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
        delete!(base_labels[(node, node)][node_labels], (0.0, 0,))
    end


    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - μ[end_node]
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
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
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
    for node in graph.N_charging
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
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
        current_fset = state[3][1:graph.n_nodes]
        if starting_node in graph.N_charging
            current_bset = state[3][graph.n_nodes+1:end]
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_fset) = ngroute_check_create_fset(
                neighborhoods, current_fset, next_node
            )
            !feasible && continue
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            !feasible && continue
            if starting_node in graph.N_depots
                new_set = new_fset
            else
                new_bset = ngroute_create_bset(
                    neighborhoods, new_subpath.nodes, current_bset,
                )
                new_set = [new_fset; new_bset]
            end
            
            new_key = (
                new_subpath.cost,
                new_subpath.time_taken, 
                new_set,
            )
            added = add_subpath_longlabel_to_collection_ngroute!(
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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        key = (0.0, 0, node_labels,)
        delete!(base_labels[(node, node)], key)
    end
    for node in graph.N_charging
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
        key = (0.0, 0, node_labels,)
        delete!(base_labels[(node, node)], key)
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels
end


function generate_base_labels_ngroute_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{NTuple{3, Int}, Float64},
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
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, BitVector, Int, Int}}()
    for node in graph.N_depots
        λ_labels = falses(length(λ))
        key = (0.0, 0, λ_labels)
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
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
    for node in graph.N_charging
        λ_labels = falses(length(λ))
        key = (0.0, 0, λ_labels)
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
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
        current_key = state[1:3]
        current_λ_labels = state[3]
        current_node_labels = state[4]
        current_fset = current_node_labels[1:graph.n_nodes]
        if starting_node in graph.N_charging
            current_bset = current_node_labels[graph.n_nodes+1:end]
        end
        if !(current_key in keys(base_labels[(starting_node, current_node)][current_node_labels]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_node_labels][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_fset) = ngroute_check_create_fset(
                neighborhoods, current_fset, next_node,
            )
            !feasible && continue
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            !feasible && continue

            new_λ_labels = compute_new_lambda_labels!(
                new_subpath, current_λ_labels, λvals, λcust,
            )

            if starting_node in graph.N_depots
                new_set = new_fset
            else
                new_bset = ngroute_create_bset(
                    neighborhoods, new_subpath.nodes, current_bset,
                )
                new_set = [new_fset; new_bset]
            end

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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        λ_labels = falses(length(λ))
        delete!(base_labels[(node, node)][node_labels], (0.0, 0, λ_labels))
    end
    for node in graph.N_charging
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
        λ_labels = falses(length(λ))
        delete!(base_labels[(node, node)][node_labels], (0.0, 0, λ_labels))
    end


    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels
end


function generate_base_labels_ngroute_alt_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{NTuple{3, Int}, Float64},
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
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
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
    for node in graph.N_charging
        λ_labels = falses(length(λ))
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
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
        current_λ_labels = state[3]
        current_fset = state[4][1:graph.n_nodes]
        if starting_node in graph.N_charging
            current_bset = state[4][graph.n_nodes+1:end]
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_fset) = ngroute_check_create_fset(
                neighborhoods, current_fset, next_node
            )
            !feasible && continue
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            !feasible && continue

            new_λ_labels = compute_new_lambda_labels!(
                new_subpath, current_λ_labels, λvals, λcust,
            )

            if starting_node in graph.N_depots
                new_set = new_fset
            else
                new_bset = ngroute_create_bset(
                    neighborhoods, new_subpath.nodes, current_bset,
                )
                new_set = [new_fset; new_bset]
            end
            
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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots
        λ_labels = falses(length(λ))
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        key = (0.0, 0, node_labels, λ_labels,)
        delete!(base_labels[(node, node)], key)
    end
    for node in graph.N_charging
        λ_labels = falses(length(λ))
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
        key = (0.0, 0, node_labels, λ_labels,)
        delete!(base_labels[(node, node)], key)
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels
end


function generate_base_labels_ngroute_lambda_lmSR3(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)

    base_labels = Dict(
        (starting_node, current_node) => Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int, BitVector},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        }()
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, BitVector, BitVector, Int, Int}}()
    for node in graph.N_depots
        λ_labels = falses(length(λ))
        key = (0.0, 0, λ_labels)
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
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
    for node in graph.N_charging
        λ_labels = falses(2 * length(λ))
        key = (0.0, 0, λ_labels)
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
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
        current_key = state[1:3]
        current_λ_labels = state[3]
        current_λ_flabels = current_λ_labels[1:length(λ)]
        current_node_labels = state[4]
        current_fset = current_node_labels[1:graph.n_nodes]
        if starting_node in graph.N_charging
            current_λ_blabels = current_λ_labels[length(λ)+1:end]
            current_bset = current_node_labels[graph.n_nodes+1:end]
        end
        if !(current_key in keys(base_labels[(starting_node, current_node)][current_node_labels]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_node_labels][current_key]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_fset) = ngroute_check_create_fset(
                neighborhoods, current_fset, next_node,
            )
            !feasible && continue
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            !feasible && continue

            new_λ_flabels = compute_new_subpath_lambda_flabels_lmSR3!(
                new_subpath, current_λ_flabels, λvals, λcust, λmemory,
            )

            if starting_node in graph.N_depots
                new_set = new_fset
                new_λ_labels = new_λ_flabels
            else
                new_bset = ngroute_create_bset(
                    neighborhoods, new_subpath.nodes, current_bset,
                )
                new_set = [new_fset; new_bset]
                new_λ_blabels = compute_new_subpath_lambda_blabels_lmSR3(
                    next_node, current_λ_blabels, λcust, λmemory,
                )
                new_λ_labels = [new_λ_flabels; new_λ_blabels]
            end

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
            added = add_subpath_longlabel_to_collection_ngroute_lambda_lmSR3!(
                base_labels[(starting_node, next_node)][new_set],
                new_key, new_subpath, λvals,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set, starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        λ_labels = falses(length(λ))
        delete!(base_labels[(node, node)][node_labels], (0.0, 0, λ_labels))
    end
    for node in graph.N_charging
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
        λ_labels = falses(2 * length(λ))
        delete!(base_labels[(node, node)][node_labels], (0.0, 0, λ_labels))
    end


    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for set in keys(base_labels[(starting_node, end_node)])
                for v in values(base_labels[(starting_node, end_node)][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels
end


function generate_base_labels_ngroute_alt_lambda_lmSR3(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix, 
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)

    base_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, BitVector, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
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
    for node in graph.N_charging
        λ_labels = falses(2 * length(λ))
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
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
        current_λ_flabels = state[3][1:length(λ)]
        current_fset = state[4][1:graph.n_nodes]
        if starting_node in graph.N_charging
            current_λ_blabels = state[3][length(λ)+1:end]
            current_bset = state[4][graph.n_nodes+1:end]
        end
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_fset) = ngroute_check_create_fset(
                neighborhoods, current_fset, next_node
            )
            !feasible && continue
            (feasible, new_subpath) = compute_new_subpath(
                current_subpath, graph,
                current_node, next_node, modified_costs,
            )
            !feasible && continue

            new_λ_flabels = compute_new_subpath_lambda_flabels_lmSR3!(
                new_subpath, current_λ_flabels, λvals, λcust, λmemory,
            )

            if starting_node in graph.N_depots
                new_set = new_fset
                new_λ_labels = new_λ_flabels
            else
                new_bset = ngroute_create_bset(
                    neighborhoods, new_subpath.nodes, current_bset,
                )
                new_set = [new_fset; new_bset]
                new_λ_blabels = compute_new_subpath_lambda_blabels_lmSR3(
                    next_node, current_λ_blabels, λcust, λmemory,
                )
                new_λ_labels = [new_λ_flabels; new_λ_blabels]
            end
            
            new_key = (
                new_subpath.cost,
                new_subpath.time_taken, 
                new_λ_labels,
                new_set,
            )
            added = add_subpath_longlabel_to_collection_ngroute_alt_lambda_lmSR3!(
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

    for starting_node in graph.N_depots_charging
        for end_node in graph.N_customers
            delete!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots
        λ_labels = falses(length(λ))
        # Forward NG-set
        node_labels = falses(graph.n_nodes)
        node_labels[node] = true
        key = (0.0, 0, node_labels, λ_labels,)
        delete!(base_labels[(node, node)], key)
    end
    for node in graph.N_charging
        λ_labels = falses(2 * length(λ))
        node_labels = falses(2 * graph.n_nodes)
        # Forward NG-set
        node_labels[node] = true
        # Backward NG-set
        node_labels[node + graph.n_nodes] = true
        key = (0.0, 0, node_labels, λ_labels,)
        delete!(base_labels[(node, node)], key)
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots_charging
            for v in values(base_labels[(starting_node, end_node)])
                v.cost = v.cost - μ[end_node]
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
)
    # don't count initial subpath again
    if (
        next_node in graph.N_depots
        && s.time_taken == 0
    )
        return (false, current_path, 0, 0)
    end

    # time horizon and charge feasibility
    (delta, end_time, end_charge) = charge_to_specified_level(
        - state[2], # current charge
        s.charge_taken, 
        state[1], # current time
    )
    if end_time + s.time_taken + graph.min_t[next_node] > graph.T
        return (false, current_path, 0, 0)
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


function compute_new_path_lambda!(
    new_path::PathLabel,
    current_path_λ_labels::BitVector,
    subpath_λ_labels::BitVector,
    λvals::Vector{Float64},
)
    new_path_λ_labels = current_path_λ_labels .⊻ subpath_λ_labels
    new_path.cost -= sum(λvals[current_path_λ_labels .& subpath_λ_labels])
    return new_path_λ_labels
end

function compute_new_path_lambda_lmSR3!(
    new_path::PathLabel,
    current_path_λ_labels::BitVector,
    nodes::Vector{Int},
    λvals::Vector{Float64},
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    new_path_λ_labels = copy(current_path_λ_labels)
    for next_node in nodes[2:end]
        (new_path_λ_labels, new_cost) = compute_flabels_cost_lmSR3(
            next_node, new_path_λ_labels,
            λvals, λcust, λmemory,
        )
        new_path.cost += new_cost
    end
    return new_path_λ_labels
end


function find_nondominated_paths_notimewindows(
    data::EVRPData,
    graph::EVRPGraph,
    base_labels::Dict{
        NTuple{2, Int},
        SortedDict{
            T,
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ;
    elementary::Bool = true,
    time_limit::Float64 = Inf,
) where {T <: Tuple{Float64, Int, Vararg{Int}}}

    start_time = time()

    if elementary
        full_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int, Vararg{Int, graph.n_customers}}, 
                PathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_depots_charging
        )
        # label key here has the following fields:
        # 0) current reduced cost
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0.0, 0, -graph.B, zeros(Int, graph.n_customers)...,)
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Vararg{Int, graph.n_customers + 2}}}()
    else
        full_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int}, 
                PathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_depots_charging
        )
        key = (0.0, 0, -graph.B,)
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Int, Int}}()
    end

    for depot in graph.N_depots
        full_labels[(depot, depot)][key] = PathLabel(
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
                depot, # current_node
            ),
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
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for s in values(base_labels[(current_node, next_node)])
                # single-service requirement
                if (
                    elementary
                    && any(s.served + current_path.served .> 1)
                )
                    continue
                end

                (feasible, new_path, end_time, end_charge) = compute_new_path(
                    current_path, s, state[2:3], next_node, 
                    data, graph,
                )
                !feasible && continue

                if elementary
                    new_key = (
                        new_path.cost,
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                        new_path.served...,
                    )
                    added = add_path_label_to_collection_nodelabels!(
                        full_labels[(starting_node, next_node)],
                        new_key, new_path,
                        ;
                    )
                else
                    new_key = (
                        new_path.cost,
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                    )
                    added = add_path_label_to_collection!(
                        full_labels[(starting_node, next_node)],
                        new_key, new_path,
                        ;
                    )
                end

                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    if elementary
        # label key here has the following fields:
        # 0) current reduced cost
        # 1) current time
        # 2) negative of current charge
        # 3) if applicable, whether i-th customer served
        key = (0.0, 0, -graph.B, zeros(Int, graph.n_customers)...,)
    else
        key = (0.0, 0, -graph.B,)
    end
    for depot in graph.N_depots
        push!(
            full_labels[(depot, depot)],
            key => PathLabel(
                - κ[depot] - μ[depot],
                [
                    BaseSubpathLabel(
                        0,
                        0,
                        - κ[depot] - μ[depot],
                        [depot, depot],
                        zeros(Int, graph.n_customers),
                    ),
                ],
                Int[],
                zeros(Int, graph.n_customers),
            )
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_ngroute(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int},
        Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()

    full_labels = Dict(
        (starting_node, current_node) => Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int, Int},
                PathLabel,
                Base.Order.ForwardOrdering,
            },
        }()
        for starting_node in graph.N_depots, 
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        key = (0.0, 0, -graph.B,)
        full_labels[(depot, depot)][node_labels] = SortedDict{
            Tuple{Float64, Int, Int},
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
                node_labels,
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
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[4]
        current_key = state[1:3]
        if !(current_key in keys(full_labels[(starting_node, current_node)][current_set]))
            continue
        end
        current_path = full_labels[(starting_node, current_node)][current_set][current_key]
        for next_node in graph.N_depots_charging
            for set in keys(base_labels[(current_node, next_node)])
                for s in values(base_labels[(current_node, next_node)][set])
                    # ngroute stitching subpaths check
                    (feasible, new_set) = ngroute_extend_partial_path_check(
                        neighborhoods, current_set, s.nodes,
                    )
                    !feasible && continue
                    (feasible, new_path, end_time, end_charge) = compute_new_path(
                        current_path, s, state[2:3], next_node, 
                        data, graph,
                    )
                    !feasible && continue

                    new_key = (
                        new_path.cost,
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                    )
                    
                    if !(new_set in keys(full_labels[(starting_node, next_node)]))
                        full_labels[(starting_node, next_node)][new_set] = SortedDict{
                            Tuple{Float64, Int, Int},
                            PathLabel,
                            Base.Order.ForwardOrdering,
                        }(Base.Order.ForwardOrdering())
                    end
                    added = add_path_label_to_collection!(
                        full_labels[(starting_node, next_node)][new_set],
                        new_key, new_path,
                        ;
                    )
                    if added && next_node in graph.N_charging
                        new_state = (new_key..., new_set, starting_node, next_node)
                        push!(unexplored_states, new_state)
                    end
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) current reduced cost
    # 1) current time
    # 2) negative of current charge
    key = (0.0, 0, -graph.B)
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        push!(
            full_labels[(depot, depot)][node_labels],
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
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_ngroute_alt(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int}, 
        SortedDict{
            Tuple{Float64, Int, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()

    full_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, Int, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        key = (0.0, 0, -graph.B)
        full_labels[(depot, depot)][(key..., node_labels)] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node_labels,
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
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[end-2]
        current_key = state[1:end-2]
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for s in values(base_labels[(current_node, next_node)])
                # ngroute stitching subpaths check
                (feasible, new_set) = ngroute_extend_partial_path_check(
                    neighborhoods, current_set, s.nodes,
                )
                !feasible && continue
                (feasible, new_path, end_time, end_charge) = compute_new_path(
                    current_path, s, state[2:3], next_node, 
                    data, graph,
                )
                !feasible && continue

                new_key = (
                    new_path.cost,
                    end_time + s.time_taken, 
                    - (end_charge - s.charge_taken),
                )
                added = add_path_label_to_collection_ngroute_alt!(
                    full_labels[(starting_node, next_node)],
                    (new_key..., new_set), new_path,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., new_set, starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) current reduced cost
    # 1) current time
    # 2) negative of current charge
    key = (0.0, 0, -graph.B)
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        full_labels[(depot, depot)][(key..., node_labels)] = PathLabel(
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
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_ngroute_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int},
        Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int, BitVector},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    λ::Dict{NTuple{3, Int}, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    λvals, _ = prepare_lambda(λ, graph.n_nodes)

    full_labels = Dict(
        (starting_node, current_node) => Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int, Int, BitVector},
                PathLabel,
                Base.Order.ForwardOrdering,
            },
        }()
        for starting_node in graph.N_depots, 
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, BitVector, Int, Int}}()
    for depot in graph.N_depots
        λ_labels = falses(length(λ))
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        key = (0.0, 0, -graph.B, λ_labels,)
        full_labels[(depot, depot)][node_labels] = SortedDict{
            Tuple{Float64, Int, Int, BitVector},
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
                node_labels,
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
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[5]
        current_key = state[1:4]
        if !(current_key in keys(full_labels[(starting_node, current_node)][current_set]))
            continue
        end
        current_path_λ_labels = state[4]
        current_path = full_labels[(starting_node, current_node)][current_set][current_key]
        for next_node in graph.N_depots_charging
            for set in keys(base_labels[(current_node, next_node)])
                for ((_, _, subpath_λ_labels), s) in pairs(base_labels[(current_node, next_node)][set])
                    # ngroute stitching subpaths check
                    (feasible, new_set) = ngroute_extend_partial_path_check(
                        neighborhoods, current_set, s.nodes,
                    )
                    !feasible && continue
                    (feasible, new_path, end_time, end_charge) = compute_new_path(
                        current_path, s, state[2:3], next_node, 
                        data, graph,
                    )
                    !feasible && continue
                    new_path_λ_labels = compute_new_path_lambda!(
                        new_path, current_path_λ_labels, subpath_λ_labels, λvals,
                    )

                    new_key = (
                        new_path.cost,
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                        new_path_λ_labels,
                    )
                    
                    if !(new_set in keys(full_labels[(starting_node, next_node)]))
                        full_labels[(starting_node, next_node)][new_set] = SortedDict{
                            Tuple{Float64, Int, Int, BitVector},
                            PathLabel,
                            Base.Order.ForwardOrdering,
                        }(Base.Order.ForwardOrdering())
                    end
                    added = add_path_label_to_collection_ngroute_lambda!(
                        full_labels[(starting_node, next_node)][new_set],
                        new_key, new_path, λvals, 
                        ;
                    )
                    if added && next_node in graph.N_charging
                        new_state = (new_key..., new_set, starting_node, next_node)
                        push!(unexplored_states, new_state)
                    end
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) current reduced cost
    # 1) current time
    # 2) negative of current charge
    # 3) one binary label per SR3 inequality added
    key = (0.0, 0, -graph.B, falses(length(λ)))
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        push!(
            full_labels[(depot, depot)][node_labels],
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
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_ngroute_alt_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int}, 
        SortedDict{
            Tuple{Float64, Int, BitVector, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    λ::Dict{NTuple{3, Int}, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    λvals, _ = prepare_lambda(λ, graph.n_nodes)

    full_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, Int, BitVector, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, BitVector, Int, Int}}()
    for depot in graph.N_depots
        λ_labels = falses(length(λ))
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        key = (0.0, 0, -graph.B)
        full_labels[(depot, depot)][(key..., λ_labels, node_labels)] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                λ_labels,
                node_labels,
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
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[5]
        current_key = state[1:5]
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        current_path_λ_labels = state[4]
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for ((_, _, subpath_λ_labels, _), s) in pairs(base_labels[(current_node, next_node)])
                # ngroute stitching subpaths check
                (feasible, new_set) = ngroute_extend_partial_path_check(
                    neighborhoods, current_set, s.nodes,
                )
                !feasible && continue
                (feasible, new_path, end_time, end_charge) = compute_new_path(
                    current_path, s, state[2:3], next_node, 
                    data, graph,
                )
                !feasible && continue
                new_path_λ_labels = compute_new_path_lambda!(
                    new_path, current_path_λ_labels, subpath_λ_labels, λvals,
                )

                new_key = (
                    new_path.cost,
                    end_time + s.time_taken, 
                    - (end_charge - s.charge_taken),
                )
                added = add_path_label_to_collection_ngroute_alt_lambda!(
                    full_labels[(starting_node, next_node)],
                    (new_key..., new_path_λ_labels, new_set), new_path, λvals,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., new_path_λ_labels, new_set, starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) current reduced cost
    # 1) current time
    # 2) negative of current charge
    key = (0.0, 0, -graph.B)
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        full_labels[(depot, depot)][(key..., falses(length(λ)), node_labels)] = PathLabel(
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
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end



function find_nondominated_paths_notimewindows_ngroute_lambda_lmSR3(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int},
        Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int, BitVector},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    λ::Dict{Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)

    full_labels = Dict(
        (starting_node, current_node) => Dict{
            BitVector,
            SortedDict{
                Tuple{Float64, Int, Int, BitVector},
                PathLabel,
                Base.Order.ForwardOrdering,
            },
        }()
        for starting_node in graph.N_depots, 
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, BitVector, Int, Int}}()
    for depot in graph.N_depots
        λ_labels = falses(length(λ))
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        key = (0.0, 0, -graph.B, λ_labels,)
        full_labels[(depot, depot)][node_labels] = SortedDict{
            Tuple{Float64, Int, Int, BitVector},
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
                node_labels,
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
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[5]
        current_key = state[1:4]
        if !(current_key in keys(full_labels[(starting_node, current_node)][current_set]))
            continue
        end
        current_path_λ_labels = state[4]
        current_path = full_labels[(starting_node, current_node)][current_set][current_key]
        for next_node in graph.N_depots_charging
            for set in keys(base_labels[(current_node, next_node)])
                for s in values(base_labels[(current_node, next_node)][set])
                    # ngroute stitching subpaths check
                    (feasible, new_set) = ngroute_extend_partial_path_check(
                        neighborhoods, current_set, s.nodes,
                    )
                    !feasible && continue
                    (feasible, new_path, end_time, end_charge) = compute_new_path(
                        current_path, s, state[2:3], next_node, 
                        data, graph,
                    )
                    !feasible && continue
                    new_path_λ_labels = compute_new_path_lambda_lmSR3!(
                        new_path, current_path_λ_labels, s.nodes, 
                        λvals, λcust, λmemory,
                    )

                    new_key = (
                        new_path.cost,
                        end_time + s.time_taken, 
                        - (end_charge - s.charge_taken),
                        new_path_λ_labels,
                    )
                    
                    if !(new_set in keys(full_labels[(starting_node, next_node)]))
                        full_labels[(starting_node, next_node)][new_set] = SortedDict{
                            Tuple{Float64, Int, Int, BitVector},
                            PathLabel,
                            Base.Order.ForwardOrdering,
                        }(Base.Order.ForwardOrdering())
                    end
                    added = add_path_label_to_collection_ngroute_lambda!(
                        full_labels[(starting_node, next_node)][new_set],
                        new_key, new_path, λvals, 
                        ;
                    )
                    if added && next_node in graph.N_charging
                        new_state = (new_key..., new_set, starting_node, next_node)
                        push!(unexplored_states, new_state)
                    end
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) current reduced cost
    # 1) current time
    # 2) negative of current charge
    # 3) one binary label per SR3 inequality added
    key = (0.0, 0, -graph.B, falses(length(λ)))
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        push!(
            full_labels[(depot, depot)][node_labels],
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
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_ngroute_alt_lambda_lmSR3(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int}, 
        SortedDict{
            Tuple{Float64, Int, BitVector, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        },
    },
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    λ::Dict{Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}, Float64},
    ;
    time_limit::Float64 = Inf,
)

    start_time = time()
    λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)

    full_labels = Dict(
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, Int, BitVector, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, BitVector, Int, Int}}()
    for depot in graph.N_depots
        λ_labels = falses(length(λ))
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        key = (0.0, 0, -graph.B)
        full_labels[(depot, depot)][(key..., λ_labels, node_labels)] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            NTuple{2, Int}[],
            zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                λ_labels,
                node_labels,
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
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[5]
        current_key = state[1:5]
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        current_path_λ_labels = state[4]
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for s in values(base_labels[(current_node, next_node)])
                # ngroute stitching subpaths check
                (feasible, new_set) = ngroute_extend_partial_path_check(
                    neighborhoods, current_set, s.nodes,
                )
                !feasible && continue
                (feasible, new_path, end_time, end_charge) = compute_new_path(
                    current_path, s, state[2:3], next_node, 
                    data, graph,
                )
                !feasible && continue
                new_path_λ_labels = compute_new_path_lambda_lmSR3!(
                    new_path, current_path_λ_labels, s.nodes,
                    λvals, λcust, λmemory,
                )

                new_key = (
                    new_path.cost,
                    end_time + s.time_taken, 
                    - (end_charge - s.charge_taken),
                )
                added = add_path_label_to_collection_ngroute_alt_lambda!(
                    full_labels[(starting_node, next_node)],
                    (new_key..., new_path_λ_labels, new_set), new_path, λvals,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., new_path_λ_labels, new_set, starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) current reduced cost
    # 1) current time
    # 2) negative of current charge
    key = (0.0, 0, -graph.B)
    for depot in graph.N_depots
        node_labels = falses(graph.n_nodes)
        node_labels[depot] = true
        full_labels[(depot, depot)][(key..., falses(length(λ)), node_labels)] = PathLabel(
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
        for end_node in graph.N_charging
            delete!(full_labels, (starting_node, end_node))
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
        NTuple{2, Int}, 
        T,
    },
) where {T <: AbstractDict}
    return PathLabel[
        path_label
        for path_label in unwrap_path_labels(path_labels)
            if path_label.cost < -1e-6
    ]
end


function subproblem_iteration_ours(
    data::EVRPData, 
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    λ::Dict{
        T,
        Float64,
    },
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
            base_labels_result = @timed generate_base_labels_ngroute(
                data, graph, neighborhoods, 
                κ, μ, ν,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            base_labels_time = base_labels_result.time
            full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute(
                data, graph, neighborhoods, 
                base_labels_result.value, κ, μ,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            full_labels_time = full_labels_result.time
        else
            if keytype(λ) == NTuple{3, Int}
                base_labels_result = @timed generate_base_labels_ngroute_lambda(
                    data, graph, neighborhoods, 
                    κ, μ, ν, λ
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                base_labels_time = base_labels_result.time
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_lambda(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            elseif keytype(λ) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
                base_labels_result = @timed generate_base_labels_ngroute_lambda_lmSR3(
                    data, graph, neighborhoods, 
                    κ, μ, ν, λ
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                base_labels_time = base_labels_result.time
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_lambda_lmSR3(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            else
                error("Unrecognized key type for λ: $(keytype(λ))")
            end
        end
    elseif ngroute && ngroute_alt
        if length(λ) == 0
            base_labels_result = @timed generate_base_labels_ngroute_alt(
                data, graph, neighborhoods, 
                κ, μ, ν,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            base_labels_time = base_labels_result.time
            full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_alt(
                data, graph, neighborhoods, 
                base_labels_result.value, κ, μ,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            full_labels_time = full_labels_result.time
        else
            if keytype(λ) == NTuple{3, Int}
                base_labels_result = @timed generate_base_labels_ngroute_alt_lambda(
                    data, graph, neighborhoods, 
                    κ, μ, ν, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                base_labels_time = base_labels_result.time
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_alt_lambda(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            elseif keytype(λ) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
                base_labels_result = @timed generate_base_labels_ngroute_alt_lambda_lmSR3(
                    data, graph, neighborhoods, 
                    κ, μ, ν, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                base_labels_time = base_labels_result.time
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_alt_lambda_lmSR3(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            else
                error("Unrecognized key type for λ: $(keytype(λ))")
            end
        end
    else 
        base_labels_result = @timed generate_base_labels(
            data, graph, κ, μ, ν,
            ;
            elementary = elementary,
            time_limit = time_limit - (time() - start_time),
        )
        base_labels_time = base_labels_result.time
        full_labels_result = @timed find_nondominated_paths_notimewindows(
            data, graph, base_labels_result.value, κ, μ,
            ;
            elementary = elementary,
            time_limit = time_limit - (time() - start_time),
        )
        full_labels_time = full_labels_result.time
    end

    negative_full_labels = get_negative_path_labels_from_path_labels(full_labels_result.value)
    negative_full_labels_count = length(negative_full_labels)
    return (negative_full_labels, negative_full_labels_count, base_labels_time, full_labels_time)
end