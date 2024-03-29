include("utils.jl")

mutable struct BaseSubpathLabel <: Label
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

mutable struct PathLabel <: Label
    cost::Float64
    subpath_labels::Vector{BaseSubpathLabel}
    charging_actions::Vector{Int}
    charge_rebalance_indexes::Vector{Int} # ωₖ
    served::Vector{Int}
end

Base.copy(p::PathLabel) = PathLabel(
    p.cost,
    [s for s in p.subpath_labels],
    [t for t in p.charging_actions],
    [i for i in p.charge_rebalance_indexes],
    copy(p.served),
)


function ngroute_extend_partial_path_check(
    path_fset::BitVector,
    node::Int,
    subpath_fset::BitVector,
    subpath_residue::BitVector,
    subpath_bset::BitVector,
)
    new_path_inf = path_fset .& subpath_bset
    new_path_inf[node] = false
    if any(new_path_inf)
        return (false, path_fset)
    else
        return (true, subpath_fset .| (path_fset .& subpath_residue))
    end
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



function compute_new_subpath_lambda_memory_blabels_lmSR3(
    next_node::Int,
    current_λ_blabels::BitVector,
    current_λmemory_labels::BitVector,
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    new_λ_blabels = current_λ_blabels .⊻ (current_λmemory_labels .& λcust[:, next_node])
    new_λmemory_labels = current_λmemory_labels .& λmemory[:, next_node]
    return (new_λ_blabels, new_λmemory_labels)
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
                Tuple{Float64, Int, Int, BitVector},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots_charging,
                current_node in graph.N_nodes
        )        
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time taken
        # 2) charge taken
        # 3) if applicable, whether i-th customer served
        key = (0.0, 0, 0, falses(graph.n_customers),)
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, Int, Int}}()
    else
        base_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots_charging,
                current_node in graph.N_nodes
        )
        key = (0.0, 0, 0,)
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Int, Int}}()
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
        (current_key..., starting_node, current_node) = state
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
                    new_subpath.charge_taken,
                    BitVector(new_subpath.served),
                )
            else
                new_key = (
                    new_subpath.cost,
                    new_subpath.time_taken,
                    new_subpath.charge_taken,
                )
            end
            added = add_label_to_collection!(
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
            pop!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots_charging
        if elementary
            pop!(base_labels[(node, node)], (0.0, 0, 0, falses(graph.n_customers),))
        else
            pop!(base_labels[(node, node)], (0.0, 0, 0,))
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
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, Int, Int, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, Int, Int}}()
    for node in graph.N_depots_charging
        # Π: Forward NG-set
        fset = falses(graph.n_nodes)
        fset[node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        residue = copy(neighborhoods[:, node])
        # Φ: Backward NG-set
        bset = falses(graph.n_nodes)
        bset[node] = true
        ng_resources = [fset; residue; bset]
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time taken
        # 2) charge taken
        # 3) whether i-th node is in forward ng-set
        key = (0.0, 0, 0, ng_resources,)
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
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(base_labels[(starting_node, current_node)]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_key]
        (_, _, _, current_ng_resources) = current_key
        current_fset = current_ng_resources[1:graph.n_nodes]
        current_residue = current_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
        current_bset = current_ng_resources[2*graph.n_nodes+1:end]
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
                new_ng_resources = [new_fset; current_residue; current_bset]
            else
                new_residue = ngroute_create_residue(
                    neighborhoods, next_node, current_residue,
                )
                new_bset = ngroute_create_bset(
                    next_node, current_bset, current_residue,
                )
                new_ng_resources = [new_fset; new_residue; new_bset]
            end
            
            new_key = (
                new_subpath.cost,
                new_subpath.time_taken,
                new_subpath.charge_taken,
                new_ng_resources,
            )
            added = add_label_to_collection!(
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
            pop!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots_charging
        # Π: Forward NG-set
        fset = falses(graph.n_nodes)
        fset[node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        residue = copy(neighborhoods[:, node])
        # Φ: Backward NG-set
        bset = falses(graph.n_nodes)
        bset[node] = true
        ng_resources = [fset; residue; bset]
        key = (0.0, 0, 0, ng_resources,)
        pop!(base_labels[(node, node)], key)
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
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, BitVector, Int, Int, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, BitVector, Int, Int}}()
    for node in graph.N_depots_charging
        # Π: Forward NG-set
        fset = falses(graph.n_nodes)
        fset[node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        residue = copy(neighborhoods[:, node])
        # Φ: Backward NG-set
        bset = falses(graph.n_nodes)
        bset[node] = true
        ng_resources = [fset; residue; bset]
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) time taken
        # 3) charge taken
        # 4a) whether i-th node is in forward ng-set
        # 4b) whether i-th node is in forward ng-set residue
        key = (0.0, falses(length(λ)), 0, 0, ng_resources,)
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
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(base_labels[(starting_node, current_node)]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_key]
        (
            _, current_λ_labels, 
            _, _, current_ng_resources,
        ) = current_key
        current_fset = current_ng_resources[1:graph.n_nodes]
        current_residue = current_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
        current_bset = current_ng_resources[2*graph.n_nodes+1:end]
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

            (new_λ_labels, λ_cost) = compute_new_lambda_labels_cost(
                next_node, current_λ_labels, λvals, λcust,
            )
            new_subpath.cost += λ_cost

            if starting_node in graph.N_depots
                new_ng_resources = [new_fset; current_residue; current_bset]
            else
                new_residue = ngroute_create_residue(
                    neighborhoods, next_node, current_residue,
                )
                new_bset = ngroute_create_bset(
                    next_node, current_bset, current_residue,
                )
                new_ng_resources = [new_fset; new_residue; new_bset]
            end
            
            new_key = (
                new_subpath.cost,
                new_λ_labels,
                new_subpath.time_taken, 
                new_subpath.charge_taken, 
                new_ng_resources,
            )
            added = add_label_to_collection_cuts!(
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
            pop!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots_charging
        # Π: Forward NG-set
        fset = falses(graph.n_nodes)
        fset[node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        residue = copy(neighborhoods[:, node])
        # Φ: Backward NG-set
        bset = falses(graph.n_nodes)
        bset[node] = true
        ng_resources = [fset; residue; bset]
        key = (0.0, falses(length(λ)), 0, 0, ng_resources,)
        pop!(base_labels[(node, node)], key)
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
        (starting_node, current_node) => SortedDict{
            Tuple{Float64, BitVector, BitVector, BitVector, Int, Int, BitVector},
            BaseSubpathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots_charging,
            current_node in graph.N_nodes
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, BitVector, BitVector, Int, Int, BitVector, Int, Int}}()
    for node in graph.N_depots_charging
        λ_flabels = falses(length(λ))
        λ_blabels = falses(length(λ))
        λmemory_labels = copy(λmemory[:, node])
        # Π: Forward NG-set
        fset = falses(graph.n_nodes)
        fset[node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        residue = copy(neighborhoods[:, node])
        # Φ: Backward NG-set
        bset = falses(graph.n_nodes)
        bset[node] = true
        ng_resources = [fset; residue; bset]
        # label key here has the following fields:
        # 0) reduced cost
        # 1a) binary cut labels (forward)
        # 1b) binary cut labels (backward)
        # 2) time taken
        # 3) charge taken
        # 4a) whether i-th node is in forward ng-set
        # 4b) whether i-th node is in forward ng-set residue
        key = (0.0, λ_flabels, λ_blabels, λmemory_labels, 0, 0, ng_resources,)
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
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(base_labels[(starting_node, current_node)]))
            continue
        end
        current_subpath = base_labels[(starting_node, current_node)][current_key]
        (
            _, current_λ_flabels, current_λ_blabels, current_λmemory_labels, 
            _, _, current_ng_resources,
        ) = current_key
        current_fset = current_ng_resources[1:graph.n_nodes]
        current_residue = current_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
        current_bset = current_ng_resources[2*graph.n_nodes+1:end]
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

            (new_λ_flabels, λ_cost) = compute_lambda_flabels_cost_lmSR3(
                next_node, current_λ_flabels, λvals, λcust, λmemory,
            )
            new_subpath.cost += λ_cost

            if starting_node in graph.N_depots
                new_ng_resources = [new_fset; current_residue; current_bset]
                # don't update cut backward labels if subpath starts at depot
                new_λ_blabels = current_λ_blabels
                new_λmemory_labels = current_λmemory_labels
            else
                new_residue = ngroute_create_residue(
                    neighborhoods, next_node, current_residue,
                )
                new_bset = ngroute_create_bset(
                    next_node, current_bset, current_residue,
                )
                new_ng_resources = [new_fset; new_residue; new_bset]
                new_λ_blabels, new_λmemory_labels = compute_new_subpath_lambda_memory_blabels_lmSR3(
                    next_node, current_λ_blabels, current_λmemory_labels, λcust, λmemory,
                )
            end

            new_key = (
                new_subpath.cost,
                new_λ_flabels,
                new_λ_blabels,
                new_λmemory_labels,
                new_subpath.time_taken, 
                new_subpath.charge_taken, 
                new_ng_resources,
            )
            added = add_label_to_collection_cuts!(
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
            pop!(base_labels, (starting_node, end_node))
        end
    end

    for node in graph.N_depots_charging
        λ_flabels = falses(length(λ))
        λ_blabels = falses(length(λ))
        λmemory_labels = copy(λmemory[:, node])
        # Π: Forward NG-set
        fset = falses(graph.n_nodes)
        fset[node] = true
        # Ω: Nodes that if in the previous forward NG-set, stay in the next forward NG-set
        residue = copy(neighborhoods[:, node])
        # Φ: Backward NG-set
        bset = falses(graph.n_nodes)
        bset[node] = true
        ng_resources = [fset; residue; bset]
        key = (0.0, λ_flabels, λ_blabels, λmemory_labels, 0, 0, ng_resources,)
        pop!(base_labels[(node, node)], key)
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
    current_node::Int,
    current_time::Int,
    current_charge::Int,
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
    (original_charge_amount, end_time, end_charge) = charge_to_specified_level(
        current_charge,
        s.charge_taken, # desired charge  
        current_time,
    )
    end_time += s.time_taken
    end_charge -= s.charge_taken
    if end_time + graph.min_t[next_node] > graph.T
        return (false, current_path, 0, 0)
    end

    new_path = copy(current_path)
    new_path.cost += s.cost
    push!(new_path.subpath_labels, s)
    new_path.served += s.served
    if length(current_path.subpath_labels) > 0
        push!(new_path.charging_actions, original_charge_amount)
        new_path.cost += data.charge_cost_coeffs[current_node] * original_charge_amount
    end

    return (true, new_path, end_time, end_charge)
end


function compute_new_path_heterogenous_charging(
    current_path::PathLabel,
    s::BaseSubpathLabel,
    current_node::Int,
    current_time::Int,
    current_charge::Int,
    rebalancing_slacks::Vector{Int},
    next_node::Int,
    data::EVRPData,
    graph::EVRPGraph,
)

    # don't count initial subpath again
    if (
        next_node in graph.N_depots
        && s.time_taken == 0
    )
        return (false, current_path, 0, 0, rebalancing_slacks)
    end

    # time horizon and charge feasibility
    (original_charge_amount, end_time, end_charge) = charge_to_specified_level(
        current_charge,
        s.charge_taken, # desired charge  
        current_time,
    )
    end_time += s.time_taken
    end_charge -= s.charge_taken
    if end_time + graph.min_t[next_node] > graph.T
        return (false, current_path, 0, 0, rebalancing_slacks)
    end

    new_path = copy(current_path)
    new_path.cost += s.cost
    push!(new_path.subpath_labels, s)
    new_path.served += s.served

    # charge rebalancing
    if length(current_path.subpath_labels) == 0    
        return (true, new_path, end_time, end_charge, rebalancing_slacks)
    end

    charge_amount = original_charge_amount
    rho = data.charge_cost_levels[current_node]
    rebalancing_slacks_diffs = rebalancing_slacks .- vcat([0], rebalancing_slacks[1:end-1])
    for k in 1:rho-1
        ω_k = new_path.charge_rebalance_indexes[k]
        if ω_k == 0
            continue
        end
        charge_amount_k = min(charge_amount, rebalancing_slacks_diffs[k])
        charge_amount -= charge_amount_k
        new_path.charging_actions[ω_k] += charge_amount_k
        new_path.cost += data.charge_cost_levelslist[k] * charge_amount_k
        if charge_amount == 0
            break
        end
    end
    push!(new_path.charging_actions, charge_amount)
    new_path.cost += data.charge_cost_coeffs[current_node] * charge_amount

    new_path.charge_rebalance_indexes[rho] = length(current_path.subpath_labels)
    new_path.charge_rebalance_indexes[rho+1:end] .= 0

    new_rebalancing_slacks = copy(rebalancing_slacks)
    new_rebalancing_slacks[1:rho-1] .- original_charge_amount
    new_rebalancing_slacks[1:rho-1] .= max.(0, new_rebalancing_slacks[1:rho-1])
    new_rebalancing_slacks[rho:end] .= min.(data.B - s.charge_taken, data.B - current_charge)

    return (true, new_path, end_time, end_charge, new_rebalancing_slacks)
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
    subpath_λ_flabels::BitVector,
    subpath_λ_blabels::BitVector,
    subpath_λmemory_labels::BitVector,
    λvals::Vector{Float64},
)
    new_path.cost -= sum(λvals[current_path_λ_labels .& subpath_λ_blabels])
    return subpath_λ_flabels .⊻ (current_path_λ_labels .& subpath_λmemory_labels)
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
) where {T <: Union{Tuple{Float64, Int, Int, BitVector}, Tuple{Float64, Int, Int}}}

    start_time = time()

    if elementary
        full_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int, BitVector}, 
                PathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_depots_charging
        )
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time
        # 2) - charge
        # 3) if applicable, whether i-th customer served
        key = (0.0, 0, -graph.B, falses(graph.n_customers),)
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, BitVector, Int, Int}}()
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
            Int[],
            zeros(Int, data.charge_cost_nlevels),
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
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_time, 
            negative_current_charge,
        ) = current_key[1:3]
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

                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                ) = compute_new_path(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue

                if elementary
                    new_key = (
                        new_path.cost,
                        end_time, 
                        - end_charge,
                        BitVector(new_path.served),
                    )
                else
                    new_key = (
                        new_path.cost,
                        end_time, 
                        - end_charge,
                    )
                end
                added = add_label_to_collection!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path,
                    ;
                )

                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    if elementary
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time
        # 2) - charge
        # 3) if applicable, whether i-th customer served
        key = (0.0, 0, -graph.B, falses(graph.n_customers),)
    else
        key = (0.0, 0, -graph.B,)
    end
    for depot in graph.N_depots
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
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
        SortedDict{
            Tuple{Float64, Int, Int, BitVector},
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
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time
        # 2) - charge
        # 3) whether i-th node is in forward ng-set
        key = (0.0, 0, -graph.B, ng_resources,)
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
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
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_time, 
            negative_current_charge, 
            current_ng_resources,
        ) = current_key
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for (
                (_, _, _, subpath_ng_resources), 
                s,
            ) in pairs(base_labels[(current_node, next_node)])
                subpath_fset = subpath_ng_resources[1:graph.n_nodes]
                subpath_residue = subpath_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
                subpath_bset = subpath_ng_resources[2*graph.n_nodes+1:end]
                # ngroute stitching subpaths check
                (feasible, new_ng_resources) = ngroute_extend_partial_path_check(
                    current_ng_resources, 
                    current_node,
                    subpath_fset,
                    subpath_residue,
                    subpath_bset,
                )
                !feasible && continue
                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                ) = compute_new_path(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue

                new_key = (
                    new_path.cost,
                    end_time, 
                    - end_charge,
                    new_ng_resources,
                )
                added = add_label_to_collection!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) reduced cost
    # 1) time
    # 2) - charge
    # 3) whether i-th node is in forward ng-set
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        key = (0.0, 0, -graph.B, ng_resources,)
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
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
        SortedDict{
            Tuple{Float64, BitVector, Int, Int, BitVector},
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
            Tuple{Float64, BitVector, Int, Int, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) time
        # 3) - charge
        # 4) whether i-th node is in forward ng-set
        key = (0.0, falses(length(λ)), 0, -graph.B, ng_resources,)
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
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
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_path_λ_labels, 
            current_time, 
            negative_current_charge, 
            current_ng_resources,
        ) = current_key
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for (
                (_, subpath_λ_labels, _, _, subpath_ng_resources), 
                s,
            ) in pairs(base_labels[(current_node, next_node)])
                subpath_fset = subpath_ng_resources[1:graph.n_nodes]
                subpath_residue = subpath_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
                subpath_bset = subpath_ng_resources[2*graph.n_nodes+1:end]
                # ngroute stitching subpaths check
                (feasible, new_ng_resources) = ngroute_extend_partial_path_check(
                    current_ng_resources, 
                    current_node,
                    subpath_fset,
                    subpath_residue,
                    subpath_bset,
                )
                !feasible && continue
                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                ) = compute_new_path(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue
                new_path_λ_labels = compute_new_path_lambda!(
                    new_path, current_path_λ_labels, subpath_λ_labels, λvals,
                )

                new_key = (
                    new_path.cost,
                    new_path_λ_labels,
                    end_time, 
                    - end_charge,
                    new_ng_resources,
                )
                added = add_label_to_collection_cuts!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path, λvals,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) reduced cost
    # 1) binary cut labels
    # 2) time
    # 3) - charge
    # 4) whether i-th node is in forward ng-set
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        key = (0.0, falses(length(λ)), 0, -graph.B, ng_resources,)
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
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
        SortedDict{
            Tuple{Float64, BitVector, BitVector, BitVector, Int, Int, BitVector},
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
            Tuple{Float64, BitVector, Int, Int, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, BitVector, Int, Int}}()
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) time
        # 3) - charge
        # 4) whether i-th node is in forward ng-set
        key = (0.0, falses(length(λ)), 0, -graph.B, ng_resources,)
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
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
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_path_λ_labels, 
            current_time, 
            negative_current_charge, 
            current_ng_resources,
        ) = current_key
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for (
                (_, subpath_λ_flabels, subpath_λ_blabels, subpath_λmemory_labels, 
                _, _, subpath_ng_resources), 
                s,
            ) in pairs(base_labels[(current_node, next_node)])
                subpath_fset = subpath_ng_resources[1:graph.n_nodes]
                subpath_residue = subpath_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
                subpath_bset = subpath_ng_resources[2*graph.n_nodes+1:end]
                # ngroute stitching subpaths check
                (feasible, new_ng_resources) = ngroute_extend_partial_path_check(
                    current_ng_resources, 
                    current_node,
                    subpath_fset,
                    subpath_residue,
                    subpath_bset,
                )
                !feasible && continue
                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                ) = compute_new_path(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue
                new_path_λ_labels = compute_new_path_lambda_lmSR3!(
                    new_path, current_path_λ_labels, 
                    subpath_λ_flabels, subpath_λ_blabels, subpath_λmemory_labels, 
                    λvals,
                )

                new_key = (
                    new_path.cost,
                    new_path_λ_labels, 
                    end_time, 
                    - end_charge,
                    new_ng_resources,
                )
                added = add_label_to_collection_cuts!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path, λvals,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) reduced cost
    # 1) binary cut labels
    # 2) time
    # 3) - charge
    # 4) whether i-th node is in forward ng-set
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        key = (0.0, falses(length(λ)), 0, -graph.B, ng_resources,)
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end



function find_nondominated_paths_notimewindows_heterogenous_charging(
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
) where {T <: Union{Tuple{Float64, Int, Int, BitVector}, Tuple{Float64, Int, Int}}}

    start_time = time()

    if elementary
        full_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int, Vector{Int}, BitVector}, 
                PathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_depots_charging
        )
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time
        # 2) - charge
        # 3) - rebalancing slacks
        # 4) if applicable, whether i-th customer served
        key = (
            0.0, 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            falses(graph.n_customers),
        )
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Vector{Int}, BitVector, Int, Int}}()
    else
        full_labels = Dict(
            (starting_node, current_node) => SortedDict{
                Tuple{Float64, Int, Int, Vector{Int}}, 
                PathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for starting_node in graph.N_depots,
                current_node in graph.N_depots_charging
        )
        key = (
            0.0, 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
        )
        unexplored_states = SortedSet{Tuple{Float64, Int, Int, Vector{Int}, Int, Int}}()
    end

    for depot in graph.N_depots
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
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
        (current_key..., starting_node, current_node) = state
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_time, 
            negative_current_charge, 
            negative_rebalancing_slacks,
        ) = current_key[1:4]
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

                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                    new_rebalancing_slacks,
                ) = compute_new_path_heterogenous_charging(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    - negative_rebalancing_slacks, 
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue

                if elementary
                    new_key = (
                        new_path.cost,
                        end_time, 
                        - end_charge,
                        - new_rebalancing_slacks,
                        BitVector(new_path.served),
                    )
                else
                    new_key = (
                        new_path.cost,
                        end_time, 
                        - end_charge,
                        - new_rebalancing_slacks,
                    )
                end
                added = add_label_to_collection!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path,
                    ;
                )

                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    if elementary
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time
        # 2) - charge
        # 3) - rebalancing slacks
        # 4) if applicable, whether i-th customer served
        key = (
            0.0, 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            falses(graph.n_customers),
        )
    else
        key = (
            0.0, 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1),
        )
    end
    for depot in graph.N_depots
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_heterogenous_charging_ngroute(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int}, 
        SortedDict{
            Tuple{Float64, Int, Int, BitVector},
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
            Tuple{Float64, Int, Int, Vector{Int}, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, Int, Int, Vector{Int}, BitVector, Int, Int}}()
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) time
        # 2) - charge
        # 3) - rebalancing slacks
        # 4) whether i-th node is in forward ng-set
        key = (
            0.0, 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            ng_resources,
        )
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
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
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_time, 
            negative_current_charge, 
            negative_rebalancing_slacks, 
            current_ng_resources,
        ) = current_key
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for (
                (_, _, _, subpath_ng_resources), 
                s,
            ) in pairs(base_labels[(current_node, next_node)])
                subpath_fset = subpath_ng_resources[1:graph.n_nodes]
                subpath_residue = subpath_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
                subpath_bset = subpath_ng_resources[2*graph.n_nodes+1:end]
                # ngroute stitching subpaths check
                (feasible, new_ng_resources) = ngroute_extend_partial_path_check(
                    current_ng_resources, 
                    current_node,
                    subpath_fset,
                    subpath_residue,
                    subpath_bset,
                )
                !feasible && continue
                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                    new_rebalancing_slacks,
                ) = compute_new_path_heterogenous_charging(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    - negative_rebalancing_slacks,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue

                new_key = (
                    new_path.cost,
                    end_time, 
                    - end_charge,
                    - new_rebalancing_slacks,
                    new_ng_resources,
                )
                added = add_label_to_collection!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) reduced cost
    # 1) time
    # 2) - charge
    # 3) - rebalancing slacks
    # 4) whether i-th node is in forward ng-set
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        key = (
            0.0, 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            ng_resources,
        )
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end


function find_nondominated_paths_notimewindows_heterogenous_charging_ngroute_lambda(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int}, 
        SortedDict{
            Tuple{Float64, BitVector, Int, Int, BitVector},
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
            Tuple{Float64, BitVector, Int, Int, Vector{Int}, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, Vector{Int}, BitVector, Int, Int}}()
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) time
        # 3) - charge
        # 4) - rebalancing slacks
        # 4) whether i-th node is in forward ng-set
        key = (
            0.0, falses(length(λ)), 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            ng_resources,
        )
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
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
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_path_λ_labels, 
            current_time, 
            negative_current_charge, 
            negative_rebalancing_slacks, 
            current_ng_resources,
        ) = current_key
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for (
                (_, subpath_λ_labels, _, _, subpath_ng_resources), 
                s,
            ) in pairs(base_labels[(current_node, next_node)])
                subpath_fset = subpath_ng_resources[1:graph.n_nodes]
                subpath_residue = subpath_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
                subpath_bset = subpath_ng_resources[2*graph.n_nodes+1:end]
                # ngroute stitching subpaths check
                (feasible, new_ng_resources) = ngroute_extend_partial_path_check(
                    current_ng_resources, 
                    current_node,
                    subpath_fset,
                    subpath_residue,
                    subpath_bset,
                )
                !feasible && continue
                (
                    feasible, 
                    new_path, 
                    end_time, 
                    end_charge, 
                    new_rebalancing_slacks,
                ) = compute_new_path_heterogenous_charging(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    - negative_rebalancing_slacks,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue
                new_path_λ_labels = compute_new_path_lambda!(
                    new_path, current_path_λ_labels, subpath_λ_labels, λvals,
                )

                new_key = (
                    new_path.cost,
                    new_path_λ_labels,
                    end_time, 
                    - end_charge,
                    - new_rebalancing_slacks,
                    new_ng_resources,
                )
                added = add_label_to_collection_cuts!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path, λvals,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) reduced cost
    # 1) binary cut labels
    # 2) time
    # 3) - charge
    # 4) - rebalancing slacks
    # 4) whether i-th node is in forward ng-set
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        key = (
            0.0, falses(length(λ)), 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            ng_resources,
        )
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
        end
    end

    return full_labels

end



function find_nondominated_paths_notimewindows_heterogenous_charging_ngroute_lambda_lmSR3(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    base_labels::Dict{
        NTuple{2, Int}, 
        SortedDict{
            Tuple{Float64, BitVector, BitVector, Int, Int, BitVector},
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
            Tuple{Float64, BitVector, Int, Int, Vector{Int}, BitVector}, 
            PathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_depots_charging
    )

    unexplored_states = SortedSet{Tuple{Float64, BitVector, Int, Int, Vector{Int}, BitVector, Int, Int}}()
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        # label key here has the following fields:
        # 0) reduced cost
        # 1) binary cut labels
        # 2) time
        # 3) - charge
        # 4) - rebalancing slacks
        # 5) whether i-th node is in forward ng-set
        key = (
            0.0, falses(length(λ)), 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            ng_resources,
        )
        full_labels[(depot, depot)][key] = PathLabel(
            0.0,
            BaseSubpathLabel[],
            Int[],
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
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
        if !(current_key in keys(full_labels[(starting_node, current_node)]))
            continue
        end
        (
            _, current_path_λ_labels, 
            current_time, 
            negative_current_charge, 
            negative_rebalancing_slacks, 
            current_ng_resources,
        ) = current_key
        current_path = full_labels[(starting_node, current_node)][current_key]
        for next_node in graph.N_depots_charging
            for (
                (_, subpath_λ_flabels, subpath_λ_blabels, subpath_λmemory_labels, 
                _, _, subpath_ng_resources), 
                s,
            ) in pairs(base_labels[(current_node, next_node)])                
                subpath_fset = subpath_ng_resources[1:graph.n_nodes]
                subpath_residue = subpath_ng_resources[graph.n_nodes+1:2*graph.n_nodes]
                subpath_bset = subpath_ng_resources[2*graph.n_nodes+1:end]
                # ngroute stitching subpaths check
                (feasible, new_ng_resources) = ngroute_extend_partial_path_check(
                    current_ng_resources, 
                    current_node,
                    subpath_fset,
                    subpath_residue,
                    subpath_bset,
                )
                !feasible && continue
                (feasible, new_path, end_time, end_charge, new_rebalancing_slacks) = compute_new_path_heterogenous_charging(
                    current_path, 
                    s, 
                    current_node,
                    current_time,
                    - negative_current_charge,
                    - negative_rebalancing_slacks,
                    next_node, 
                    data, 
                    graph,
                )
                !feasible && continue
                new_path_λ_labels = compute_new_path_lambda_lmSR3!(
                    new_path, current_path_λ_labels, 
                    subpath_λ_flabels, subpath_λ_blabels, subpath_λmemory_labels, 
                    λvals,
                )

                new_key = (
                    new_path.cost,
                    new_path_λ_labels, 
                    end_time, 
                    - end_charge,
                    - new_rebalancing_slacks,
                    new_ng_resources,
                )
                added = add_label_to_collection_cuts!(
                    full_labels[(starting_node, next_node)],
                    new_key, new_path, λvals,
                    ;
                )
                if added && next_node in graph.N_charging
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    # label key here has the following fields:
    # 0) reduced cost
    # 1) binary cut labels
    # 2) time
    # 3) - charge
    # 4) - rebalancing slacks
    # 5) whether i-th node is in forward ng-set
    for depot in graph.N_depots
        ng_resources = falses(graph.n_nodes)
        ng_resources[depot] = true
        key = (
            0.0, falses(length(λ)), 0, -graph.B, 
            zeros(Int, data.charge_cost_nlevels - 1), 
            ng_resources,
        )
        full_labels[(depot, depot)][key] = PathLabel(
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
            zeros(Int, data.charge_cost_nlevels),
            zeros(Int, graph.n_customers),
        )
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_charging
            pop!(full_labels, (starting_node, end_node))
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
    charge_cost_heterogenous::Bool = false,
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    elementary::Bool = true,
    time_limit::Float64 = Inf,
) where {T}
    start_time = time()
    if ngroute
        if length(λ) == 0
            base_labels_result = @timed generate_base_labels_ngroute(
                data, graph, neighborhoods, 
                κ, μ, ν,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            base_labels_time = base_labels_result.time
            if charge_cost_heterogenous
                full_labels_result = @timed find_nondominated_paths_notimewindows_heterogenous_charging_ngroute(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            else
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            end
        elseif keytype(λ) == NTuple{3, Int}
            base_labels_result = @timed generate_base_labels_ngroute_lambda(
                data, graph, neighborhoods, 
                κ, μ, ν, λ,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            base_labels_time = base_labels_result.time
            if charge_cost_heterogenous
                full_labels_result = @timed find_nondominated_paths_notimewindows_heterogenous_charging_ngroute_lambda(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            else
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_lambda(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            end
        elseif keytype(λ) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
            base_labels_result = @timed generate_base_labels_ngroute_lambda_lmSR3(
                data, graph, neighborhoods, 
                κ, μ, ν, λ,
                ;
                time_limit = time_limit - (time() - start_time),
            )
            base_labels_time = base_labels_result.time
            if charge_cost_heterogenous
                full_labels_result = @timed find_nondominated_paths_notimewindows_heterogenous_charging_ngroute_lambda_lmSR3(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            else
                full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_lambda_lmSR3(
                    data, graph, neighborhoods, 
                    base_labels_result.value, κ, μ, λ,
                    ;
                    time_limit = time_limit - (time() - start_time),
                )
                full_labels_time = full_labels_result.time
            end
        else
            error("Unrecognized key type for λ: $(keytype(λ))")
        end
    else 
        base_labels_result = @timed generate_base_labels(
            data, graph, κ, μ, ν,
            ;
            elementary = elementary,
            time_limit = time_limit - (time() - start_time),
        )
        base_labels_time = base_labels_result.time
        if charge_cost_heterogenous
            full_labels_result = @timed find_nondominated_paths_notimewindows_heterogenous_charging(
                data, graph, base_labels_result.value, κ, μ,
                ;
                elementary = elementary,
                time_limit = time_limit - (time() - start_time),
            )
            full_labels_time = full_labels_result.time
        else
            full_labels_result = @timed find_nondominated_paths_notimewindows(
                data, graph, base_labels_result.value, κ, μ,
                ;
                elementary = elementary,
                time_limit = time_limit - (time() - start_time),
            )
            full_labels_time = full_labels_result.time
        end
    end

    negative_full_labels = get_negative_path_labels_from_path_labels(full_labels_result.value)
    negative_full_labels_count = length(negative_full_labels)
    return (negative_full_labels, negative_full_labels_count, base_labels_time, full_labels_time)
end