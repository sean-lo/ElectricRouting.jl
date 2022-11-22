using CompositeStructs
using Suppressor

using JuMP
using Gurobi

Base.@kwdef mutable struct Subpath
    n_customers::Int
    starting_node::Int
    starting_time::Float64
    starting_charge::Float64
    current_node::Int = starting_node
    arcs::Vector{Tuple} = []
    time::Float64 = starting_time
    charge::Float64 = starting_charge
    served::BitVector = falses(n_customers)
    delta_time::Float64 = 0.0
    delta_charge::Float64 = 0.0
    end_time::Float64 = starting_time
    end_charge::Float64 = starting_charge
    round_time::Float64 = starting_time
    round_charge::Float64 = starting_charge
end

Base.copy(s::Subpath) = Subpath(
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    time = s.time,
    charge = s.charge,
    served = copy(s.served),
    delta_time = s.delta_time,
    delta_charge = s.delta_charge,
    end_time = s.end_time,
    end_charge = s.end_charge,
    round_time = s.round_time,
    round_charge = s.round_charge,
)

Base.show(io::IO, s::Subpath) = print(io, """Subpath:
($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.round_time), $(s.round_charge))=
arcs:   $(s.arcs)
served: $(s.served)
now:    ($(s.time), $(s.charge))
delta:  ($(s.delta_time), $(s.delta_charge))
end:    ($(s.end_time), $(s.end_charge))
round:  ($(s.round_time), $(s.round_charge))
""")

Base.isequal(s1::Subpath, s2::Subpath) = begin 
    (
        s1.n_customers == s2.n_customers
        && s1.starting_node == s2.starting_node
        && s1.starting_time == s2.starting_time
        && s1.starting_charge == s2.starting_charge
        && s1.current_node == s2.current_node
        && s1.arcs == s2.arcs
        && s1.time == s2.time
        && s1.charge == s2.charge
        && s1.served == s2.served
        && s1.delta_time == s2.delta_time
        && s1.delta_charge == s2.delta_charge
        && s1.end_time == s2.end_time
        && s1.end_charge == s2.end_charge
        && s1.round_time == s2.round_time
        && s1.round_charge == s2.round_charge
    )
end

@composite Base.@kwdef mutable struct SubpathWithCost 
    Subpath...
    cost::Float64 = 0.0
end

Base.show(io::IO, s::SubpathWithCost) = print(io, """SubpathWithCost:
($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.round_time), $(s.round_charge))
cost:   $(s.cost)
arcs:   $(s.arcs)
served: $(s.served)
now:    ($(s.time), $(s.charge))
delta:  ($(s.delta_time), $(s.delta_charge))
end:    ($(s.end_time), $(s.end_charge))
round:  ($(s.round_time), $(s.round_charge))
""")

Base.copy(s::SubpathWithCost) = SubpathWithCost(
    cost = s.cost,
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    time = s.time,
    charge = s.charge,
    served = copy(s.served),
    delta_time = s.delta_time,
    delta_charge = s.delta_charge,
    end_time = s.end_time,
    end_charge = s.end_charge,
    round_time = s.round_time,
    round_charge = s.round_charge,
)

Base.isequal(s1::SubpathWithCost, s2::SubpathWithCost) = begin 
    (
        s1.cost == s2.cost
        && s1.n_customers == s2.n_customers
        && s1.starting_node == s2.starting_node
        && s1.starting_time == s2.starting_time
        && s1.starting_charge == s2.starting_charge
        && s1.current_node == s2.current_node
        && s1.arcs == s2.arcs
        && s1.time == s2.time
        && s1.charge == s2.charge
        && s1.served == s2.served
        && s1.delta_time == s2.delta_time
        && s1.delta_charge == s2.delta_charge
        && s1.end_time == s2.end_time
        && s1.end_charge == s2.end_charge
        && s1.round_time == s2.round_time
        && s1.round_charge == s2.round_charge
    )
end

function construct_graph(data)
    G = SimpleWeightedDiGraph(data["n_nodes"])
    for (i, j) in keys(data["A"])
        if i == j
            add_edge!(G, i, i, 1)
        else 
            add_edge!(G, i, j, data["t"][i,j])
        end
    end
    return G
end

function enumerate_subpaths(
    starting_node, 
    starting_time, 
    starting_charge,
    G,
    data,
)
    """
    Enumerates all subpaths starting from `starting_node`
    (which must be a charging station or depot), 
    at time `starting_time` with charge `starting_charge`,
    visiting some pickups and dropoffs, and arriving at a
    charging station or depot (i) empty, (ii) with feasible charge 
    and (iii) in feasible time.
    """

    # Initializes queue of subpaths
    subpaths = [
        Subpath(
            n_customers = data["n_customers"],
            starting_node = starting_node,
            starting_time = starting_time,
            starting_charge = starting_charge,
        )
    ]
    
    qmin = minimum(
        data["q"][
            union(data["N_depots"], data["N_charging"]),
            1:end
        ],
        dims = 1
    )

    # Processes subpath s in subpaths one at a time
    completed_subpaths = []
    while length(subpaths) > 0
        s = popfirst!(subpaths)
        # If current node of subpath is a charging station or depot, end the subpath
        if (s.current_node in union(data["N_charging"], data["N_depots"])
            && length(s.arcs) > 0
        )
            # Cleanup actions
            s.end_time = s.time
            s.round_time = s.time
            s.end_charge = s.charge
            s.round_charge = s.charge
            push!(completed_subpaths, s)
            continue
        end
        # Iterate over out-neighbors of current node
        for j in outneighbors(G, s.current_node)
            if j in data["N_pickups"] && s.served[j]
                continue
            end
            if j in data["N_pickups"]
                new_node = j + data["n_customers"]
                # two-step lookahead if j is a pickup
                new_time_1 = max(
                    s.time + data["t"][s.current_node,j],
                    data["α"][j],
                )
                feasible_timewindow = (new_time_1 ≤ data["β"][j])
                new_time = max(
                    new_time_1 + data["t"][j,new_node],
                    data["α"][new_node],
                )
                feasible_timewindow = feasible_timewindow && (new_time ≤ data["β"][new_node])
                new_charge = s.charge - data["q"][s.current_node,j] - data["q"][j,new_node]
                feasible_charge = (new_charge ≥ qmin[new_node])
            else
                new_node = j
                new_time = max(
                    s.time + data["t"][s.current_node,new_node],
                    data["α"][new_node],
                )
                feasible_timewindow = (new_time ≤ data["β"][new_node])
                new_charge = s.charge - data["q"][s.current_node,new_node]
                feasible_charge = (new_charge ≥ qmin[new_node])
            end

            if !feasible_timewindow
                continue
            end
            if !feasible_charge
                continue
            end

            s_j = copy(s)
            s_j.current_node = new_node
            if j in data["N_pickups"]
                push!(s_j.arcs, (s.current_node, j), (j, new_node))
                s_j.served[j] = true
            else 
                push!(s_j.arcs, (s.current_node, new_node))
            end
            s_j.time = new_time
            s_j.charge = new_charge
            push!(subpaths, s_j)
        end
    end
    println("Enumerated $(length(completed_subpaths)) subpaths.")

    # Split completed subpaths into groups according to:
    # (i) starting_node (ii) ending_node (iii) customers served
    split_completed_subpaths = Dict(
        (starting_node, j, BitVector(x)) => []
        for j in union(data["N_charging"], data["N_depots"])
        for x in Iterators.product(repeat([[true, false]], data["n_customers"])...)
    )
    for s in completed_subpaths 
        push!(
            split_completed_subpaths[(s.starting_node, s.current_node, s.served)],
            s
        )
    end

    # Within each group: remove dominated subpaths
    # Domination criterion: 
    # (i)   s1.time ≤ s2.time, and 
    # (ii)  s1.charge ≥ s2.charge
    for ((starting_node, ending_node, served), group_subpaths) in pairs(split_completed_subpaths)
        if length(group_subpaths) ≤ 1
            continue
        end
        keep = trues(length(group_subpaths))
        for (i1, i2) in combinations(1:length(group_subpaths), 2)
            if !(keep[i1] && keep[i2])
                continue
            end
            s1 = group_subpaths[i1]
            s2 = group_subpaths[i2]
            if s1.time ≤ s2.time && s1.charge ≥ s2.charge
                keep[i2] = false
            elseif s1.time ≥ s2.time && s1.charge ≤ s2.charge
                keep[i1] = false
            end
        end
        split_completed_subpaths[(starting_node, ending_node, served)] = group_subpaths[keep] 
    end

    # Across two groups: remove dominated subpaths
    # (starting_node, ending_node, served_1),
    # (starting_node, ending_node, served_2),
    # with served_2 ≥ served_1, served_2 ≂̸ served_1
    # and for some ending_node
    # Question: is this step ever necessary?
    for diff in 1:data["n_customers"]
        for n1 in 0:(data["n_customers"]-diff)
            n2 = n1 + diff
            for n1_cust in combinations(1:data["n_customers"], n1)
                n1_served = BitVector([i in n1_cust for i in 1:data["n_customers"]])
                for n2_new_cust in combinations(setdiff(1:data["n_customers"], n1_cust), diff)
                    n2_cust = union(n1_cust, n2_new_cust)
                    n2_served = BitVector([i in n2_cust for i in 1:data["n_customers"]])
                    for j in union(data["N_depots"], data["N_charging"])
                        large_group = split_completed_subpaths[(starting_node, j, n2_served)]
                        small_group = split_completed_subpaths[(starting_node, j, n1_served)]
                        if !(length(large_group) ≥ 1 && length(small_group) ≥ 1)
                            continue
                        end
                        small_group_keep = trues(length(small_group))
                        for s_large in large_group
                            for (ind, s_small) in enumerate(small_group)
                                if !small_group_keep[ind]
                                    continue
                                end
                                if s_small.time ≥ s_large.time && s_small.charge ≤ s_large.charge
                                    small_group_keep[ind] = false
                                end
                            end
                        end
                        split_completed_subpaths[(starting_node, j, n1_served)] = small_group[small_group_keep]
                    end
                end
            end
        end
    end

    nondominated_subpaths = vcat(collect(values(split_completed_subpaths))...)
    println("Retaining $(length(nondominated_subpaths)) non-dominated subpaths.")
    return completed_subpaths, nondominated_subpaths
end

function dceil(
    x::Float64,
    points,
)
    return points[searchsortedfirst(points, x)]
end

function dfloor(
    x::Float64,
    points,
)
    return points[searchsortedlast(points, x)]
end

function dceilall(
    x::Float64,
    points,
)
    return points[searchsortedfirst(points, x):end]
end

function dfloorall(
    x::Float64,
    points,
)
    return points[1:searchsortedlast(points, x)]
end

function generate_charging_options(
    starting_time,
    starting_charge,
    data,
    T_range,
    B_range,
)
    # Version 2: charge by fixed time intervals inside T_range
    # If maximum time reached, charge to maximum time and get the corresponding charge
    # If maximum charge reached, charge to the next higher time (in T_range) and get to max charge
    max_end_time_unbounded = starting_time + (data["B"] - starting_charge) / data["μ"]
    if max_end_time_unbounded > data["T"]
        max_end_time = data["T"]
    else
        # This dceil means that if you charge to full charge,
        # you have to wait until the next discretized time
        max_end_time = dceil(max_end_time_unbounded, T_range)
    end
    end_times = [
        t for t in T_range
        if starting_time < t ≤ max_end_time
    ]
    round_times = end_times
    delta_times = end_times .- starting_time

    delta_charges = delta_times .* data["μ"]
    end_charges = delta_charges .+ starting_charge
    # the possible charges at which you end charging
    round_charges = [dfloor(b, B_range) for b in end_charges]
    return collect(Iterators.zip(
        delta_times, delta_charges, 
        end_times, end_charges, 
        round_times, round_charges,
    ))
end

function enumerate_subpaths_withcharge(
    starting_node, 
    starting_time, 
    starting_charge, 
    G, 
    data,
    T_range,
    B_range,
)
    """
    `nondominated_subpaths_withcharge`: Vector{Subpath}, 
    generated from a single Subpath with additional discretized charging options
    """
    _, nondominated_subpaths = enumerate_subpaths(starting_node, starting_time, starting_charge, G, data)
    
    nondominated_subpaths_withcharge = []
    for s in nondominated_subpaths
        if s.current_node in data["N_depots"]
            push!(nondominated_subpaths_withcharge, s)
        else
            for (delta_time, delta_charge, 
                end_time, end_charge, 
                round_time, round_charge) in generate_charging_options(
                s.time, s.charge, data, T_range, B_range,
            )
                s_withcharge = copy(s)
                s_withcharge.delta_time = delta_time
                s_withcharge.delta_charge = delta_charge
                s_withcharge.end_time = end_time
                s_withcharge.end_charge = end_charge
                s_withcharge.round_time = round_time
                s_withcharge.round_charge = round_charge
                push!(nondominated_subpaths_withcharge, s_withcharge)
            end
        end 
    end
    return nondominated_subpaths_withcharge
end

function enumerate_all_subpaths(
    G, 
    data, 
    T_range, 
    B_range,
    ;
    charging_in_subpath::Bool = false,
)
    """
    `all_subpaths`: Dictionary mapping (I,J)-pairs of states
    to a Vector{Subpath}.
    """
    start_time = time()
    all_subpaths = Dict()
    for (starting_node, starting_time, starting_charge) in Iterators.flatten((
        Iterators.product(
            data["N_charging"],
            T_range,
            B_range,
        ),
        Iterators.product(
            data["N_depots"],
            [0.0],
            [data["B"]],
        ),
    ))
        if charging_in_subpath
            subpaths = @suppress enumerate_subpaths_withcharge(
                starting_node, starting_time, starting_charge,
                G, data, T_range, B_range,
            )
        else
            _, subpaths = @suppress enumerate_subpaths(
                starting_node, starting_time, starting_charge,
                G, data,
            )
        end
        for s in subpaths
            # perform rounding here
            s.round_time = dceil(s.end_time, T_range)
            s.round_charge = dfloor(s.end_charge, B_range)
            key = (
                (starting_node, starting_time, starting_charge),
                (s.current_node, s.round_time, s.round_charge)
            )
            if !(key in keys(all_subpaths))
                all_subpaths[key] = []
            end
            push!(all_subpaths[key], s)
        end
    end
    end_time = time()
    return all_subpaths, end_time - start_time
    return all_subpaths
end

function compute_subpath_costs(
    data,
    all_subpaths,
)
    subpath_costs = Dict(
        key => [
            sum(data["c"][a...] for a in s.arcs)
            for s in all_subpaths[key]
        ]
        for key in keys(all_subpaths)
    )
    return subpath_costs
end

function compute_subpath_service(
    data, 
    all_subpaths,
)
    subpath_service = Dict(
        (key, i) => [
            s.served[i]
            for s in all_subpaths[key]
        ]
        for key in keys(all_subpaths), i in 1:data["n_customers"]
    )
    return subpath_service
end

function enumerate_all_charging_arcs(
    data,
    all_subpaths,
    T_range,
    B_range,
)
    states = Set()
    for k in keys(all_subpaths)
        push!(states, k[1], k[2])
    end
    all_charging_arcs = []
    for starting_node in data["N_charging"]
        for t1 in T_range
            for b1 in B_range
                if !((starting_node, t1, b1) in states)
                    continue 
                end
                append!(
                    all_charging_arcs, 
                    [
                        ((starting_node, t1, b1), (starting_node, t2, b2))
                        for (_, _, _, _, t2, b2) in generate_charging_options(
                            t1,
                            b1,
                            data,
                            T_range,
                            B_range,
                        )
                        if (starting_node, t2, b2) in states
                    ]
                )
            end
        end
    end
    charging_arcs_costs = Dict(state_pair => 0.0 for state_pair in all_charging_arcs)
    return all_charging_arcs, charging_arcs_costs
end

function subpath_formulation(
    data, 
    all_subpaths,
    subpath_costs, 
    subpath_service,
    T_range,
    B_range,
    ;
    charging_in_subpath::Bool = false,
    all_charging_arcs = [],
    charging_arcs_costs = Dict(),
)
    # Charging always present

    start_time = time()

    model = Model(Gurobi.Optimizer)

    states = Set()
    for k in keys(all_subpaths)
        push!(states, k[1], k[2])
    end

    m = maximum(length(x) for x in values(all_subpaths))
    @variable(model, z[
        k=keys(all_subpaths), 
        p=1:m; p ≤ length(all_subpaths[k])
    ], Bin);
    if !charging_in_subpath
        @variable(model, y[all_charging_arcs], Bin)
    end

    # Constraint (4b): number of subpaths starting at depot i
    # at time 0 with full charge to be v_startᵢ
    @constraint(
        model,
        [i in data["N_depots"]],
        sum(
            sum(
                z[((i,0,data["B"]),state),p]
                for p in 1:length(all_subpaths[((i,0,data["B"]),state)])
            )        
            for state in states
                if ((i,0,data["B"]),state) in keys(all_subpaths)
        )
        == data["v_start"][findfirst(x -> (x == i), data["N_depots"])]
    )

    flow_conservation_exprs = Dict()
    for n1 in data["N_charging"]
        for t1 in T_range
            t1_floor = dfloorall(t1, T_range)
            t1_ceil = dceilall(t1, T_range)
            for b1 in B_range
                state1 = (n1, t1, b1)
                b1_floor = dfloorall(b1, B_range)
                b1_ceil = dceilall(b1, B_range)
                out_sp_neighbors = [state for state in states if (state1, state) in keys(all_subpaths)]
                in_sp_neighbors = [state for state in states if (state, state1) in keys(all_subpaths)]
                flow_conservation_exprs[(state1,"out_sp")] = @expression(
                    model, 
                    sum(
                        sum(
                            z[(state1, state2),p] 
                            for p in 1:length(all_subpaths[(state1, state2)])
                        )
                        for state2 in out_sp_neighbors
                    )
                )
                flow_conservation_exprs[(state1,"in_sp")] = @expression(
                    model, 
                    sum(
                        sum(
                            z[(state2, state1),p] 
                            for p in 1:length(all_subpaths[(state2, state1)])
                        )
                        for state2 in in_sp_neighbors
                    )
                )
                if charging_in_subpath
                    # Constraint (4c) Flow conservation at charging stations
                    # (here the other node can be either a depot or charging station)
                    @constraint(
                        model,
                        flow_conservation_exprs[(state1,"out_sp")]
                        == flow_conservation_exprs[(state1,"in_sp")]
                    )
                else
                    out_a_neighbors = [(n1, t2, b2) for t2 in t1_ceil, b2 in b1_ceil if (state1, (n1, t2, b2)) in all_charging_arcs]
                    in_a_neighbors = [(n1, t2, b2) for t2 in t1_floor, b2 in b1_floor if ((n1, t2, b2), state1) in all_charging_arcs]
                    flow_conservation_exprs[(state1,"out_a")] = @expression(
                        model, 
                        sum(y[(state1, state2)] for state2 in out_a_neighbors)
                    )
                    flow_conservation_exprs[(state1,"in_a")] = @expression(
                        model, 
                        sum(y[(state2, state1)] for state2 in in_a_neighbors)
                    )
                    # FIXME: make more efficient
                    # Constraint (7b) Flow conservation at charging stations
                    # (here the other node can be either a depot or charging station)
                    # (here the arcs can be subpaths or charging arcs)
                    @constraint(
                        model,
                        flow_conservation_exprs[(state1,"out_sp")]
                        + flow_conservation_exprs[(state1,"out_a")]
                        == flow_conservation_exprs[(state1,"in_sp")]
                        + flow_conservation_exprs[(state1,"in_a")]
                    )
                    # FIXME: make more efficient
                    # Constraint (7c): having at most 1 charging arc incident to any charging node
                    @constraint(
                        model,
                        flow_conservation_exprs[(state1,"out_a")]
                        + flow_conservation_exprs[(state1,"in_a")]
                        ≤ 1
                    )
                end
            end
        end
    end

    # Constraint (4d) Serving all customers exactly once across all subpaths
    @constraint(
        model,
        [j in data["N_pickups"]],
        sum(
            sum(
                subpath_service[((state1, state2),j)][p] * z[(state1, state2),p]
                for p in 1:length(all_subpaths[(state1, state2)])
            )
            for (state1, state2) in keys(all_subpaths)
        ) == 1
    )

    # Constraint (4e) Number of subpaths ending at each depot n2 
    # is at least the number of vehicles required
    @constraint(
        model,
        [n2 in data["N_depots"]],
        sum(
            sum(
                z[((n1,t1,b1),(n2,t2,b2)),p]
                for p in 1:length(all_subpaths[((n1,t1,b1),(n2,t2,b2))])
            )
            for (n1, t1, b1) in states,
                t2 in dceilall(t1, T_range), 
                b2 in dfloorall(b1, B_range)
            if ((n1,t1,b1),(n2,t2,b2)) in keys(all_subpaths)
        ) ≥ data["v_end"][n2]
    )

    @expression(
        model,
        subpath_costs_expr,
        sum(
            sum(
                subpath_costs[state_pair][p] * z[state_pair,p]
                for p in 1:length(all_subpaths[state_pair])
            )
            for state_pair in keys(all_subpaths)
        )
    )
    if charging_in_subpath
        # Objective: minimize the cost of all chosen subpaths, 
        # and the cost of all charging arcs
        @objective(model, Min, subpath_costs_expr)
    else
        @expression(
            model, 
            charging_arcs_costs_expr,
            sum(
                y[state_pair] * charging_arcs_costs[state_pair] 
                for state_pair in all_charging_arcs
            )
        )
        # Objective: minimize the cost of all chosen subpaths, 
        # and the cost of all charging arcs
        @objective(model, Min, subpath_costs_expr + charging_arcs_costs_expr)
    end

    constraint_end_time = time()

    optimize!(model)
    end_time = time()
    time_taken = end_time - start_time
    constraint_time_taken = constraint_end_time - start_time
    solution_time_taken = end_time - constraint_end_time

    params = Dict(
        "time_taken" => time_taken,
        "constraint_time_taken" => constraint_time_taken,
        "solution_time_taken" => solution_time_taken,
    )
    results = Dict(
        "objective" => objective_value(model),
        "z" => value.(z),
    )
    if !charging_in_subpath
        results["y"] = value.(y)
    end
    return results, params
end

function construct_paths_from_subpath_solution(
    results, data, all_subpaths, 
    ;
    charging_in_subpath::Bool = false,
    all_charging_arcs::Vector{Any} = [],
)
    results_subpaths = []
    for key in keys(all_subpaths)
        for p in 1:length(all_subpaths[key])
            if results["z"][key,p] != 0
                push!(results_subpaths, all_subpaths[key][p])
            end
        end
    end
    if !charging_in_subpath
        results_charging_arcs = []
        for key in all_charging_arcs
            if results["y"][key] != 0
                push!(results_charging_arcs, key)
            end
        end
    end

    paths = Dict(v => [] for v in 1:data["n_vehicles"])
    schedules = Dict(v => [] for v in 1:data["n_vehicles"])
    for v in 1:data["n_vehicles"]
        t = 0
        current_state = nothing
        while true
            if t == 0
                i = findfirst(
                    s -> (
                        s.starting_node in data["N_depots"]
                        && s.starting_time == 0.0
                        && s.starting_charge == data["B"]
                    ), 
                    results_subpaths
                )
            else
                i = findfirst(
                    s -> (
                        (
                            s.starting_node, 
                            s.starting_time,
                            s.starting_charge
                        ) == current_state
                    ),
                    results_subpaths
                )
            end
            if !charging_in_subpath && isnothing(i) && t != 0
                i = findfirst(
                    s -> (s[1] == current_state),
                    results_charging_arcs
                )
                current_a = popat!(results_charging_arcs, i)
                push!(schedules[v], current_a)
                push!(paths[v], (current_a[1][1], current_a[2][1]))
                current_state = current_a[2]
            else
                current_s = popat!(results_subpaths, i)
                push!(schedules[v], current_s)
                append!(paths[v], current_s.arcs)
                current_state = (current_s.current_node, current_s.round_time, current_s.round_charge)
            end
            if current_state[1] in data["N_depots"]
                break
            end        
            t += 1
        end
    end
    return paths, schedules
end

function subpath_results_printout(
    results,
    params,
    data,
    all_subpaths,
    ;
    charging_in_subpath::Bool = true,
    all_charging_arcs::Vector{Any} = [],
)

    @printf("Objective:                 %7.1f\n", results["objective"])
    println("")
    @printf("Time taken:                %7.1f s\n", params["time_taken"])
    @printf("Time taken (formulation):  %7.1f s\n", params["constraint_time_taken"])
    @printf("Time taken (solving):      %7.1f s\n", params["solution_time_taken"])
    println("")

    println("Vehicles         From      time   charge             To      time   charge  |   sp_time |   sp_charge ")
    println("----------------------------------------------------------------------------|-----------|-------------")

    paths, schedules = construct_paths_from_subpath_solution(
        results, data, all_subpaths;
        charging_in_subpath = charging_in_subpath,
        all_charging_arcs = all_charging_arcs,
    )
    
    for k in data["N_vehicles"]
        for s in schedules[k]
            if isa(s, Subpath)
                nodes = []
                times = []
                charges = []
                push!(nodes, data["node_labels"][s.starting_node])
                current_time = s.starting_time
                push!(times, current_time)
                current_charge = s.starting_charge   
                push!(charges, current_charge)
                for (i, a) in enumerate(s.arcs)
                    push!(nodes, data["node_labels"][a[2]])
                    current_time = max(current_time + data["t"][a...], data["α"][a[2]])
                    push!(times, current_time)
                    current_charge = current_charge - data["q"][a...]
                    push!(charges, current_charge)
                end
                @printf(
                    "Vehicle %s:                                                                  |    %6.1f |      %6.1f \n", 
                    k, 
                    s.starting_time, 
                    s.starting_charge,
                )
                for i in 1:length(nodes)-2
                    @printf(
                        "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) |           | \n", 
                        k, 
                        nodes[i], 
                        times[i], 
                        charges[i],
                        nodes[i+1],
                        times[i+1],
                        charges[i+1],
                    )
                end
                @printf(
                    "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | +  %6.1f | -    %6.1f \n", 
                    k, 
                    nodes[length(nodes)-1], 
                    times[length(nodes)-1], 
                    charges[length(nodes)-1],
                    nodes[length(nodes)],
                    times[length(nodes)],
                    charges[length(nodes)],
                    s.time - s.starting_time, 
                    s.starting_charge - s.charge,
                )
                if charging_in_subpath && s.current_node in data["N_charging"]
                    @printf(
                        "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | +  %6.1f | +    %6.1f \n", 
                        k, 
                        data["node_labels"][s.current_node],
                        s.time,
                        s.charge,
                        data["node_labels"][s.current_node],
                        s.end_time,
                        s.end_charge,
                        s.delta_time, 
                        s.delta_charge,
                    )
                end
                # Rounding
                @printf(
                    "Vehicle %s:                                -> %12s (%6.1f,  %6.1f) | +  %6.1f | -    %6.1f \n", 
                    k, 
                    data["node_labels"][s.current_node],
                    s.round_time,
                    s.round_charge,
                    s.round_time - s.end_time,
                    s.end_charge - s.round_charge,
                )
            else
                # Charging arc 
                @printf(
                    "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | +  %6.1f | +    %6.1f \n", 
                    k, 
                    data["node_labels"][s[1][1]],
                    s[1][2],
                    s[1][3],
                    data["node_labels"][s[2][1]],
                    s[2][2],
                    s[2][3],
                    s[2][2] - s[1][2],
                    s[2][3] - s[1][3],
                )
            end
        end
        println("------------------------------------------------------------------------------------------------------")
    end
end