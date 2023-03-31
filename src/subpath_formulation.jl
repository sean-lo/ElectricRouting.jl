using Suppressor
using Printf
using DataStructures

using JuMP
using Gurobi

include("utils.jl")

struct TripleStateOrdering <: Base.Order.Ordering
end
struct DoubleStateOrdering <: Base.Order.Ordering
end

import Base.Order.lt
import DataStructures.eq
lt(o::TripleStateOrdering, a, b) = isless((a[2], -a[3], a[1]), (b[2], -b[3], b[1]))
lt(o::DoubleStateOrdering, a, b) = isless((a[1], -a[2]), (b[1], -b[2]))
eq(o::TripleStateOrdering, a, b) = isequal(a, b)
eq(o::DoubleStateOrdering, a, b) = isequal(a, b)

function construct_heuristic_paths(
    data, 
    T_range, 
    B_range,
)
    groups = sort(
        Dict(
            [i] => data["c"][i, i + data["n_customers"]]
            for i in data["N_pickups"]
        ),
        byvalue = true,
    )
    arclist_cc = sort(
        Dict(
            (k1, k2) => (
                groups[k1]
                + data["c"][k1[end] + data["n_customers"], k2[1]]
                + groups[k2]
            )
            for (k1, k2) in collect(permutations(collect(keys(groups)), 2))
        ),
        byvalue = true
    )
    nearest_charge = merge(
        Dict(
            (i, "in") => minimum(
                data["c"][c,i] for c in data["N_charging"]
            )
            for i in data["N_pickups"]
        ),
        Dict(
            (i - data["n_customers"], "out") => minimum(
                data["c"][i,c] for c in data["N_charging"]
            )
            for i in data["N_dropoffs"]
        ),
    )

    while true
        myk1, myk2, myval = nothing, nothing, nothing
        for ((k1, k2), val) in pairs(arclist_cc)
            safety_cost = (
                val 
                + nearest_charge[(k1[1], "in")] 
                + nearest_charge[(k2[end], "out")]
            )
            if safety_cost ≤ data["B"]
                myk1, myk2, myval = k1, k2, val
                break
            end
        end
        if isnothing(myval)
            break
        end
    
        pop!(groups, myk1)
        pop!(groups, myk2)
        groups[vcat(myk1,myk2)] = myval
        arclist_cc = sort(
            Dict(
                (k1, k2) => (
                    groups[k1]
                    + data["c"][k1[end] + data["n_customers"], k2[1]]
                    + groups[k2]
                )
                for (k1, k2) in collect(permutations(collect(keys(groups)), 2))
            ),
            byvalue = true
        )
    end

    base_labels = Dict(
        start_node => Dict(
            end_node => Dict{Tuple{Float64, Float64}, Vector{Subpath}}()
            for end_node in union(data["N_depots"], data["N_charging"])
        )
        for start_node in union(data["N_depots"], data["N_charging"])
    )
    for k in keys(groups)
        for start_node in union(data["N_depots"], data["N_charging"])
            for end_node in union(data["N_depots"], data["N_charging"])
                nodelist = vcat(
                    start_node,
                    reduce(vcat, [[i, i+data["n_customers"]] for i in k]),
                    end_node,
                )
                arcs = collect(zip(nodelist[1:end-1], nodelist[2:end]))
                
                starting_time = 0.0
                starting_charge = data["B"]
                time_taken = sum(data["c"][a...] for a in arcs)
                charge_taken = sum(data["q"][a...] for a in arcs)
                
                current_time = starting_time + time_taken
                current_charge = starting_charge - charge_taken
                if !(current_time ≤ data["T"])
                    continue
                end
                if !(0.0 ≤ current_charge ≤ data["B"])
                    continue
                end
                served = falses(data["n_customers"])
                served[k] .= true
                for i in 1:(length(arcs) ÷ 2)
                    customer = nodelist[2*i]
                    cumulative_time_taken = sum(data["c"][a...] for a in arcs[1:2*i])
                end
                if end_node in data["N_charging"]
                    (dt, db, et, eb, rt, rb) = generate_charging_options(
                        starting_time + time_taken, 
                        starting_charge - charge_taken, 
                        data["μ"],
                        T_range, B_range,
                        charge_to_full_only = true,
                    )[1]
                    s = Subpath(
                        n_customers = data["n_customers"],
                        starting_node = start_node,
                        starting_time = starting_time,
                        starting_charge = starting_charge,
                        current_node = end_node,
                        arcs = arcs,
                        time = current_time,
                        charge = current_charge,
                        served = served,
                        delta_time = dt, 
                        delta_charge = db,
                        end_time = et, 
                        end_charge = eb,
                        round_time = rt, 
                        round_charge = rb,
                    )
                    if !((rt, rb) in keys(base_labels[start_node][end_node]))
                        base_labels[start_node][end_node][(rt, rb)] = []
                    end
                    push!(base_labels[start_node][end_node][(rt, rb)], s)
                else
                    s = Subpath(
                        n_customers = data["n_customers"],
                        starting_node = start_node,
                        starting_time = starting_time,
                        starting_charge = starting_charge,
                        current_node = end_node,
                        arcs = arcs,
                        time = current_time,
                        charge = current_charge,
                        served = served,
                        delta_time = 0.0, 
                        delta_charge = 0.0,
                        end_time = current_time, 
                        end_charge = current_charge,
                        round_time = current_time, 
                        round_charge = current_charge,
                    )
                    if !((current_time, current_charge) in keys(base_labels[start_node][end_node]))
                        base_labels[start_node][end_node][(current_time, current_charge)] = []
                    end
                    push!(base_labels[start_node][end_node][(current_time, current_charge)], s)
                end
            end
        end
    end

    unexplored_states = SortedSet(
        TripleStateOrdering(), 
        [
            (depot, 0.0, data["B"])
            for depot in data["N_depots"]
        ]
    )

    full_labels = Dict(
        start_node => Dict(
            end_node => Dict{Tuple{Float64, Float64}, Vector{Path}}()
            for end_node in union(data["N_depots"], data["N_charging"])
        )
        for start_node in data["N_depots"]
    )
    for depot in data["N_depots"]
        full_labels[depot][depot][(0.0, data["B"])] = [Path(
            subpaths = [], served = zeros(Int, data["n_customers"])
        )]
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        # where do you want to go next
        for next_node in union(data["N_depots"], data["N_charging"])
            # what subpaths bring you there
            for s_list in values(base_labels[state[1]][next_node])
                for s in s_list
                    if !(
                        s.round_time + state[2] ≤ data["T"]
                    )
                        continue
                    end
                    # translate subpath
                    s_new = copy(s)
                    s_new.starting_time = state[2]
                    s_new.starting_charge = data["B"]
                    s_new.time = s.time + state[2]
                    s_new.charge = s.charge
                    s_new.end_time = s.end_time + state[2]
                    s_new.end_charge = s.end_charge
                    s_new.round_time = s.round_time + state[2]
                    s_new.round_charge = s.round_charge
                    if next_node in data["N_depots"]
                        key = (s_new.end_time, s_new.end_charge)
                    else
                        key = (s_new.round_time, s_new.round_charge)
                    end
                    add_next_state = false
                    for starting_node in data["N_depots"]
                        if !((state[2], state[3]) in keys(full_labels[starting_node][state[1]]))
                            continue
                        end
                        for current_path in full_labels[starting_node][state[1]][(state[2], state[3])]
                            new_path = Path(
                                subpaths = vcat(
                                    current_path.subpaths, 
                                    [s_new],
                                ),
                            )
                            add_next_state = true
                            if !(key in keys(full_labels[starting_node][next_node]))
                                full_labels[starting_node][next_node][key] = []
                            end
                            push!(full_labels[starting_node][next_node][key], new_path)
                        end
                    end
                    next_state = (next_node, key...)
                    if (
                        add_next_state  
                        && next_node in data["N_charging"] 
                        && !(next_state in unexplored_states)
                    )
                        insert!(unexplored_states, next_state)
                    end
                end
            end
        end
    end

    for depot in data["N_depots"]
        full_labels[depot][depot][(0.0, data["B"])] = [Path(
            subpaths = [Subpath(
                n_customers = data["n_customers"],
                starting_node = depot,
                starting_time = 0.0,
                starting_charge = data["B"],
                arcs = [(depot, depot)],
            )], 
            served = zeros(Int, data["n_customers"])
        )]
    end

    heuristic_paths = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Path},
    }()
    for start_node in data["N_depots"], end_node in data["N_depots"]
        for ((end_time, end_charge), p_list) in pairs(full_labels[start_node][end_node])
            for p in p_list
                start_state = (start_node, 0.0, data["B"])
                end_state = (end_node, end_time, end_charge)
                state_pair = (start_state, end_state)
                if !(state_pair in keys(heuristic_paths))
                    heuristic_paths[state_pair] = []
                end
                if !any(isequal(p, p1) for p1 in heuristic_paths[state_pair])
                    push!(heuristic_paths[state_pair], p)
                end
            end
        end
    end

    return heuristic_paths

end

function construct_heuristic_subpaths(
    data, 
    T_range, 
    B_range,
)
    
    heuristic_paths = construct_heuristic_paths(data, T_range, B_range)

    heuristic_subpaths = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Subpath},
    }()

    for path_list in values(heuristic_paths)
        for path in path_list
            for s in path.subpaths
                state_pair = (
                    (s.starting_node, s.starting_time, s.starting_charge),
                    (s.current_node, s.round_time, s.round_charge),
                )
                if !(state_pair in keys(heuristic_subpaths))
                    heuristic_subpaths[state_pair] = []
                end
                push!(heuristic_subpaths[state_pair], s)
            end
        end
    end

    return heuristic_subpaths

end

function enumerate_subpaths(
    starting_node, 
    starting_time, 
    starting_charge,
    G,
    data,
    ;
    time_windows::Bool = true,
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
                if time_windows
                    current_time = max(
                        s.time + data["t"][s.current_node,j],
                        data["α"][j],
                    )
                    feasible_timewindow = (current_time ≤ data["β"][j])
                    current_time = max(
                        current_time + data["t"][j,new_node],
                        data["α"][new_node],
                    )
                    feasible_timewindow = feasible_timewindow && (current_time ≤ data["β"][new_node])
                else
                    current_time = s.time + data["t"][s.current_node,j] + data["t"][j,new_node]
                    feasible_timewindow = (current_time ≤ data["T"])
                end
                new_charge = s.charge - data["q"][s.current_node,j] - data["q"][j,new_node]
                feasible_charge = (new_charge ≥ qmin[new_node])
            else
                new_node = j
                if time_windows
                    current_time = max(
                        s.time + data["t"][s.current_node,new_node],
                        data["α"][new_node],
                    )
                    feasible_timewindow = (current_time ≤ data["β"][new_node])
                else
                    current_time = s.time + data["t"][s.current_node,new_node]
                    feasible_timewindow = (current_time ≤ data["T"])
                end
                new_charge = s.charge - data["q"][s.current_node,new_node]
                feasible_charge = (new_charge ≥ qmin[new_node])
            end

            feasible = (
                feasible_timewindow
                && feasible_charge
            )
            if !feasible
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
            s_j.time = current_time
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

function generate_charging_options(
    starting_time,
    starting_charge,
    charging_rate,
    T_range,
    B_range,
    ;
    charge_to_full_only::Bool = false,
)
    # Version 2: charge by fixed time intervals inside T_range
    # If maximum time reached, charge to maximum time and get the corresponding charge
    # If maximum charge reached, charge to the next higher time (in T_range) and get to max charge
    max_delta_time = min(
        (B_range[end] - starting_charge) / charging_rate,
        T_range[end] - starting_time,
    )
    max_end_time = starting_time + max_delta_time
    max_end_time = dceil(max_end_time, T_range)
    if charge_to_full_only
        end_times = [max_end_time]
    else
        end_times = [
            t for t in T_range
            if starting_time < t ≤ max_end_time
        ]
    end
    round_times = end_times
    delta_times = end_times .- starting_time
    delta_charges = delta_times .* charging_rate
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
    ;
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
)
    """
    `nondominated_subpaths_withcharge`: Vector{Subpath}, 
    generated from a single Subpath with additional discretized charging options
    """
    _, nondominated_subpaths = enumerate_subpaths(
        starting_node, starting_time, starting_charge, G, data;
        time_windows = time_windows,
    )
    
    nondominated_subpaths_withcharge = []
    for s in nondominated_subpaths
        if s.current_node in data["N_depots"]
            s.round_time = dceil(s.end_time, T_range)
            s.round_charge = dfloor(s.end_charge, T_range)
            push!(nondominated_subpaths_withcharge, s)
        else
            append!(nondominated_subpaths_withcharge,
                [
                    Subpath(
                        n_customers = data["n_customers"],
                        starting_node = s.starting_node,
                        starting_time = s.starting_time,
                        starting_charge = s.starting_charge,
                        current_node = s.current_node,
                        arcs = copy(s.arcs),
                        time = s.time,
                        charge = s.charge,
                        served = copy(s.served),
                        delta_time = delta_time,
                        delta_charge = delta_charge,
                        end_time = end_time,
                        end_charge = end_charge,
                        round_time = round_time,
                        round_charge = round_charge,
                    )
                    for (delta_time, delta_charge, 
                        end_time, end_charge, 
                        round_time, round_charge) in generate_charging_options(
                        s.time, s.charge, data["μ"], T_range, B_range,
                        ;
                        
                        charge_to_full_only = charge_to_full_only,
                    )
                ]
            )
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
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
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
                ;
                
                charge_to_full_only = charge_to_full_only,
                time_windows = time_windows,
            )
        else
            _, subpaths = @suppress enumerate_subpaths(
                starting_node, starting_time, starting_charge,
                G, data,
                ;
                time_windows = time_windows,
            )
            for s in subpaths
                # perform rounding here
                s.round_time = dceil(s.end_time, T_range)
                s.round_charge = dfloor(s.end_charge, B_range)
            end
        end
        for s in subpaths
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
    return all_subpaths, round(end_time - start_time, digits=3)
end

function enumerate_all_subpaths_faster(
    G, 
    data, 
    T_range, 
    B_range,
    ;
    charging_in_subpath::Bool = false,
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
)
    """
    `all_subpaths`: Dictionary mapping (I,J)-pairs of states
    to a Vector{Subpath}.
    """
    start_time = time()
    all_subpaths = Dict()
    unexplored_states = Set(
        Iterators.product(
            data["N_depots"],
            [0.0],
            [data["B"]],
        )
    )
    explored_states = Set()
    while length(unexplored_states) > 0
        next_unexplored_states = Set()
        for state in unexplored_states
            (starting_node, starting_time, starting_charge) = state
            if charging_in_subpath
                subpaths = @suppress enumerate_subpaths_withcharge(
                    starting_node, starting_time, starting_charge,
                    G, data, T_range, B_range,
                    ;
                    
                    charge_to_full_only = charge_to_full_only,
                    time_windows = time_windows,
                )
            else
                _, subpaths = @suppress enumerate_subpaths(
                    starting_node, starting_time, starting_charge,
                    G, data,
                    ;
                    time_windows = time_windows,
                )
                for s in subpaths
                    # perform rounding here
                    s.round_time = dceil(s.end_time, T_range)
                    s.round_charge = dfloor(s.end_charge, B_range)
                end
            end
            for s in subpaths
                next_state = (s.current_node, s.round_time, s.round_charge)
                key = (state, next_state)
                if !(key in keys(all_subpaths))
                    all_subpaths[key] = []
                end
                push!(all_subpaths[key], s)
                push!(next_unexplored_states, next_state)
            end
            push!(explored_states, state)
        end
        unexplored_states = setdiff(next_unexplored_states, explored_states)
    end
    end_time = time()
    return all_subpaths, end_time - start_time
end

function generate_artificial_subpaths(data)
    artificial_subpaths = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}},
        Vector{Subpath},
    }()
    start_depots = zeros(Int, data["n_vehicles"])
    for (k, v_list) in pairs(data["V"])
        for v in v_list
            start_depots[v] = k
        end
    end
    end_depots = []
    for k in data["N_depots"]
        append!(end_depots, repeat([k], data["v_end"][k]))
    end
    append!(end_depots, 
        repeat(
            [data["N_depots"][1]], 
            outer = data["n_vehicles"] - sum(values(data["v_end"]))
        )
    )
    for (v, (starting_node, ending_node)) in enumerate(zip(start_depots, end_depots))
        starting_time = 0.0
        starting_charge = data["B"]
        ending_time = 0.0
        ending_charge = data["B"]
        key = (
            (starting_node, starting_time, starting_charge),  
            (ending_node, ending_time, ending_charge)
        )
        # initialise a proportion of the customers to be served
        served = falses(data["n_customers"])
        for i in 1:length(served)
            if mod1(i, data["n_vehicles"]) == v
                served[i] = true
            end
        end
        s = Subpath(
            n_customers = data["n_customers"],
            starting_node = starting_node,
            starting_time = starting_time,
            starting_charge = starting_charge,
            current_node = ending_node,
            arcs = [],
            time = ending_time,
            charge = ending_charge,
            served = served,
            artificial = true,
        )
        if !(key in keys(artificial_subpaths))
            artificial_subpaths[key] = []
        end
        push!(artificial_subpaths[key], s)
    end
    return artificial_subpaths
end

function compute_subpath_reduced_cost(
    s::Subpath,
    data,
    κ,
    λ,
    μ,
    ν,
    ; 
    with_lambda::Bool = true,
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
    verbose::Bool = false,
)
    reduced_cost = compute_subpath_cost(
        data, s,
        ;
        with_charging_cost = with_charging_cost,
        time_windows = time_windows,
        verbose = verbose,
    )

    for j in findall(s.served)
        reduced_cost = reduced_cost - ν[j]
    end

    if s.starting_node in data["N_depots"]
        if s.starting_time == 0.0 && s.starting_charge == data["B"]
            reduced_cost = reduced_cost - κ[s.starting_node]
        end
    elseif s.starting_node in data["N_charging"] && with_lambda
        reduced_cost = reduced_cost - λ[(s.starting_node, s.starting_time, s.starting_charge)]
    end

    if s.current_node in data["N_depots"]
        reduced_cost = reduced_cost - μ[s.current_node]
    elseif s.current_node in data["N_charging"] && with_lambda
        reduced_cost = reduced_cost + λ[(s.current_node, s.round_time, s.round_charge)]
    end

    return reduced_cost
end

function find_smallest_reduced_cost_paths(
    starting_node,
    G, 
    data, 
    T_range, 
    B_range,
    # dual values associated with flow conservation constraints
    κ,
    μ,
    # dual values associated with customer service constraints
    ν,
    ;
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
    verbose::Bool = false,
)

    # initialize modified arc costs, subtracting values of dual variables
    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for i in 1:data["n_customers"]
        j = data["n_customers"] + i
        modified_costs[i,j] -= ν[i]
        modified_costs[j,i] -= ν[i]
    end
    # initialize set of labels
    initial_cost = - κ[starting_node]
    starting_time = 0.0
    starting_charge = data["B"]
    labels = Dict()
    labels[starting_node] = Dict(
        (starting_time, starting_charge) => PathWithCost(
            subpaths = [
                SubpathWithCost(
                    cost = initial_cost,
                    n_customers = data["n_customers"],
                    starting_node = starting_node,
                    starting_time = starting_time,
                    starting_charge = starting_charge,
                )
            ],
            cost = initial_cost,
        )
    )
    Q = [starting_node]

    while !isempty(Q)
        i = popfirst!(Q)
        for ((now_time, now_charge), path) in pairs(labels[i])
            if path.explored
                continue
            else 
                path.explored = true
            end
            for j in setdiff(outneighbors(G, i), i)
                add_to_queue = false
                if j in data["N_pickups"]
                    # feasible_service = !any(s.served[j] for s in path.subpaths)
                    feasible_service = !path.subpaths[end].served[j]
                else
                    feasible_service = true
                end
                if !feasible_service 
                    continue
                end

                if j in data["N_pickups"]
                    current_node = j + data["n_customers"]
                    if time_windows
                        current_time = max(
                            data["α"][j], 
                            now_time + data["t"][i,j],
                        )
                        feasible_timewindow = (current_time ≤ data["β"][j])
                        current_time = max(
                            data["α"][current_node], 
                            current_time + data["t"][j, current_node],
                        )
                        feasible_timewindow = feasible_timewindow && (current_time ≤ data["β"][current_node])
                    else
                        current_time = now_time + data["t"][i,j] + data["t"][j, current_node]
                        feasible_timewindow = (current_time ≤ data["T"])
                    end
                    current_charge = now_charge - data["q"][i,j]
                    feasible_charge = (current_charge ≥ 0.0)
                    current_charge = current_charge - data["q"][j, current_node]
                    feasible_charge = feasible_charge && (current_charge ≥ 0.0)
                else
                    current_node = j
                    if time_windows
                        current_time = max(
                            data["α"][j], 
                            now_time + data["t"][i,j],
                        )
                        feasible_timewindow = (current_time ≤ data["β"][j])
                    else
                        current_time = now_time + data["t"][i,j]
                        feasible_timewindow = (current_time ≤ data["T"])
                    end
                    current_charge = now_charge - data["q"][i,j]
                    feasible_charge = (current_charge ≥ 0.0)
                end

                feasible = (
                    feasible_timewindow
                    && feasible_charge
                )
                if !feasible
                    continue
                end
                
                current_cost = path.subpaths[end].cost
                cumulative_cost = path.cost
                current_cost = current_cost + modified_costs[i,j]
                cumulative_cost = cumulative_cost + modified_costs[i,j]
                if j in data["N_pickups"]
                    current_cost = current_cost + modified_costs[j, current_node]
                    cumulative_cost = cumulative_cost + modified_costs[j, current_node]
                elseif j in data["N_depots"]
                    current_cost = current_cost - μ[j]
                    cumulative_cost = cumulative_cost - μ[j]
                end

                if j in data["N_charging"]
                    charging_options = generate_charging_options(
                        current_time, current_charge, 
                        data["μ"], T_range, B_range,
                        ;
                        
                        charge_to_full_only = charge_to_full_only,    
                    )
                    if length(charging_options) == 0
                        continue
                    end
                else
                    charging_options = [(
                        0, 0, current_time, current_charge, 
                        current_time, current_charge
                    )]
                end
    
                
                for (delta_time, delta_charge, end_time, end_charge, store_time, store_charge) in charging_options
                    # determine if label ought to be updated
                    this_current_cost = current_cost
                    this_cumulative_cost = cumulative_cost 
                    if with_charging_cost
                        this_current_cost += data["charge_cost_coeff"] * delta_time
                        this_cumulative_cost += data["charge_cost_coeff"] * delta_time
                    end

                    if current_node in keys(labels)
                        if (store_time, store_charge) in keys(labels[current_node])
                            if labels[current_node][(store_time, store_charge)].cost ≤ this_cumulative_cost
                                if verbose
                                    println("$i, $now_time, $now_charge dominated by $current_node, $store_time, $store_charge")
                                end
                                continue
                            end
                        end
                    end
    
                    # update labels
                    add_to_queue = true
                    s_j = copy(path.subpaths[end])
                    s_j.cost = this_current_cost
                    s_j.current_node = current_node
                    if j in data["N_pickups"]
                        push!(s_j.arcs, (i, j))
                        push!(s_j.arcs, (j, current_node))
                    else
                        push!(s_j.arcs, (i, current_node))
                    end
                    s_j.time = current_time
                    s_j.charge = current_charge
                    if current_node in data["N_dropoffs"]
                        s_j.served[j] = true
                    end
                    s_j.delta_time = delta_time
                    s_j.delta_charge = delta_charge
                    s_j.end_time = end_time
                    s_j.end_charge = end_charge
                    s_j.round_time = dceil(end_time, T_range)
                    s_j.round_charge = dfloor(end_charge, B_range)
    
                    s_list_new = vcat(path.subpaths[1:end-1], [s_j])
                    if current_node in data["N_charging"]
                        push!(
                            s_list_new, 
                            # new null subpath starting at current_node
                            SubpathWithCost(
                                cost = 0,
                                n_customers = data["n_customers"],
                                starting_node = current_node,
                                starting_time = store_time,
                                starting_charge = store_charge,
                            )
                        )
                    end
                    if !(current_node in keys(labels))
                        labels[current_node] = Dict{Tuple{Float64, Float64}, PathWithCost}()
                    end

                    labels[current_node][(store_time, store_charge)] = PathWithCost(
                        subpaths = s_list_new,
                        cost = this_cumulative_cost,
                    )
                    if verbose
                        println("Added $current_node, $store_time, $store_charge from $i, $now_time, $now_charge")
                    end
                end
                if add_to_queue && !(current_node in Q) && !(current_node in data["N_depots"])
                    push!(Q, current_node)
                end
            end
        end
    end
    # remove labels corresponding to pickup, dropoff, and charging stations
    for k in union(data["N_pickups"], data["N_dropoffs"], data["N_charging"])
        delete!(labels, k)
    end
    # include subpath from a depot to itself
    null_subpath_cost = - κ[starting_node] - μ[starting_node]
    labels[starting_node][(starting_time, starting_charge)] = PathWithCost(
        subpaths = [
            SubpathWithCost(
                cost = null_subpath_cost,
                n_customers = data["n_customers"],
                starting_node = starting_node,
                starting_time = starting_time,
                starting_charge = starting_charge,
                arcs = [(starting_node, starting_node)],
            )
        ],
        cost = null_subpath_cost,
        explored = true,
    )
    return labels
end

function update_generated_subpaths_withcharge_from_path!(
    generated_subpaths_withcharge,
    path::PathWithCost,
)
    ## VERSION 1: include all subpaths in path if cumulative cost is negative
    if path.cost ≥ -1e-6
        return generated_subpaths_withcharge
    end
    ## VERSION 2: include all subpaths if any r.c. is negative
    # if all(s.cost ≥ -1e-6 for s in path.subpaths)
    #     return generated_subpaths_withcharge
    # end
    states = vcat(
        [(path.subpaths[1].starting_node, path.subpaths[1].starting_time, path.subpaths[1].starting_charge)],
        [
            (s.current_node, s.round_time, s.round_charge)
            for s in path.subpaths
        ]
    )
    for (ind, swc) in enumerate(path.subpaths)
        state_pair = (states[ind], states[ind+1])
        s = Subpath(swc)
        if state_pair in keys(generated_subpaths_withcharge)
            if !any(isequal(s, s1) for s1 in generated_subpaths_withcharge[state_pair])
                push!(generated_subpaths_withcharge[state_pair], s)
            end
        else
            generated_subpaths_withcharge[state_pair] = [s]
        end
    end
    
    return generated_subpaths_withcharge
end

function generate_subpaths_withcharge_from_paths(
    G,
    data,
    T_range,
    B_range,
    κ,
    μ,
    ν,
    ;
    charging_in_subpath::Bool = true,
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
)
    if !charging_in_subpath
        error()
    end
    generated_subpaths_withcharge = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Subpath},
    }()
    generated_paths_withcharge = Dict{Int, Any}()
    for starting_node in data["N_depots"]
        r = @timed find_smallest_reduced_cost_paths(
            starting_node, G, data, T_range, B_range, 
            κ, μ, ν,
            ;
            
            charge_to_full_only = charge_to_full_only,
            time_windows = time_windows,
            with_charging_cost = with_charging_cost,
        )
        labels = r.value
        generated_paths_withcharge[starting_node] = labels
        # remove those corresponding to positive reduced cost
        for (end_node, path_dict) in pairs(labels)
            for path in values(path_dict)
                update_generated_subpaths_withcharge_from_path!(
                    generated_subpaths_withcharge, path, 
                )
            end
        end
    end
    return generated_subpaths_withcharge, nothing, nothing, generated_paths_withcharge
end

function find_smallest_reduced_cost_subpaths_nocharge_notimewindows(
    starting_node, 
    G,
    data, 
    κ,
    μ,
    ν,
    ;
)
    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for i in 1:data["n_customers"]
        j = data["n_customers"] + i
        modified_costs[i,j] -= ν[i]
        modified_costs[j,i] -= ν[i]
    end
    # initialize set of labels
    if starting_node in data["N_depots"]
        initial_cost = - κ[starting_node]
    else
        initial_cost = 0.0
    end
    labels = Dict(
        n => OrderedDict{Float64, SubpathWithCost}()
        for n in union(
            data["N_pickups"], 
            data["N_charging"], 
            data["N_dropoffs"], 
            data["N_depots"],
        )
    )
    labels[starting_node][0.0] = SubpathWithCost(
        cost = initial_cost,
        n_customers = data["n_customers"],
        starting_node = starting_node,
        starting_time = 0.0, # sensible default
        starting_charge = data["B"],
    )
    # initialize queue of unexplored nodes
    Q = [starting_node]
    
    while length(Q) > 0
        i = popfirst!(Q)
        for s in values(labels[i])
            if s.explored
                continue
            else
                s.explored = true
            end
            for j in setdiff(outneighbors(G, i), i)
                add_to_queue = false
            
                # feasibility check
                if j in data["N_pickups"]
                    feasible_service = !(s.served[j])
                else
                    feasible_service = true
                end
                if !feasible_service 
                    continue
                end

                current_time = s.time + data["t"][i,j]
                current_charge = s.charge - data["q"][i,j]
                feasible = (
                    current_time ≤ data["T"] # feasible in terms of time horizon
                    && current_charge ≥ 0.0 # feasible in terms of charge horizon
                )
                if !feasible
                    continue
                end

                current_cost = s.cost + modified_costs[i,j]
                if j in data["N_depots"]
                    current_cost += - μ[j]
                end
               
                add_subpath = true
                if current_time in keys(labels[j])
                    if labels[j][current_time].cost ≤ current_cost
                        add_subpath = false
                    end
                end
                
                if add_subpath
                    add_to_queue = true
                    s_j = copy(s)
                    s_j.cost = current_cost
                    s_j.current_node = j
                    push!(s_j.arcs, (i, j))
                    s_j.time = current_time
                    s_j.charge = current_charge
                    if j in data["N_dropoffs"]
                        s_j.served[j-data["n_customers"]] = true
                    end
                    s_j.delta_time = 0.0
                    s_j.delta_charge = 0.0
                    s_j.end_time = current_time
                    s_j.end_charge = current_charge
                    s_j.round_time = current_time
                    s_j.round_charge = current_charge

                    labels[j][current_time] = s_j
                end
                if add_to_queue && !(j in Q) && !(j in union(data["N_charging"], data["N_depots"]))
                    push!(Q, j)
                end            
            end
        end
    end
    

    for node in union(data["N_pickups"], data["N_dropoffs"])
        delete!(labels, node)
    end
    for node in union(data["N_depots"], data["N_charging"])
        labels[node] = sort(labels[node])
    end

    return labels
end

function find_smallest_reduced_cost_subpaths_notimewindows(
    starting_node, 
    G,
    data, 
    T_range,
    B_range, 
    κ,
    μ,
    ν,
    ;
    charge_to_full_only::Bool = false,
    with_charging_cost::Bool = false,
)

    labels = find_smallest_reduced_cost_subpaths_nocharge_notimewindows(
        starting_node, G, data, 
        κ, μ, ν,
        ;
    )

    # generate charging options and perform dominance criteria
    new_labels = Dict{Int, SortedDict{Tuple{Float64, Float64}, SubpathWithCost}}(
        node => SortedDict{Tuple{Float64, Float64}, SubpathWithCost}(DoubleStateOrdering())
        for node in union(data["N_depots"], data["N_charging"])
    )
    for node in data["N_charging"]
        for (time, s) in labels[node]
            # TODO: speed up by sorting labels[node]
            for (dt, db, et, eb, rt, rb) in generate_charging_options(
                s.time, s.charge, data["μ"], T_range, B_range,
                ;
                charge_to_full_only = charge_to_full_only,
            )
                key = (rt, rb)
                new_cost = (
                    with_charging_cost ?
                        s.cost + data["charge_cost_coeff"] * dt :
                        s.cost 
                )
                if key in keys(new_labels[node])
                    if new_labels[node][key].cost ≤ new_cost
                        continue
                    end
                end
                s_copy = copy(s)
                s_copy.cost = new_cost
                s_copy.delta_time = dt
                s_copy.delta_charge = db
                s_copy.end_time = et
                s_copy.end_charge = eb
                s_copy.round_time = rt
                s_copy.round_charge = rb
                new_labels[node][key] = s_copy
            end
        end
    end
    for node in data["N_depots"]
        for (time, s) in labels[node]
            key = (s.end_time, s.end_charge)
            if key in keys(new_labels[node])
                if new_labels[node][key].cost ≤ s.cost
                    continue
                end
            end
            new_labels[node][key] = s
        end
    end

    return new_labels
end

function remove_dominated_subpaths_paths_withcharge!(
    collection::SortedDict{
        Tuple{Float64, Float64}, 
        Union{PathWithCost, SubpathWithCost},
        DoubleStateOrdering
    },
)
    st1 = beforestartsemitoken(collection)
    # while deleted
    while true
        st1 = advance((collection, st1))
        if status((collection, st1)) == 3 # past-end token
            break
        end
        k1, v1 = deref((collection, st1))
        for (st2, k2, v2) in semitokens(inclusive(collection, st1, lastindex(collection)))
            if (
                lt(dso, k1, k2)
                && v1.cost ≤ v2.cost
            )
                # println("$(k1[1]), $(k1[2]), $(v1.cost) dominates $(k2[1]), $(k2[2]), $(v2.cost)")
                if status((collection, st2)) == 1
                    delete!((collection, st2))
                end
            end
        end
    end
end

function generate_subpaths_withcharge_from_paths_notimewindows_V3(
    G,
    data, 
    T_range, 
    B_range,
    κ, 
    μ, 
    ν, 
    ;
    charge_to_full_only::Bool = false,
    with_charging_cost::Bool = false,
)

    function charge_to_specified_level(start_charge, end_charge, start_time, charging_rate)
        if end_charge ≤ start_charge
            return (0, 0, start_time, start_charge)
        end
        delta_charge = (end_charge - start_charge)
        delta_time = delta_charge / charging_rate
        end_time = start_time + delta_time
        return (
            round(delta_time, digits = 1), 
            delta_charge,
            round(end_time, digits = 1), 
            end_charge,
        )
    end

    function add_subpath_path_withcost_to_collection!(
        collection::Union{
            SortedDict{
                Tuple{Float64, Float64}, 
                PathWithCost,
                DoubleStateOrdering
            },
            SortedDict{
                Tuple{Float64, Float64}, 
                SubpathWithCost,
                DoubleStateOrdering
            },
        },
        k1::Tuple{Float64, Float64}, 
        v1::Union{PathWithCost, SubpathWithCost},
        ;
        on_timecharge::Bool = true,
        verbose::Bool = false,
    )
        added = true
        for (k2, v2) in collection
            # check if v2 dominates v1
            if v2.cost ≤ v1.cost
                if !on_timecharge || (on_timecharge && k2[1] ≤ k1[1] && k2[2] ≥ k1[2])
                    added = false
                    if verbose
                        println("$(k1[1]), $(k1[2]), $(v1.cost) dominated by $(k2[1]), $(k2[2]), $(v2.cost)")
                    end
                    break
                end
            # check if v1 dominates v2
            elseif v1.cost ≤ v2.cost
                if !on_timecharge || (on_timecharge && k1[1] ≤ k2[1] && k1[2] ≥ k2[2])
                    if verbose
                        println("$(k1[1]), $(k1[2]), $(v1.cost) dominates $(k2[1]), $(k2[2]), $(v2.cost)")
                    end
                    pop!(collection, k2)
                end
            end
        end
        if added
            if verbose
                println("$(k1[1]), $(k1[2]), $(v1.cost) added!")
            end
            insert!(collection, k1, v1)
        end
        return added
    end

    start_time = time()

    tso = TripleStateOrdering()
    dso = DoubleStateOrdering()

    base_labels = Dict{Int, Dict{Int64, Dict{Float64, SubpathWithCost}}}()
    for node in union(data["N_depots"], data["N_charging"])   
        base_labels[node] = find_smallest_reduced_cost_subpaths_nocharge_notimewindows(
            node, G, data, 
            κ, μ, ν,
            ;
        )
    end
    base_labels_wc = Dict(
        node1 => Dict(
            node2 => SortedDict(
                DoubleStateOrdering(),
                (time, s.charge) => s
                for (time, s) in pairs(base_labels[node1][node2])
            )
            for node2 in union(data["N_depots"], data["N_charging"])
        )
        for node1 in union(data["N_depots"], data["N_charging"])
    )
    time1 = time()

    nondom_base_labels = Dict(
        node1 => Dict(
            node2 => SortedDict{Tuple{Float64, Float64}, SubpathWithCost}(dso)
            for node2 in union(data["N_depots"], data["N_charging"])
        )
        for node1 in union(data["N_depots"], data["N_charging"])
    )
    for node1 in union(data["N_depots"], data["N_charging"]), 
        node2 in union(data["N_depots"], data["N_charging"])
        for (k1, v1) in base_labels_wc[node1][node2]
            # 1. check if v1 is not dominated
            # 2. check if v1 dominates anyone
            add_subpath_path_withcost_to_collection!(
                nondom_base_labels[node1][node2],
                k1, v1,
                ;
            )
        end
    end

    time2 = time()
    full_labels = Dict(
        # starting depot
        starting_node => Dict(
            # current state
            current_node => SortedDict{Tuple{Float64, Float64}, PathWithCost}(dso)
            for current_node in union(data["N_charging"], data["N_depots"])
        )
        for starting_node in data["N_depots"]
    )
    for depot in data["N_depots"]
        full_labels[depot][depot][(0.0, data["B"])] = PathWithCost(
            subpaths = [], served = zeros(Int, data["n_customers"]), cost = 0.0
        )
    end
    unexplored_states = SortedSet(
        tso, 
        [
            (depot, 0.0, data["B"]) 
            for depot in data["N_depots"]
        ]
    )
    time3 = time()

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        for starting_node in data["N_depots"]
            if !((state[2], state[3]) in keys(full_labels[starting_node][state[1]]))
                continue
            end
            current_path = full_labels[starting_node][state[1]][(state[2], state[3])]
            # where do you want to go next
            for next_node in union(data["N_depots"], data["N_charging"])
                # what subpaths bring you there
                for s in values(nondom_base_labels[state[1]][next_node])
                    if ( # empty subpath
                        s.current_node in data["N_depots"] 
                        && s.round_time == 0.0
                        && s.round_charge == data["B"]
                    )
                        continue
                    end
                    if ( # each customer served at most once
                        any(s.served + current_path.served .> 1)
                    )
                        continue
                    end
                    # compute amount of charge required, charging required
                    charge_required = data["B"] - s.charge 
                    (
                        delta_time, 
                        delta_charge,
                        end_time, 
                        end_charge,
                    ) = charge_to_specified_level(state[3], charge_required, state[2], data["μ"])
                    if end_time + s.time > data["T"]
                        continue
                    end
                    if end_charge - charge_required < 0
                        continue
                    end

                    current_cost = current_path.cost
                    if length(current_path.subpaths) > 0
                        s_prev = copy(current_path.subpaths[end])
                        if with_charging_cost
                            s_prev.cost += data["charge_cost_coeff"] * delta_time
                            current_cost += data["charge_cost_coeff"] * delta_time
                        end
                        s_prev.delta_time = delta_time
                        s_prev.delta_charge = delta_charge
                        s_prev.end_time = end_time
                        s_prev.end_charge = end_charge
                        s_prev.round_time = end_time
                        s_prev.round_charge = end_charge
                    end
                    s_new = copy(s)
                    current_cost += s_new.cost
                    s_new.starting_time = end_time
                    s_new.starting_charge = end_charge
                    s_new.time = s.time + end_time
                    s_new.charge = end_charge - (data["B"] - s.charge)
                    s_new.end_time = s.end_time + end_time
                    s_new.end_charge = end_charge - (data["B"] - s.end_charge)
                    s_new.round_time = s.end_time + end_time
                    s_new.round_charge = end_charge - (data["B"] - s.end_charge)
                    
                    key = (s_new.end_time, s_new.end_charge)
                    
                    if length(current_path.subpaths) > 0
                        new_path = PathWithCost(
                            subpaths = vcat(
                                current_path.subpaths[1:end-1], 
                                [s_prev, s_new],
                            ),
                            cost = current_cost,
                        )
                    else 
                        new_path = PathWithCost(
                            subpaths = [s_new],
                            cost = current_cost,
                        )
                    end
    
                    add_next_state = false
                    if add_subpath_path_withcost_to_collection!(
                        full_labels[starting_node][next_node],
                        key, new_path,
                        ;
                        on_timecharge = true,
                    )
                        # println("Extending: $(starting_node) -> $(state[1]), $(state[2]), $(state[3]) along $(state[1]) -> $(next_node), $(key[1]), $(key[2]); cost $(new_path.cost)")
                        add_next_state = true
                        # println()
                    else
                        # println("Dominated: $(starting_node) -> $(state[1]), $(state[2]), $(state[3]) along $(state[1]) -> $(next_node), $(key[1]), $(key[2]); cost $(new_path.cost)")
                        nothing
                    end
                    
                    next_state = (next_node, key[1], key[2])
                    if (
                        add_next_state  
                        && next_node in data["N_charging"] 
                        && !(next_state in unexplored_states)
                    )
                        insert!(unexplored_states, next_state)
                    end
                end
            end
        end
    end
    for depot in data["N_depots"]
        push!(
            full_labels[depot][depot], 
            (0.0, data["B"]) => PathWithCost(
                subpaths = [
                    nondom_base_labels[depot][depot][(0.0, data["B"])]
                ],
                cost = nondom_base_labels[depot][depot][(0.0, data["B"])].cost,
            )
        )
    end

    for starting_node in data["N_depots"]
        for end_node in data["N_charging"]
            delete!(full_labels[starting_node], end_node)
        end
    end

    time4 = time()

    generated_subpaths_withcharge = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Subpath},
    }()
    for starting_node in data["N_depots"], end_node in data["N_depots"]
        for path in values(full_labels[starting_node][end_node])
            update_generated_subpaths_withcharge_from_path!(
                generated_subpaths_withcharge, path,
            )
        end
    end

    time5 = time()

    return generated_subpaths_withcharge, nothing, (time1 - start_time, time2 - time1, time3 - time2, time4 - time3, time5 - time4), full_labels
end

function generate_subpaths_withcharge_from_paths_notimewindows_V2(
    G,
    data, 
    T_range, 
    B_range,
    κ, 
    μ, 
    ν, 
    ;
    charge_to_full_only::Bool = false,
    with_charging_cost::Bool = false,
)

    function add_subpath_path_withcost_to_collection!(
        collection::Union{
            SortedDict{
                Tuple{Float64, Float64}, 
                PathWithCost,
                DoubleStateOrdering
            },
            SortedDict{
                Tuple{Float64, Float64}, 
                SubpathWithCost,
                DoubleStateOrdering
            },
        },
        k1::Tuple{Float64, Float64}, 
        v1::Union{PathWithCost, SubpathWithCost},
        ;
        on_timecharge::Bool = true,
        verbose::Bool = false,
    )
        added = true
        for (k2, v2) in collection
            # check if v2 dominates v1
            if v2.cost ≤ v1.cost
                if !on_timecharge || (on_timecharge && k2[1] ≤ k1[1] && k2[2] ≥ k1[2])
                    added = false
                    if verbose
                        println("$(k1[1]), $(k1[2]), $(v1.cost) dominated by $(k2[1]), $(k2[2]), $(v2.cost)")
                    end
                    break
                end
            # check if v1 dominates v2
            elseif v1.cost ≤ v2.cost
                if !on_timecharge || (on_timecharge && k1[1] ≤ k2[1] && k1[2] ≥ k2[2])
                    if verbose
                        println("$(k1[1]), $(k1[2]), $(v1.cost) dominates $(k2[1]), $(k2[2]), $(v2.cost)")
                    end
                    pop!(collection, k2)
                end
            end
        end
        if added
            if verbose
                println("$(k1[1]), $(k1[2]), $(v1.cost) added!")
            end
            insert!(collection, k1, v1)
        end
        return added
    end

    start_time = time()

    tso = TripleStateOrdering()
    dso = DoubleStateOrdering()
    base_labels = Dict{Int, Dict{Int64, SortedDict{Tuple{Float64, Float64}, SubpathWithCost}}}()
    for node in union(data["N_depots"], data["N_charging"])   
        base_labels[node] = find_smallest_reduced_cost_subpaths_notimewindows(
            node, G, data, T_range, B_range, κ, μ, ν,
            ;
            charge_to_full_only = charge_to_full_only,
            with_charging_cost = with_charging_cost,
        )
    end
    time1 = time()

    nondom_base_labels = Dict(
        node1 => Dict(
            node2 => SortedDict{Tuple{Float64, Float64}, SubpathWithCost}(dso)
            for node2 in union(data["N_depots"], data["N_charging"])
        )
        for node1 in union(data["N_depots"], data["N_charging"])
    )
    for node1 in union(data["N_depots"], data["N_charging"]), node2 in union(data["N_depots"], data["N_charging"])
        for (k1, v1) in base_labels[node1][node2]
            # 1. check if v1 is not dominated
            # 2. check if v1 dominates anyone
            add_subpath_path_withcost_to_collection!(
                nondom_base_labels[node1][node2],
                k1, v1,
                ;
            )
        end
    end

    time2 = time()
    full_labels = Dict(
        # starting depot
        starting_node => Dict(
            # current state
            current_node => SortedDict{Tuple{Float64, Float64}, PathWithCost}(dso)
            for current_node in union(data["N_charging"], data["N_depots"])
        )
        for starting_node in data["N_depots"]
    )
    for depot in data["N_depots"]
        full_labels[depot][depot][(0.0, data["B"])] = PathWithCost(
            subpaths = [], served = zeros(Int, data["n_customers"]), cost = 0.0
        )
    end
    unexplored_states = SortedSet(
        tso, 
        [
            (depot, 0.0, data["B"]) 
            for depot in data["N_depots"]
        ]
    )
    time3 = time()

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        # where do you want to go next
        for next_node in union(data["N_depots"], data["N_charging"])
            # what subpaths bring you there
            for s in values(nondom_base_labels[state[1]][next_node])
                if !(
                    # time required + time you start at ≤ time horizon
                    s.round_time + state[2] ≤ data["T"]
                    # charge required (note: this is different from charge taken)
                    # ≤ charge you start at
                    && (data["B"]) - s.charge ≤ state[3]
                )
                    # our list of subpaths is sorted by increasing time and decreasing charge
                    break
                end
                if ( # empty subpath
                    s.current_node in data["N_depots"] 
                    && s.round_time == 0.0
                    && s.round_charge == data["B"]
                )
                    continue
                end
                # translate subpath
                s_new = copy(s)
                s_new.starting_time = state[2]
                s_new.starting_charge = state[3]
                s_new.time = s.time + state[2]
                s_new.charge = s.charge + (state[3] - data["B"])
                s_new.end_time = s.end_time + state[2]
                s_new.end_charge = s.end_charge + (state[3] - data["B"])
                s_new.round_time = s.round_time + state[2]
                s_new.round_charge = s.round_charge + (state[3] - data["B"])
    
                if next_node in data["N_depots"]
                    key = (s_new.end_time, s_new.end_charge)
                else
                    key = (s_new.round_time, s_new.round_charge)
                end
                add_next_state = false
                for starting_node in data["N_depots"]
                    if !((state[2], state[3]) in keys(full_labels[starting_node][state[1]]))
                        continue
                    end
                    current_cost = full_labels[starting_node][state[1]][(state[2], state[3])].cost + s_new.cost
                    new_path = PathWithCost(
                        subpaths = vcat(
                            full_labels[starting_node][state[1]][(state[2], state[3])].subpaths, 
                            [s_new],
                        ),
                        cost = current_cost,
                    )
                    if add_subpath_path_withcost_to_collection!(
                        full_labels[starting_node][next_node],
                        key, new_path,
                        ;
                        on_timecharge = true,
                    )
                        # println("Extending: $(starting_node) -> $(state[1]), $(state[2]), $(state[3]) along $(state[1]) -> $(next_node), $(key[1]), $(key[2])")
                        add_next_state = true
                    end
                end
                next_state = (next_node, key[1], key[2])
                if (
                    add_next_state  
                    && next_node in data["N_charging"] 
                    && !(next_state in unexplored_states)
                )
                    insert!(unexplored_states, next_state)
                end
            end
        end
        # TODO: prevent long number of redundant operations by sorting states.
    end
    for depot in data["N_depots"]
        push!(
            full_labels[depot][depot], 
            (0.0, data["B"]) => PathWithCost(
                subpaths = [
                    nondom_base_labels[depot][depot][(0.0, data["B"])]
                ],
                cost = nondom_base_labels[depot][depot][(0.0, data["B"])].cost,
            )
        )
    end

    for starting_node in data["N_depots"]
        for end_node in data["N_charging"]
            delete!(full_labels[starting_node], end_node)
        end
    end

    time4 = time()

    generated_subpaths_withcharge = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Subpath},
    }()
    for starting_node in data["N_depots"], end_node in data["N_depots"]
        for path in values(full_labels[starting_node][end_node])
            update_generated_subpaths_withcharge_from_path!(
                generated_subpaths_withcharge, path,
            )
        end
    end

    time5 = time()

    return generated_subpaths_withcharge, nothing, (time1 - start_time, time2 - time1, time3 - time2, time4 - time3, time5 - time4), full_labels
end

function generate_subpaths_withcharge_from_paths_notimewindows(
    G,
    data, 
    T_range, 
    B_range,
    κ, 
    μ, 
    ν, 
    ;
    charge_to_full_only::Bool = false,
    with_charging_cost::Bool = false,
)
    start_time = time()

    tso = TripleStateOrdering()
    dso = DoubleStateOrdering()
    base_labels = Dict()
    for node in union(data["N_depots"], data["N_charging"])   
        base_labels[node] = find_smallest_reduced_cost_subpaths_notimewindows(
            node, G, data, T_range, B_range, κ, μ, ν,
            ;
            charge_to_full_only = charge_to_full_only,
            with_charging_cost = with_charging_cost,
        )
    end

    time1 = time()
    nondom_base_labels = base_labels
    time2 = time()

    full_labels = Dict(
        # starting depot
        starting_node => Dict(
            # current state
            current_node => SortedDict{Tuple{Float64, Float64}, PathWithCost}(dso)
            for current_node in union(data["N_charging"], data["N_depots"])
        )
        for starting_node in data["N_depots"]
    )
    for depot in data["N_depots"]
        full_labels[depot][depot][(0.0, data["B"])] = PathWithCost(
            subpaths = [], served = zeros(Int, data["n_customers"]), cost = 0.0
        )
    end
    unexplored_states = SortedSet(
        tso, 
        [
            (depot, 0.0, data["B"]) 
            for depot in data["N_depots"]
        ]
    )
    time3 = time()

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        # where do you want to go next
        for next_node in union(data["N_depots"], data["N_charging"])
            # what subpaths bring you there
            for s in values(nondom_base_labels[state[1]][next_node])
                if !(
                    # time required + time you start at ≤ time horizon
                    s.round_time + state[2] ≤ data["T"]
                    # charge required (note: this is different from charge taken)
                    # ≤ charge you start at
                    && (data["B"]) - s.charge ≤ state[3]
                )
                    # our list of subpaths is sorted by increasing time and decreasing charge
                    break
                end
                if ( # empty subpath
                    s.current_node in data["N_depots"] 
                    && s.round_time == 0.0
                    && s.round_charge == data["B"]
                )
                    continue
                end
                # translate subpath
                s_new = copy(s)
                s_new.starting_time = state[2]
                s_new.starting_charge = state[3]
                s_new.time = s.time + state[2]
                s_new.charge = s.charge + (state[3] - data["B"])
                s_new.end_time = s.end_time + state[2]
                s_new.end_charge = s.end_charge + (state[3] - data["B"])
                s_new.round_time = s.round_time + state[2]
                s_new.round_charge = s.round_charge + (state[3] - data["B"])
    
                if next_node in data["N_depots"]
                    key = (s_new.end_time, s_new.end_charge)
                else
                    key = (s_new.round_time, s_new.round_charge)
                end
                add_next_state = false
                for starting_node in data["N_depots"]
                    if !((state[2], state[3]) in keys(full_labels[starting_node][state[1]]))
                        continue
                    end
                    current_cost = full_labels[starting_node][state[1]][(state[2], state[3])].cost + s_new.cost
                    new_path = PathWithCost(
                        subpaths = vcat(
                            full_labels[starting_node][state[1]][(state[2], state[3])].subpaths, 
                            [s_new],
                        ),
                        cost = current_cost,
                    )
                    if key in keys(full_labels[starting_node][next_node])
                        if current_cost ≥ full_labels[starting_node][next_node][key].cost
                            continue
                        end
                    end
                    full_labels[starting_node][next_node][key] = new_path
                    add_next_state = true
                end
                next_state = (next_node, key[1], key[2])
                if (
                    add_next_state 
                    && next_node in data["N_charging"] 
                    && !(next_state in unexplored_states)
                )
                    insert!(unexplored_states, next_state)
                end
            end
        end
        # TODO: prevent long number of redundant operations by sorting states.
    end
    for depot in data["N_depots"]
        push!(
            full_labels[depot][depot], 
            (0.0, data["B"]) => PathWithCost(
                subpaths = [
                    nondom_base_labels[depot][depot][(0.0, data["B"])]
                ],
                cost = nondom_base_labels[depot][depot][(0.0, data["B"])].cost,
            )
        )
    end

    for starting_node in data["N_depots"]
        for end_node in data["N_charging"]
            delete!(full_labels[starting_node], end_node)
        end
    end

    time4 = time()

    generated_subpaths_withcharge = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Subpath},
    }()
    for starting_node in data["N_depots"], end_node in data["N_depots"]
        for path in values(full_labels[starting_node][end_node])
            update_generated_subpaths_withcharge_from_path!(
                generated_subpaths_withcharge, path,
            )
        end
    end

    time5 = time()

    return generated_subpaths_withcharge, nothing, (time1 - start_time, time2 - time1, time3 - time2, time4 - time3, time5 - time4), full_labels
end

function find_smallest_reduced_cost_subpaths(
    starting_node,
    starting_time,
    starting_charge,
    G,
    data,
    T_range,
    B_range,
    # dual values associated with flow conservation constraints
    κ,
    λ,
    μ,
    # dual values associated with customer service constraints
    ν,
    ;
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
)
    """
    Generates feasible subpaths from a state 
    (`starting_node`, `starting_time`, `starting_charge`) to all nodes;
    for each node, generates the one with the smallest reduced cost 
    based on the reduced costs (from the dual variables).
    """
    # initialize modified arc costs, subtracting values of dual variables
    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for i in 1:data["n_customers"]
        j = data["n_customers"] + i
        modified_costs[i,j] -= ν[i]
        modified_costs[j,i] -= ν[i]
    end
    # initialize set of labels
    initial_cost = 0.0
    if starting_node in data["N_depots"]
        if starting_time == 0.0 && starting_charge == data["B"]
            initial_cost = initial_cost - κ[starting_node]
        end
    elseif starting_node in data["N_charging"]
        initial_cost = initial_cost - λ[(starting_node, starting_time, starting_charge)]
    end
    labels = Dict()
    labels[starting_node] = Dict(
        (starting_time, starting_charge) => SubpathWithCost(
            cost = initial_cost,
            n_customers = data["n_customers"],
            starting_node = starting_node,
            starting_time = starting_time,
            starting_charge = starting_charge,
        )
    )
    # initialize queue of unexplored nodes
    Q = [starting_node]

    while !isempty(Q)
        i = popfirst!(Q)
        # iterate over all out-neighbors of i
        for ((now_time, now_charge), s) in pairs(labels[i])
            if s.explored
                continue
            else 
                s.explored = true
            end
            for j in setdiff(outneighbors(G, i), i)
                add_to_queue = false
                # feasibility check
                if j in data["N_pickups"]
                    feasible_service = !(s.served[j])
                else
                    feasible_service = true
                end
                if !feasible_service 
                    continue
                end
                if time_windows
                    current_time = max(
                        data["α"][j], 
                        now_time + data["t"][i,j],
                    )
                    feasible_timewindow = (current_time ≤ data["β"][j])
                else
                    current_time = now_time + data["t"][i,j]
                    feasible_timewindow = (current_time ≤ data["T"])
                end
                current_charge = now_charge - data["q"][i,j]
                feasible_charge = (current_charge ≥ 0.0)
                feasible = (
                    feasible_timewindow
                    && feasible_charge
                )
                if !feasible
                    continue
                end
                current_cost = s.cost + modified_costs[i,j]

                if j in data["N_charging"]
                    charging_options = generate_charging_options(
                        current_time, current_charge, 
                        data["μ"], T_range, B_range,
                        ;
                        
                        charge_to_full_only = charge_to_full_only,
                    )
                    if length(charging_options) == 0
                        continue
                    end
                else
                    if j in data["N_depots"]
                        current_cost = current_cost - μ[j]
                    end
                    charging_options = [(
                        0, 0, current_time, current_charge, 
                        current_time, current_charge
                    )]
                end

                for (delta_time, delta_charge, end_time, end_charge, store_time, store_charge) in charging_options
                    # determine if label ought to be updated
                    this_subpath_cost = current_cost 
                    if j in data["N_charging"]
                        this_subpath_cost += λ[(j, store_time, store_charge)]
                        if with_charging_cost
                            this_subpath_cost += data["charge_cost_coeff"] * delta_time
                        end
                    end
                    
                    if j in keys(labels)
                        if (store_time, store_charge) in keys(labels[j])
                            if labels[j][(store_time, store_charge)].cost ≤ this_subpath_cost
                                continue
                            end
                        end
                    end
                    # update labels 
                    add_to_queue = true
                    s_j = copy(s)
                    s_j.cost = this_subpath_cost
                    s_j.current_node = j
                    push!(s_j.arcs, (i, j))
                    s_j.time = current_time
                    s_j.charge = current_charge
                    if j in data["N_dropoffs"]
                        s_j.served[j-data["n_customers"]] = true
                    end
                    s_j.delta_time = delta_time
                    s_j.delta_charge = delta_charge
                    s_j.end_time = end_time
                    s_j.end_charge = end_charge
                    s_j.round_time = dceil(end_time, T_range)
                    s_j.round_charge = dfloor(end_charge, B_range)
                    if !(j in keys(labels))
                        labels[j] = Dict{Tuple{Float64, Float64}, SubpathWithCost}()
                    end
                    labels[j][(store_time, store_charge)] = s_j
                end
                # add node j to the list of unexplored nodes
                if add_to_queue && !(j in Q) && !(j in union(data["N_charging"], data["N_depots"]))
                    push!(Q, j)
                end
            end
        end
    end
    # remove labels corresponding to pickup and dropoff nodes
    for k in union(data["N_pickups"], data["N_dropoffs"])
        delete!(labels, k)
    end
    # include subpath from a depot to itself, if it is better
    if starting_node in data["N_depots"]
        null_subpath_cost = initial_cost - μ[starting_node]
        if !(
            (starting_time, starting_charge) in keys(labels[starting_node])
            && labels[starting_node][(starting_time, starting_charge)].cost 
            ≤ null_subpath_cost
        )
            labels[starting_node][(starting_time, starting_charge)] = SubpathWithCost(
                cost = null_subpath_cost,
                n_customers = data["n_customers"],
                starting_node = starting_node,
                starting_time = starting_time,
                starting_charge = starting_charge,
                arcs = [(starting_node, starting_node)],
            )
        end
    end
    return labels
end

function generate_subpaths_withcharge(
    G,
    data,
    T_range,
    B_range,
    κ,
    λ,
    μ,
    ν,
    ;
    charging_in_subpath::Bool = true,
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
)
    """
    Generates subpaths (inclusive of charging) with 
    negative reduced cost (determined by dual variables)
    for all (I, j) pairs.
    """
    # FIXME
    if !charging_in_subpath
        error()
    end
    generated_subpaths_withcharge = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Subpath},
    }()
    smallest_reduced_costs = Dict{Tuple, Float64}()
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
        # generate subpaths from this starting state
        state1 = (starting_node, starting_time, starting_charge)
        r = @timed find_smallest_reduced_cost_subpaths(
            starting_node, starting_time, starting_charge,
            G, data, T_range, B_range,
            κ, λ, μ, ν,
            ;
            
            charge_to_full_only = charge_to_full_only,
            time_windows = time_windows,
            with_charging_cost = with_charging_cost,
        )
        labels = r.value
        # remove those corresponding to positive reduced cost
        smallest_reduced_cost = Inf
        for (end_node, subpath_dict) in pairs(labels)
            for s in values(subpath_dict)
                # Toss out subpaths that have no arcs
                if length(s.arcs) == 0
                    if !(starting_node in data["N_depots"]
                        && starting_time == 0.0
                        && starting_charge == data["B"]
                    )
                        continue
                    end
                end
                if s.cost ≥ -1e-6
                    continue
                end
                if s.cost < smallest_reduced_cost
                    smallest_reduced_cost = s.cost
                end
                state2 = (s.current_node, s.round_time, s.round_charge)
                generated_subpaths_withcharge[(state1, state2)] = [Subpath(s)]
            end
        end
        smallest_reduced_costs[state1] = smallest_reduced_cost
    end
    return generated_subpaths_withcharge, smallest_reduced_costs, nothing
end

function compute_subpath_cost(
    data,
    s::Subpath,
    M::Float64 = 1e8,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
    verbose::Bool = false,
)
    if s.artificial 
        return M
    elseif length(s.arcs) == 0
        return 0
    end

    travel_cost = data["travel_cost_coeff"] * sum(
        data["c"][a...] for a in s.arcs
    )
    charge_cost = 0.0
    if with_charging_cost
        charge_cost += data["charge_cost_coeff"] * s.delta_time
    end
    if verbose
        println("Travel cost:          $(travel_cost)")
        println("Charge cost:          $(charge_cost)")
    end
    return travel_cost + charge_cost
end

function compute_subpath_costs(
    data,
    all_subpaths,
    M::Float64 = 1e8,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
)
    subpath_costs = Dict(
        key => [
            compute_subpath_cost(
                data, s, M
                ;
                with_charging_cost = with_charging_cost,
                time_windows = time_windows,
            )
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
    ;
    charge_to_full_only::Bool = false,
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
                            data["μ"],
                            T_range,
                            B_range,
                            ;
                            
                            charge_to_full_only = charge_to_full_only,    
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
    integral::Bool = true,
    charging_in_subpath::Bool = false,
    monotonic = false,
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
    if integral
        @variable(model, z[
            k=keys(all_subpaths), 
            p=1:m; p ≤ length(all_subpaths[k])
        ], Bin);
        if !charging_in_subpath
            @variable(model, y[all_charging_arcs], Bin)
        end
    else 
        @variable(model, z[
            k=keys(all_subpaths), 
            p=1:m; p ≤ length(all_subpaths[k])
        ] ≥ 0);
        if !charging_in_subpath
            @variable(model, y[all_charging_arcs] ≥ 0)
        end
    end

    if monotonic
        @variable(model, α_upper[
            data["N_charging"],
            B_range,
        ] ≥ 0)
        @variable(model, α_lower[
            data["N_charging"],
            B_range,
        ] ≥ 0)
        @variable(model, β_upper[
            data["N_charging"],
            T_range,
        ] ≥ 0)
        @variable(model, β_lower[
            data["N_charging"],
            T_range,
        ] ≥ 0)
    end

    # Constraint (4b): number of subpaths starting at depot i
    # at time 0 with full charge to be v_startᵢ
    @constraint(
        model,
        κ[i in data["N_depots"]],
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
    flow_conservation_constrs = Dict()
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
                    if monotonic
                        # Constraint (19c) Relaxed flow conservation at charging stations
                        # (here the other node can be either a depot or charging station)
                        flow_conservation_constrs[state1] = @constraint(
                            model,
                            flow_conservation_exprs[(state1,"out_sp")]
                            - flow_conservation_exprs[(state1,"in_sp")]
                            + α_lower[n1,b1] - α_upper[n1,b1]
                            + β_upper[n1,t1] - β_lower[n1,t1]
                            == 0.0
                        )
                    else
                        # Constraint (4c) Flow conservation at charging stations
                        # (here the other node can be either a depot or charging station)
                        flow_conservation_constrs[state1] = @constraint(
                            model,
                            flow_conservation_exprs[(state1,"out_sp")]
                            == flow_conservation_exprs[(state1,"in_sp")]
                        )
                    end
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
                    if monotonic # FIXME: untested
                        flow_conservation_constrs[state1] = @constraint(
                            model,
                            flow_conservation_exprs[(state1,"out_sp")]
                            + flow_conservation_exprs[(state1,"out_a")]
                            - flow_conservation_exprs[(state1,"in_sp")]
                            - flow_conservation_exprs[(state1,"in_a")]
                            + α_lower[n1,b1] - α_upper[n1,b1]
                            + β_upper[n1,t1] - β_lower[n1,t1]
                            == 0.0
                        )
                    else
                        # FIXME: make more efficient
                        # Constraint (7b) Flow conservation at charging stations
                        # (here the other node can be either a depot or charging station)
                        # (here the arcs can be subpaths or charging arcs)
                        flow_conservation_constrs[state1] = @constraint(
                            model,
                            flow_conservation_exprs[(state1,"out_sp")]
                            + flow_conservation_exprs[(state1,"out_a")]
                            == flow_conservation_exprs[(state1,"in_sp")]
                            + flow_conservation_exprs[(state1,"in_a")]
                        )
                    end
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

    # Constraint (4d) Number of subpaths ending at each depot n2 
    # is at least the number of vehicles required
    @constraint(
        model,
        μ[n2 in data["N_depots"]],
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

    # Constraint (4e) Serving all customers exactly once across all subpaths
    @constraint(
        model,
        ν[j in data["N_pickups"]],
        sum(
            sum(
                subpath_service[((state1, state2),j)][p] * z[(state1, state2),p]
                for p in 1:length(all_subpaths[(state1, state2)])
            )
            for (state1, state2) in keys(all_subpaths)
        ) == 1
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
        "time_taken" => round(time_taken, digits=3),
        "constraint_time_taken" => round(constraint_time_taken, digits=3),
        "solution_time_taken" => round(solution_time_taken, digits=3),
    )
    results = Dict(
        "objective" => objective_value(model),
        "z" => value.(z),
    )
    if !integral
        # retrieve dual solutions
        results["κ"] = Dict(zip(data["N_depots"], dual.(model[:κ]).data))
        results["λ"] = Dict(
            k => dual(v) 
            for (k, v) in pairs(flow_conservation_constrs)
        )
        results["μ"] = Dict(zip(data["N_depots"], dual.(model[:μ]).data))
        results["ν"] = dual.(model[:ν]).data
    end
    if !charging_in_subpath
        results["y"] = value.(y)
    end
    if monotonic
        results["α_upper"] = value.(model[:α_upper]).data
        results["α_lower"] = value.(model[:α_lower]).data
        results["β_upper"] = value.(model[:β_upper]).data
        results["β_lower"] = value.(model[:β_lower]).data
    end
    return results, params
end

function subpath_formulation_column_generation_from_paths(
    G,
    data,
    T_range,
    B_range,
    ;
    charging_in_subpath::Bool = true,
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
    with_heuristic::Bool = true,
    verbose::Bool = false,
    time_limit::Float64 = Inf,
)
    function add_message!(
        printlist::Vector, 
        message::String, 
        verbose::Bool,
    )
        push!(printlist, message)
        if verbose
            print(message)
        end
    end

    start_time = time()

    artificial_subpaths = generate_artificial_subpaths(data)
    if with_heuristic
        some_subpaths = construct_heuristic_subpaths(data, T_range, B_range)
        for (key, artificial_subpath_list) in pairs(artificial_subpaths)
            if !(key in keys(some_subpaths))
                some_subpaths[key] = []
            end
            append!(some_subpaths[key], artificial_subpath_list)
        end
    else
        some_subpaths = deepcopy(artificial_subpaths)
    end
    mp_results = Dict()
    params = Dict()
    params["number_of_subpaths"] = [sum(length(v) for v in values(some_subpaths))]
    params["objective"] = Float64[]
    params["κ"] = Dict{Int, Float64}[]
    params["μ"] = Dict{Int, Float64}[]
    params["ν"] = Vector{Float64}[]
    params["lp_relaxation_time_taken"] = Float64[]
    params["lp_relaxation_solution_time_taken"] = Float64[]
    params["lp_relaxation_constraint_time_taken"] = Float64[]
    params["sp_total_time_taken"] = Float64[]
    params["number_of_current_subpaths"] = Int[]
    printlist = []

    counter = 0
    converged = false
    add_message!(
        printlist,
        @sprintf(
            """
            Starting column generation.
            # customers:                %2d
            # depots:                   %2d
            # charging stations:        %2d
            # vehicles:                 %2d

            """,
            data["n_customers"],
            data["n_depots"],
            data["n_charging"],
            data["n_vehicles"],
        ),
        verbose,
    )
    add_message!(
        printlist,
        @sprintf("              |  Objective | # subpaths | Time (LP) | Time (SP) | # new subpaths \n"),
        verbose,
    )
    while (
        !converged
        && time_limit > time() - start_time
    )
        counter += 1
        subpath_costs = compute_subpath_costs(
            data, 
            some_subpaths,
            ;
            with_charging_cost = with_charging_cost,
            time_windows = time_windows,
        )
        subpath_service = compute_subpath_service(
            data, 
            some_subpaths,
        )
        mp_results, mp_params = @suppress subpath_formulation(
            data, 
            some_subpaths,
            subpath_costs, 
            subpath_service,
            T_range,
            B_range,
            ;
            integral = false,
            charging_in_subpath = charging_in_subpath,
        )

        push!(params["objective"], mp_results["objective"])
        push!(params["κ"], mp_results["κ"])
        push!(params["μ"], mp_results["μ"])
        push!(params["ν"], mp_results["ν"])
        push!(params["lp_relaxation_time_taken"], mp_params["time_taken"])
        push!(params["lp_relaxation_constraint_time_taken"], mp_params["constraint_time_taken"])
        push!(params["lp_relaxation_solution_time_taken"], mp_params["solution_time_taken"])

        # generate subpaths
        generate_subpaths_result = @timed generate_subpaths_withcharge_from_paths(
            G, data, T_range, B_range,
            mp_results["κ"], mp_results["μ"], mp_results["ν"],
            ;
            charging_in_subpath = charging_in_subpath,
            
            charge_to_full_only = charge_to_full_only,
            time_windows = time_windows,
            with_charging_cost = with_charging_cost,
        )
        (current_subpaths, _, _) = generate_subpaths_result.value

        push!(
            params["sp_total_time_taken"],
            round(generate_subpaths_result.time, digits=3)
        )
        if length(current_subpaths) == 0
            push!(params["number_of_current_subpaths"], 0)
            converged = true
        else
            push!(
                params["number_of_current_subpaths"],
                sum(length(v) for v in values(current_subpaths))
            )
        end
        for key in keys(current_subpaths)
            if !(key in keys(some_subpaths))
                some_subpaths[key] = current_subpaths[key]
            else
                for s_new in current_subpaths[key]
                    if !any(isequal(s_new, s) for s in some_subpaths[key])
                        push!(some_subpaths[key], s_new)
                    end
                end
            end
        end
        push!(
            params["number_of_subpaths"], 
            sum(length(v) for v in values(some_subpaths))
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %9.3f | %9.3f | %14d \n", 
                counter,
                params["objective"][counter],
                params["number_of_subpaths"][counter],
                params["lp_relaxation_time_taken"][counter],
                params["sp_total_time_taken"][counter],
                params["number_of_current_subpaths"][counter],
            ),
            verbose,
        )
        if length(params["number_of_subpaths"]) > 1
            if params["number_of_subpaths"][end-1] == params["number_of_subpaths"][end]
                break
            end
        end
    end
    results = Dict(
        "objective" => mp_results["objective"],
        "z" => mp_results["z"],
        "κ" => mp_results["κ"],
        "μ" => mp_results["μ"],
        "ν" => mp_results["ν"],
    )
    params["counter"] = counter
    params["converged"] = converged
    end_time = time()
    time_taken = round(end_time - start_time, digits=3)
    params["time_taken"] = time_taken
    params["time_limit_reached"] = (time_taken > time_limit)

    for message in [
        @sprintf("\n")
        @sprintf("Objective: \t\t%.4e\n", mp_results["objective"])
        @sprintf("Total time (LP): \t%10.3f s\n", sum(params["lp_relaxation_time_taken"]))
        @sprintf("Total time (SP): \t%10.3f s\n", sum(params["sp_total_time_taken"]))
        @sprintf("Total time: \t\t%10.3f s\n", time_taken)
    ]
        add_message!(printlist, message, verbose)
    end
    return results, params, printlist, some_subpaths
end

function subpath_formulation_column_generation(
    G,
    data,
    T_range,
    B_range,
    ;
    charging_in_subpath::Bool = true,
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
    with_heuristic::Bool = true,
    verbose::Bool = false,
)
    function add_message!(
        printlist::Vector, 
        message::String, 
        verbose::Bool,
    )
        push!(printlist, message)
        if verbose
            print(message)
        end
    end

    start_time = time()

    artificial_subpaths = generate_artificial_subpaths(data)
    if with_heuristic
        some_subpaths = construct_heuristic_subpaths(data, T_range, B_range)
        for (key, artificial_subpath_list) in pairs(artificial_subpaths)
            if !(key in keys(some_subpaths))
                some_subpaths[key] = []
            end
            append!(some_subpaths[key], artificial_subpath_list)
        end
    else
        some_subpaths = deepcopy(artificial_subpaths)
    end
    some_subpaths = deepcopy(artificial_subpaths)
    mp_results = Dict()
    params = Dict()
    params["number_of_subpaths"] = [sum(length(v) for v in values(some_subpaths))]
    params["objective"] = Float64[]
    params["κ"] = Dict{Int, Float64}[]
    params["λ"] = []
    params["μ"] = Dict{Int, Float64}[]
    params["ν"] = Vector{Float64}[]
    params["lp_relaxation_time_taken"] = Float64[]
    params["sp_total_time_taken"] = Float64[]
    params["number_of_current_subpaths"] = Int[]
    params["smallest_reduced_costs"] = Float64[]
    
    printlist = String[]
    counter = 0
    converged = false

    add_message!(
        printlist,
        @sprintf(
            """
            Starting column generation.
            # customers:                %2d
            # depots:                   %2d
            # charging stations:        %2d
            # vehicles:                 %2d

            """,
            data["n_customers"],
            data["n_depots"],
            data["n_charging"],
            data["n_vehicles"],
        ),
        verbose,
    )
    add_message!(
        printlist,
        @sprintf("              |  Objective | # subpaths | Time (LP) | Time (SP) | # new subpaths | smallest r.c.\n"),
        verbose,
    )
    while !converged
        counter += 1
        subpath_costs = compute_subpath_costs(
            data, 
            some_subpaths,
            ;
            with_charging_cost = with_charging_cost,
            time_windows = time_windows,
        )
        subpath_service = compute_subpath_service(
            data, 
            some_subpaths,
        )
        mp_results, mp_params = @suppress subpath_formulation(
            data, 
            some_subpaths,
            subpath_costs, 
            subpath_service,
            T_range,
            B_range,
            ;
            integral = false,
            charging_in_subpath = charging_in_subpath,
        )

        push!(params["objective"], mp_results["objective"])
        push!(params["κ"], mp_results["κ"])
        push!(params["λ"], mp_results["λ"])
        push!(params["μ"], mp_results["μ"])
        push!(params["ν"], mp_results["ν"])
        push!(params["lp_relaxation_time_taken"], mp_params["time_taken"])

        # generate subpaths
        generate_subpaths_result = @timed generate_subpaths_withcharge(
            G, data, T_range, B_range,
            mp_results["κ"], mp_results["λ"], mp_results["μ"], mp_results["ν"],
            ;
            charging_in_subpath = charging_in_subpath,
            charge_to_full_only = charge_to_full_only,
            time_windows = time_windows,
            with_charging_cost = with_charging_cost,
        )
        (current_subpaths, smallest_reduced_costs, _) = generate_subpaths_result.value

        push!(
            params["sp_total_time_taken"],
            round(generate_subpaths_result.time, digits=3)
        )
        if length(current_subpaths) == 0
            push!(params["number_of_current_subpaths"], 0)
            converged = true
        else
            push!(
                params["number_of_current_subpaths"],
                sum(length(v) for v in values(current_subpaths))
            )
        end
        push!(
            params["smallest_reduced_costs"],
            minimum(values(smallest_reduced_costs)),
        )
        for key in keys(current_subpaths)
            if !(key in keys(some_subpaths))
                some_subpaths[key] = current_subpaths[key]
            else
                for s_new in current_subpaths[key]
                    if !any(isequal(s_new, s) for s in some_subpaths[key])
                        push!(some_subpaths[key], s_new)
                    end
                end
            end
        end
        push!(
            params["number_of_subpaths"], 
            sum(length(v) for v in values(some_subpaths))
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %9.3f | %9.3f | %14d | %.4e\n", 
                counter,
                params["objective"][counter],
                params["number_of_subpaths"][counter],
                params["lp_relaxation_time_taken"][counter],
                params["sp_total_time_taken"][counter],
                params["number_of_current_subpaths"][counter],
                params["smallest_reduced_costs"][counter],
            ),
            verbose,
        )
        if length(params["number_of_subpaths"]) > 1
            if params["number_of_subpaths"][end-1] == params["number_of_subpaths"][end]
                break
            end
        end
    end
    results = Dict(
        "objective" => mp_results["objective"],
        "z" => mp_results["z"],
    )
    params["counter"] = counter
    params["converged"] = converged
    end_time = time()
    time_taken = round(end_time - start_time, digits = 3)
    params["time_taken"] = time_taken

    for message in [
        @sprintf("\n")
        @sprintf("Objective: \t\t%.4e\n", mp_results["objective"])
        @sprintf("Total time (LP): \t%10.3f s\n", sum(params["lp_relaxation_time_taken"]))
        @sprintf("Total time (SP): \t%10.3f s\n", sum(params["sp_total_time_taken"]))
        @sprintf("Total time: \t\t%10.3f s\n", time_taken)
    ]
        add_message!(printlist, message, verbose)
    end
    return results, params, printlist, some_subpaths
end

function subpath_formulation_column_generation_integrated_from_paths(
    G,
    data, 
    T_range,
    B_range,
    ;
    charge_to_full_only::Bool = false,
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
    with_heuristic::Bool = true,
    verbose::Bool = true,
    time_limit::Float64 = Inf,
)
    function add_message!(
        printlist::Vector, 
        message::String, 
        verbose::Bool,
    )
        push!(printlist, message)
        if verbose
            print(message)
        end
    end

    start_time = time()

    artificial_subpaths = generate_artificial_subpaths(data)
    if with_heuristic
        some_subpaths = construct_heuristic_subpaths(data, T_range, B_range)
        for (key, artificial_subpath_list) in pairs(artificial_subpaths)
            if !(key in keys(some_subpaths))
                some_subpaths[key] = []
            end
            append!(some_subpaths[key], artificial_subpath_list)
        end
    else
        some_subpaths = deepcopy(artificial_subpaths)
    end
    subpath_costs = compute_subpath_costs(
        data, 
        some_subpaths,
        ;
        with_charging_cost = with_charging_cost,
        time_windows = time_windows,
    )
    subpath_service = compute_subpath_service(
        data, 
        some_subpaths,
    )
    mp_results = Dict()
    params = Dict()
    params["number_of_subpaths"] = [sum(length(v) for v in values(some_subpaths))]
    params["objective"] = Float64[]
    params["κ"] = Dict{Int, Float64}[]
    params["μ"] = Dict{Int, Float64}[]
    params["ν"] = Vector{Float64}[]
    params["lp_relaxation_solution_time_taken"] = Float64[]
    params["sp_total_time_taken"] = Float64[]
    params["lp_relaxation_constraint_time_taken"] = Float64[]
    params["number_of_current_subpaths"] = Int[]
    params["number_of_new_charging_states"] = Int[]
    params["number_of_charging_states"] = Int[]
    charging_states = Set()

    printlist = String[]
    counter = 0
    converged = false

    add_message!(
        printlist,
        @sprintf(
            """
            Starting column generation.
            # customers:                %2d
            # depots:                   %2d
            # charging stations:        %2d
            # vehicles:                 %2d

            """,
            data["n_customers"],
            data["n_depots"],
            data["n_charging"],
            data["n_vehicles"],
        ),
        verbose,
    )
    add_message!(
        printlist,
        @sprintf("              |  Objective | # subpaths | Time (LP) | Time (SP) | Time (LP) | # subpaths \n"),
        verbose,
    )
    add_message!(
        printlist,
        @sprintf("              |            |            |           |           |    (cons) |      (new) \n"),
        verbose,
    )

    mp_model = @suppress Model(Gurobi.Optimizer)
    JuMP.set_string_names_on_creation(mp_model, false)
    z = Dict{Tuple{Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, Int}, VariableRef}(
        (key, p) => @variable(mp_model, lower_bound = 0)
        for key in keys(some_subpaths)
            for p in 1:length(some_subpaths[key])
    )
    @constraint(
        mp_model,
        κ[i in data["N_depots"]],
        sum(
            sum(
                z[((i,0,data["B"]),state2),p]
                for p in 1:length(some_subpaths[((i,0,data["B"]),state2)])
            )        
            for (state1, state2) in keys(some_subpaths)
                if state1[1] == i && state1[2] == 0 && state1[3] == data["B"]
        )
        == data["v_start"][findfirst(x -> (x == i), data["N_depots"])]
    )
    
    flow_conservation_exprs_out = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
    flow_conservation_exprs_in = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
    flow_conservation_constrs = Dict{Tuple{Int, Float64, Float64}, ConstraintRef}()

    if with_heuristic
        for key in keys(some_subpaths)
            for (ind, s_new) in enumerate(some_subpaths[key])
                # modify flow conservation constraints
                if key[1][1] in data["N_charging"]
                    push!(charging_states, key[1])
                    if !(key[1] in keys(flow_conservation_exprs_out))
                        flow_conservation_exprs_out[key[1]] = @expression(mp_model, 0)
                    end
                    if !(key[1] in keys(flow_conservation_exprs_in))
                        flow_conservation_exprs_in[key[1]] = @expression(mp_model, 0)
                    end
                    add_to_expression!(flow_conservation_exprs_out[key[1]], z[key, ind])
                end
                if key[2][1] in data["N_charging"]
                    push!(charging_states, key[2])
                    if !(key[2] in keys(flow_conservation_exprs_out))
                        flow_conservation_exprs_out[key[2]] = @expression(mp_model, 0)
                    end
                    if !(key[2] in keys(flow_conservation_exprs_in))
                        flow_conservation_exprs_in[key[2]] = @expression(mp_model, 0)
                    end
                    add_to_expression!(flow_conservation_exprs_in[key[2]], z[key, ind])
                end
            end
        end
        
        for state in charging_states
            flow_conservation_constrs[state] = @constraint(
                mp_model,
                flow_conservation_exprs_out[state]
                == flow_conservation_exprs_in[state]
            )
        end
    end


    @constraint(
        mp_model,
        μ[n2 in data["N_depots"]],
        sum(
            sum(
                z[(state1, state2),p]
                for p in 1:length(some_subpaths[(state1, state2)])
            )
            for (state1, state2) in keys(some_subpaths)
                if state2[1] == n2
        ) ≥ data["v_end"][n2]
    )
    @constraint(
        mp_model,
        ν[j in data["N_pickups"]],
        sum(
            sum(
                subpath_service[((state1, state2),j)][p] * z[(state1, state2),p]
                for p in 1:length(some_subpaths[(state1, state2)])
            )
            for (state1, state2) in keys(some_subpaths)
        ) == 1
    )
    @expression(
        mp_model,
        subpath_costs_expr,
        sum(
            sum(
                subpath_costs[state_pair][p] * z[state_pair,p]
                for p in 1:length(some_subpaths[state_pair])
            )
            for state_pair in keys(some_subpaths)
        )
    )
    @objective(mp_model, Min, subpath_costs_expr)

    while (
        !converged
        && time_limit ≥ (time() - start_time)
    )
        counter += 1
        mp_solution_start_time = time()
        @suppress optimize!(mp_model)
        mp_solution_end_time = time()
        mp_results = Dict(
            "model" => mp_model,
            "objective" => objective_value(mp_model),
            "z" => Dict(
                (key, p) => value.(z[(key, p)])
                for (key, p) in keys(z)
            ),
            "κ" => Dict(zip(data["N_depots"], dual.(mp_model[:κ]).data)),
            "μ" => Dict(zip(data["N_depots"], dual.(mp_model[:μ]).data)),
            "ν" => dual.(mp_model[:ν]).data,
            "solution_time_taken" => round(mp_solution_end_time - mp_solution_start_time, digits = 3),
        )
        push!(params["objective"], mp_results["objective"])
        push!(params["κ"], mp_results["κ"])
        push!(params["μ"], mp_results["μ"])
        push!(params["ν"], mp_results["ν"])
        push!(params["lp_relaxation_solution_time_taken"], mp_results["solution_time_taken"])

        if time_windows
            generate_subpaths_result = @timed generate_subpaths_withcharge_from_paths(
                G, data, T_range, B_range,
                mp_results["κ"],
                mp_results["μ"], 
                mp_results["ν"],
                ;
                charging_in_subpath = true,
                charge_to_full_only = charge_to_full_only,
                time_windows = time_windows,
                with_charging_cost = with_charging_cost,
            )
        else
            generate_subpaths_result = @timed generate_subpaths_withcharge_from_paths_notimewindows_V3(
                G, data, T_range, B_range,
                mp_results["κ"],
                mp_results["μ"], 
                mp_results["ν"],
                ;
                charge_to_full_only = charge_to_full_only,
                with_charging_cost = with_charging_cost,
            )
        end
        (current_subpaths, _, _) = generate_subpaths_result.value


        push!(
            params["sp_total_time_taken"],
            round(generate_subpaths_result.time, digits=3)
        )
        if length(current_subpaths) == 0
            push!(params["number_of_current_subpaths"], 0)
            converged = true
        else
            push!(
                params["number_of_current_subpaths"],
                sum(length(v) for v in values(current_subpaths))
            )
        end

        new_charging_states = Set()
        mp_constraint_start_time = time()
        for key in keys(current_subpaths)
            if !(key in keys(some_subpaths))
                some_subpaths[key] = []
                subpath_costs[key] = []
                for i in 1:data["n_customers"]
                    subpath_service[(key, i)] = []
                end
                count = 0
            else
                count = length(some_subpaths[key])
            end
            for s_new in current_subpaths[key]
                if key in keys(some_subpaths)
                    add = !any(isequal(s_new, s) for s in some_subpaths[key])
                else
                    add = true
                end
                if add
                    # 1: include in some_subpaths
                    push!(some_subpaths[key], s_new)
                    # 2: add subpath cost
                    push!(
                        subpath_costs[key], 
                        compute_subpath_cost(
                            data, s_new,
                            ;
                            with_charging_cost = with_charging_cost,
                            time_windows = time_windows,
                        )
                    )
                    # 3: add subpath service
                    for i in 1:data["n_customers"]
                        push!(subpath_service[(key, i)], s_new.served[i])
                    end
                    # 4: create variable
                    count += 1
                    z[(key, count)] = @variable(mp_model, lower_bound = 0)
                    # 5: modify constraints starting from depot, ending at depot, and flow conservation
                    if key[1][1] in data["N_depots"] && key[1][2] == 0.0 && key[1][3] == data["B"]
                        set_normalized_coefficient(κ[key[1][1]], z[key,count], 1)
                    elseif key[1][1] in data["N_charging"]
                        push!(new_charging_states, key[1])
                        if !(key[1] in keys(flow_conservation_exprs_out))
                            flow_conservation_exprs_out[key[1]] = @expression(mp_model, 0)
                        end
                        if !(key[1] in keys(flow_conservation_exprs_in))
                            flow_conservation_exprs_in[key[1]] = @expression(mp_model, 0)
                        end
                        add_to_expression!(flow_conservation_exprs_out[key[1]], z[key, count])
                    end
                    if key[2][1] in data["N_depots"]
                        set_normalized_coefficient(μ[key[2][1]], z[key,count], 1)
                    elseif key[2][1] in data["N_charging"]
                        push!(new_charging_states, key[2])
                        if !(key[2] in keys(flow_conservation_exprs_out))
                            flow_conservation_exprs_out[key[2]] = @expression(mp_model, 0)
                        end
                        if !(key[2] in keys(flow_conservation_exprs_in))
                            flow_conservation_exprs_in[key[2]] = @expression(mp_model, 0)
                        end
                        add_to_expression!(flow_conservation_exprs_in[key[2]], z[key, count])
                    end
                    # 6: modify customer service constraints
                    for l in data["N_pickups"]
                        if subpath_service[(key, l)][count] == 1
                            set_normalized_coefficient(ν[l], z[key, count], 1)
                        end
                    end
                    # 7: modify objective
                    set_objective_coefficient(mp_model, z[key, count], subpath_costs[key][count])
                end
            end
        end
        
        for state in new_charging_states
            if state in charging_states
                con = pop!(flow_conservation_constrs, state)
                delete(mp_model, con)
            end
            flow_conservation_constrs[state] = @constraint(
                mp_model,
                flow_conservation_exprs_out[state]
                == flow_conservation_exprs_in[state]
            )
        end
        union!(charging_states, new_charging_states)
        mp_constraint_end_time = time()

        push!(
            params["number_of_new_charging_states"],
            length(new_charging_states)
        )
        push!(
            params["number_of_charging_states"],
            length(charging_states)
        )
        push!(
            params["number_of_subpaths"], 
            sum(length(v) for v in values(some_subpaths))
        )
        push!(
            params["lp_relaxation_constraint_time_taken"],
            round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %9.3f | %9.3f | %9.3f | %10d \n", 
                counter,
                params["objective"][counter],
                params["number_of_subpaths"][counter],
                params["lp_relaxation_solution_time_taken"][counter],
                params["sp_total_time_taken"][counter],
                params["lp_relaxation_constraint_time_taken"][counter],
                params["number_of_current_subpaths"][counter],
            ),
            verbose,
        )
    end

    results = Dict(
        "objective" => mp_results["objective"],
        "z" => mp_results["z"],
        "κ" => mp_results["κ"],
        "μ" => mp_results["μ"],
        "ν" => mp_results["ν"],
    )
    params["converged"] = converged
    params["counter"] = counter
    end_time = time() 
    time_taken = round(end_time - start_time, digits = 3)
    params["time_taken"] = time_taken
    params["time_limit_reached"] = (time_taken > time_limit)
    params["lp_relaxation_time_taken"] = params["lp_relaxation_constraint_time_taken"] .+ params["lp_relaxation_solution_time_taken"]

    for message in [
        @sprintf("\n")
        @sprintf("Objective: \t\t\t%.4e\n", mp_results["objective"])
        @sprintf("Total time (LP): \t\t%10.3f s\n", sum(params["lp_relaxation_solution_time_taken"]))
        @sprintf("Total time (SP): \t\t%10.3f s\n", sum(params["sp_total_time_taken"]))
        @sprintf("Total time (LP construction): \t%10.3f s\n", sum(params["lp_relaxation_constraint_time_taken"]))
        @sprintf("Total time: \t\t\t%10.3f s\n", time_taken)
    ]
        add_message!(printlist, message, verbose)
    end
    return results, params, printlist, some_subpaths

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