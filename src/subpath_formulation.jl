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
    for (v, (starting_node, current_node)) in enumerate(zip(start_depots, end_depots))
        starting_time = 0.0
        starting_charge = data["B"]
        current_time = 0.0
        current_charge = data["B"]
        key = (
            (starting_node, starting_time, starting_charge),  
            (current_node, current_time, current_charge)
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
            current_node = current_node,
            arcs = [],
            current_time = current_time,
            current_charge = current_charge,
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

function compute_subpath_modified_cost(
    s::Subpath,
    data,
    κ,
    μ,
    ν,
)
    reduced_cost = compute_subpath_cost(data, s)

    for j in findall(s.served)
        reduced_cost = reduced_cost - ν[j]
    end

    if s.starting_node in data["N_depots"]
        if s.starting_time == 0.0 && s.starting_charge == data["B"]
            reduced_cost = reduced_cost - κ[s.starting_node]
        end
    end

    if s.current_node in data["N_depots"]
        reduced_cost = reduced_cost - μ[s.current_node]
    end

    return reduced_cost
end

function compute_subpath_cost(
    data,
    s::Subpath,
    M::Float64 = 1e8,
    ;
)
    if s.artificial 
        return M
    elseif length(s.arcs) == 0
        return 0
    end

    travel_cost = data["travel_cost_coeff"] * sum(
        data["c"][a...] for a in s.arcs
    )
    return travel_cost
end

function compute_subpath_costs(
    data,
    all_subpaths,
    M::Float64 = 1e8,
    ;
)
    subpath_costs = Dict(
        key => [
            compute_subpath_cost(data, s, M;)
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

function compute_charging_arc_cost(
    data,
    a::ChargingArc,
)
    return data["charge_cost_coeff"] * a.delta_time
end

function generate_base_subpaths(
    G, 
    data, 
    κ,
    μ,
    ν,
)
    function add_subpath_to_collection!(
        collection::SortedDict{
            Float64,
            SubpathWithCost,
        },
        k1::Float64,
        v1::SubpathWithCost,
        ;
    )
        added = true
        for (k2, v2) in collection
            if k2 ≤ k1
                if v2.cost ≤ v1.cost
                    added = false
                    break
                end
            elseif k1 ≤ k2
                if v1.cost ≤ v2.cost
                    pop!(collection, k2)
                end
            end
        end
        if added
            insert!(collection, k1, v1)
        end
        return added
    end

    function connect(v1::SubpathWithCost, v2::SubpathWithCost, data)
        if v1.current_node != v2.starting_node
            return
        end
        if any(v1.served .&& v2.served)
            return
        end
        new_current_time = v1.current_time + (v2.current_time - v2.starting_time)
        if new_current_time > data["T"]
            return
        end
        new_current_charge = v1.current_charge - v2.starting_charge + v2.current_charge
        if new_current_charge < 0.0
            return
        end
        v = SubpathWithCost(
            cost = v1.cost + v2.cost,
            n_customers = data["n_customers"],
            starting_node = v1.starting_node,
            starting_time = v1.starting_time, 
            starting_charge = v1.starting_charge,
            current_node = v2.current_node,
            current_time = new_current_time,
            current_charge = new_current_charge,
            arcs = vcat(v1.arcs, v2.arcs),
            served = v1.served .|| v2.served,
        )
        return v
    end

    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for j in data["N_customers"]
        for i in data["N_nodes"]
            modified_costs[i,j] -= ν[j]
        end
    end

    base_labels = Dict(
        start_node => Dict(
            current_node => SortedDict{Float64, SubpathWithCost}()
            for current_node in data["N_nodes"]
        )
        for start_node in data["N_nodes"]
    )
    for (start_node, current_node) in keys(data["A"])
        current_time = data["t"][start_node, current_node]
        current_charge = data["B"] - data["q"][start_node, current_node]
        served = falses(data["n_customers"])
        if current_node in data["N_customers"]
            served[current_node] = true
        end
        base_labels[start_node][current_node][current_time] = SubpathWithCost(
            cost = modified_costs[start_node, current_node],
            n_customers = data["n_customers"],
            starting_node = start_node,
            starting_time = 0.0,
            starting_charge = data["B"],
            arcs = [(start_node, current_node)],
            current_node = current_node,
            current_time = current_time,
            current_charge = current_charge,
            served = served,
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
                ## TODO: find a better way to update collection with pairwise sum of two collections
                for (k1, v1) in pairs(base_labels[start_node][new_node])
                    for (k2, v2) in pairs(base_labels[new_node][end_node])
                        k = k1 + k2
                        v = connect(v1, v2, data)
                        if !isnothing(v)
                            add_subpath_to_collection!(
                                base_labels[start_node][end_node],
                                k, v,
                            )
                        end
                    end
                end
            end
        end
    end

    for start_node in data["N_depots"]
        for end_node in data["N_nodes"]
            for (k, v) in pairs(base_labels[start_node][end_node])
                v.cost = v.cost - κ[start_node]
            end
        end
    end
    for end_node in data["N_depots"]
        for start_node in data["N_nodes"]
            for (k, v) in pairs(base_labels[start_node][end_node])
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

function find_nondominated_paths(
    data,
    base_labels,
    κ,
    μ,
)
    function charge_to_specified_level(
        start_charge, end_charge, 
        start_time, charging_rate,
    )
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

    function add_path_withcost_to_collection!(
        collection::SortedDict{
            Tuple{Float64, Float64}, 
            PathWithCost,
            DoubleStateOrdering
        },
        k1::Tuple{Float64, Float64}, 
        v1::PathWithCost,
        ;
        verbose::Bool = false,
    )
        added = true
        for (k2, v2) in collection
            # check if v2 dominates v1
            if v2.cost ≤ v1.cost
                if k2[1] ≤ k1[1] && k2[2] ≥ k1[2]
                    added = false
                    if verbose
                        println("$(k1[1]), $(k1[2]), $(v1.cost) dominated by $(k2[1]), $(k2[2]), $(v2.cost)")
                    end
                    break
                end
            # check if v1 dominates v2
            elseif v1.cost ≤ v2.cost
                if k1[1] ≤ k2[1] && k1[2] ≥ k2[2]
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

    tso = TripleStateOrdering()
    dso = DoubleStateOrdering()

    full_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Float64, Float64}, 
                PathWithCost
            }(dso)
            for current_node in union(data["N_charging"], data["N_depots"])
        )
        for starting_node in data["N_depots"]
    )
    for depot in data["N_depots"]
        full_labels[depot][depot][(0.0, data["B"])] = PathWithCost(
            subpaths = [],
            charging_arcs = [],
            cost = 0.0,
            served = zeros(Int, data["n_customers"]),
        )
    end
    unexplored_states = SortedSet(
        tso,
        [
            (depot, 0.0, data["B"])
            for depot in data["N_depots"]
        ]
    )

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        for starting_node in data["N_depots"]
            if !((state[2], state[3]) in keys(full_labels[starting_node][state[1]]))
                continue
            end
            current_path = full_labels[starting_node][state[1]][(state[2], state[3])]
            for next_node in union(data["N_depots"], data["N_charging"])
                for s in values(base_labels[state[1]][next_node])
                    if (
                        s.current_node in data["N_depots"]
                        && s.current_time == 0.0
                        && s.current_charge == data["B"]
                    )
                        continue
                    end
                    if (
                        any(s.served + current_path.served .> 1)
                    )
                        continue
                    end
                    charge_required = data["B"] - s.current_charge
                    (
                        delta_time, delta_charge,
                        end_time, end_charge,
                    ) = charge_to_specified_level(
                        state[3], charge_required, 
                        state[2], data["μ"],
                    )
                    if end_time + s.current_time > data["T"]
                        continue
                    end
                    if end_charge - charge_required < 0
                        continue
                    end
                    
                    s_new = copy(s)
                    s_new.starting_time = end_time
                    s_new.starting_charge = end_charge
                    s_new.current_time = s.current_time + end_time
                    s_new.current_charge = end_charge - (data["B"] - s.current_charge)
                    
                    new_path = copy(current_path)
                    new_path.cost += s_new.cost
                    push!(new_path.subpaths, s_new)
                    if length(current_path.subpaths) > 0
                        a_new = ChargingArc(
                            starting_node = state[1],
                            starting_time = state[2],
                            starting_charge = state[3],
                            delta_time = delta_time,
                            delta_charge = delta_charge,
                            current_time = end_time,
                            current_charge = end_charge,
                        )
                        new_path.cost += compute_charging_arc_cost(data, a_new)
                        push!(new_path.charging_arcs, a_new)
                    end

                    key = (s_new.current_time, s_new.current_charge)

                    if add_path_withcost_to_collection!(
                        full_labels[starting_node][next_node],
                        key, new_path,
                    )
                        next_state = (next_node, key[1], key[2])
                        if (
                            next_node in data["N_charging"] 
                            && !(next_state in unexplored_states)
                        )
                            insert!(unexplored_states, next_state)
                        end
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
                    SubpathWithCost(
                        cost = - κ[depot] - μ[depot],
                        n_customers = data["n_customers"],
                        starting_node = depot,
                        starting_time = 0.0, 
                        starting_charge = data["B"],
                    )
                ],
                charging_arcs = [],
                cost = - κ[depot] - μ[depot],
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

function get_subpaths_charging_arcs_from_negative_paths(
    data, full_labels,
)
    generated_subpaths = Dict{
        Tuple{
            Tuple{Int, Float64, Float64}, 
            Tuple{Int, Float64, Float64}
        }, 
        Vector{Subpath},
    }()
    generated_charging_arcs = Dict{
        Tuple{
            Tuple{Int, Float64, Float64}, 
            Tuple{Int, Float64, Float64}
        }, 
        Vector{ChargingArc},
    }()
    for starting_node in data["N_depots"], end_node in data["N_depots"]
        for path in values(full_labels[starting_node][end_node])
            if path.cost ≥ -1e-6
                continue
            end
            for s in path.subpaths
                state_pair = (
                    (s.starting_node, s.starting_time, s.starting_charge),
                    (s.current_node, s.current_time, s.current_charge),
                )
                if state_pair in keys(generated_subpaths)
                    if !any(isequal(s, s1) for s1 in generated_subpaths[state_pair])
                        push!(generated_subpaths[state_pair], Subpath(s))
                    end
                else
                    generated_subpaths[state_pair] = [Subpath(s)]
                end
            end
            for a in path.charging_arcs
                state_pair = (
                    (a.starting_node, a.starting_time, a.starting_charge),
                    (a.starting_node, a.current_time, a.current_charge),
                )
                if state_pair in keys(generated_charging_arcs)
                    if !any(isequal(a, a1) for a1 in generated_charging_arcs[state_pair])
                        push!(generated_charging_arcs[state_pair], a)
                    end
                else
                    generated_charging_arcs[state_pair] = [a]
                end
            end
        end
    end
    return generated_subpaths, generated_charging_arcs
end

function subpath_formulation_column_generation_from_paths(
    G,
    data,
    T_range,
    B_range,
    ;
    charging_in_subpath::Bool = true,
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

    return (generated_subpaths, generated_charging_arcs)

end

function subpath_formulation_column_generation_integrated_from_paths(
    G,
    data, 
    ;
    time_windows::Bool = true,
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

    some_subpaths = generate_artificial_subpaths(data)
    subpath_costs = compute_subpath_costs(
        data, 
        some_subpaths,
    )
    subpath_service = compute_subpath_service(
        data, 
        some_subpaths,
    )
    some_charging_arcs = Dict{
        Tuple{
            Tuple{Int64, Float64, Float64}, 
            Tuple{Int64, Float64, Float64}
        }, 
        Vector{ChargingArc}
    }()
    charging_arc_costs = Dict{
        Tuple{
            Tuple{Int64, Float64, Float64}, 
            Tuple{Int64, Float64, Float64}
        }, 
        Vector{Float64}
    }()
    charging_states_a_s = Set()
    charging_states_s_a = Set()
    mp_results = Dict()
    params = Dict()
    params["number_of_subpaths"] = [sum(length(v) for v in values(some_subpaths))]
    params["number_of_charging_arcs"] = [0]
    params["objective"] = Float64[]
    params["κ"] = Dict{Int, Float64}[]
    params["μ"] = Dict{Int, Float64}[]
    params["ν"] = Vector{Float64}[]
    params["lp_relaxation_solution_time_taken"] = Float64[]
    params["sp_base_time_taken"] = Float64[]
    params["sp_full_time_taken"] = Float64[]
    params["sp_total_time_taken"] = Float64[]
    params["lp_relaxation_constraint_time_taken"] = Float64[]
    params["number_of_new_subpaths"] = Int[]
    params["number_of_new_charging_arcs"] = Int[]
    params["number_of_new_charging_states"] = Int[]
    params["number_of_charging_states"] = Int[]

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
        @sprintf("              |  Objective | # subpaths |  # c. arcs | Time (LP) | Time (SP) | Time (SP) | Time (LP) | # subpaths |  # c. arcs \n"),
        verbose,
    )
    add_message!(
        printlist,
        @sprintf("              |            |            |            |           |    (base) |    (full) |    (cons) |      (new) |      (new) \n"),
        verbose,
    )

    mp_model = @suppress Model(Gurobi.Optimizer)
    set_attribute(mp_model, "MIPGap", 1e-10)
    JuMP.set_string_names_on_creation(mp_model, false)
    z = Dict{
        Tuple{
            Tuple{
                Tuple{Int, Float64, Float64}, 
                Tuple{Int, Float64, Float64}
            }, 
            Int
        }, 
        VariableRef
    }(
        (key, p) => @variable(mp_model, lower_bound = 0)
        for key in keys(some_subpaths)
            for p in 1:length(some_subpaths[key])
    )
    w = Dict{
        Tuple{
            Tuple{
                Tuple{Int, Float64, Float64}, 
                Tuple{Int, Float64, Float64}
            }, 
            Int
        }, 
        VariableRef
    }()
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
    
    flow_conservation_exprs_s_out = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
    flow_conservation_exprs_s_in = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
    flow_conservation_exprs_a_out = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
    flow_conservation_exprs_a_in = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
    flow_conservation_constrs_a_s = Dict{Tuple{Int, Float64, Float64}, ConstraintRef}()
    flow_conservation_constrs_s_a = Dict{Tuple{Int, Float64, Float64}, ConstraintRef}()

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
        ν[j in data["N_customers"]],
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
    @expression(
        mp_model,
        charging_arc_costs_expr,
        0,
    )
    @objective(mp_model, Min, subpath_costs_expr + charging_arc_costs_expr)

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
            "w" => Dict(
                (key, p) => value.(w[(key, p)])
                for (key, p) in keys(w)
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

        base_labels_result = @timed generate_base_subpaths(
            G, 
            data,
            mp_results["κ"],
            mp_results["μ"],
            mp_results["ν"],
        )
        full_labels_result = @timed find_nondominated_paths(
            data, base_labels_result.value,
            mp_results["κ"],
            mp_results["μ"],
        )
        (generated_subpaths, generated_charging_arcs) = get_subpaths_charging_arcs_from_negative_paths(
            data,
            full_labels_result.value,
        )

        push!(
            params["sp_base_time_taken"],
            round(base_labels_result.time, digits=3)
        )
        push!(
            params["sp_full_time_taken"],
            round(full_labels_result.time, digits=3)
        )
        push!(
            params["sp_total_time_taken"],
            round(base_labels_result.time + full_labels_result.time, digits=3)
        )
        if length(generated_subpaths) == 0
            push!(params["number_of_new_subpaths"], 0)
            converged = true
        else
            push!(
                params["number_of_new_subpaths"],
                sum(length(v) for v in values(generated_subpaths))
            )
        end

        if length(generated_charging_arcs) == 0
            push!(params["number_of_new_charging_arcs"], 0)
        else
            push!(
                params["number_of_new_charging_arcs"],
                sum(length(v) for v in values(generated_charging_arcs))
            )
        end

        new_charging_states_a_s = Set()
        new_charging_states_s_a = Set()
        mp_constraint_start_time = time()
        for state_pair in keys(generated_subpaths)
            if !(state_pair in keys(some_subpaths))
                some_subpaths[state_pair] = []
                subpath_costs[state_pair] = []
                for i in 1:data["n_customers"]
                    subpath_service[(state_pair, i)] = []
                end
                count = 0
            else
                count = length(some_subpaths[state_pair])
            end
            for s_new in generated_subpaths[state_pair]
                if state_pair in keys(some_subpaths)
                    add = !any(isequal(s_new, s) for s in some_subpaths[state_pair])
                else
                    add = true
                end
                if add
                    # 1: include in some_subpaths
                    push!(some_subpaths[state_pair], s_new)
                    # 2: add subpath cost
                    push!(
                        subpath_costs[state_pair], 
                        compute_subpath_cost(data, s_new)
                    )
                    # 3: add subpath service
                    for i in 1:data["n_customers"]
                        push!(subpath_service[(state_pair, i)], s_new.served[i])
                    end
                    # 4: create variable
                    count += 1
                    z[(state_pair, count)] = @variable(mp_model, lower_bound = 0)
                    (state1, state2) = state_pair
                    # 5: modify constraints starting from depot, ending at depot, and flow conservation
                    if state1[1] in data["N_depots"] && state1[2] == 0.0 && state1[3] == data["B"]
                        set_normalized_coefficient(κ[state1[1]], z[state_pair,count], 1)
                    elseif state1[1] in data["N_charging"]
                        push!(new_charging_states_a_s, state1)
                        if !(state1 in keys(flow_conservation_exprs_s_out))
                            flow_conservation_exprs_s_out[state1] = @expression(mp_model, 0)
                        end
                        add_to_expression!(flow_conservation_exprs_s_out[state1], z[state_pair, count])
                    end
                    if state2[1] in data["N_depots"]
                        set_normalized_coefficient(μ[state2[1]], z[state_pair,count], 1)
                    elseif state2[1] in data["N_charging"]
                        push!(new_charging_states_s_a, state2)
                        if !(state2 in keys(flow_conservation_exprs_s_in))
                            flow_conservation_exprs_s_in[state2] = @expression(mp_model, 0)
                        end
                        add_to_expression!(flow_conservation_exprs_s_in[state2], z[state_pair, count])
                    end
                    # 6: modify customer service constraints
                    for l in data["N_customers"]
                        if subpath_service[(state_pair, l)][count] == 1
                            set_normalized_coefficient(ν[l], z[state_pair, count], 1)
                        end
                    end
                    # 7: modify objective
                    set_objective_coefficient(mp_model, z[state_pair, count], subpath_costs[state_pair][count])
                end
            end
        end

        for state_pair in keys(generated_charging_arcs)
            if !(state_pair in keys(some_charging_arcs))
                some_charging_arcs[state_pair] = []
                charging_arc_costs[state_pair] = []
                count = 0
            else
                count = length(some_charging_arcs[state_pair])
            end
            for a_new in generated_charging_arcs[state_pair]
                if state_pair in keys(some_charging_arcs)
                    add = !any(isequal(a_new, a) for a in some_charging_arcs[state_pair])
                else
                    add = true
                end
                if add
                    # 1: include in some_charging_arcs
                    push!(some_charging_arcs[state_pair], a_new)
                    # 2: add charging arc cost
                    push!(
                        charging_arc_costs[state_pair], 
                        compute_charging_arc_cost(data, a_new)
                    )
                    # 4: create variable
                    count += 1
                    w[(state_pair, count)] = @variable(mp_model, lower_bound = 0)
                    (state1, state2) = state_pair
                    # 5: modify constraints starting from depot, ending at depot, and flow conservation
                    if state1[1] in data["N_charging"]
                        push!(new_charging_states_s_a, state1)
                        if !(state1 in keys(flow_conservation_exprs_a_out))
                            flow_conservation_exprs_a_out[state1] = @expression(mp_model, 0)
                        end
                        add_to_expression!(flow_conservation_exprs_a_out[state1], w[state_pair, count])
                    end
                    if state2[1] in data["N_charging"]
                        push!(new_charging_states_a_s, state2)
                        if !(state2 in keys(flow_conservation_exprs_a_in))
                            flow_conservation_exprs_a_in[state2] = @expression(mp_model, 0)
                        end
                        add_to_expression!(flow_conservation_exprs_a_in[state2], w[state_pair, count])
                    end
                    # 7: modify objective
                    set_objective_coefficient(
                        mp_model, 
                        w[state_pair, count], 
                        charging_arc_costs[state_pair][count],
                    )
                end
            end
        end
        
        for state in new_charging_states_a_s
            if state in charging_states_a_s
                con = pop!(flow_conservation_constrs_a_s, state)
                delete(mp_model, con)
            end
            flow_conservation_constrs_a_s[state] = @constraint(
                mp_model,
                flow_conservation_exprs_s_out[state]
                == flow_conservation_exprs_a_in[state]
            )
        end
        for state in new_charging_states_s_a
            if state in charging_states_s_a
                con = pop!(flow_conservation_constrs_s_a, state)
                delete(mp_model, con)
            end
            flow_conservation_constrs_s_a[state] = @constraint(
                mp_model,
                flow_conservation_exprs_a_out[state]
                == flow_conservation_exprs_s_in[state]
            )
        end
        union!(charging_states_s_a, new_charging_states_s_a)
        union!(charging_states_a_s, new_charging_states_a_s)
        mp_constraint_end_time = time()

        push!(
            params["number_of_new_charging_states"],
            length(new_charging_states_s_a) + length(new_charging_states_a_s)
        )
        push!(
            params["number_of_charging_states"],
            length(charging_states_s_a) + length(charging_states_a_s)
        )
        push!(
            params["number_of_subpaths"], 
            sum(length(v) for v in values(some_subpaths))
        )
        push!(
            params["number_of_charging_arcs"], 
            sum(length(v) for v in values(some_charging_arcs))
        )
        push!(
            params["lp_relaxation_constraint_time_taken"],
            round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %10d | %9.3f | %9.3f | %9.3f | %9.3f | %10d | %10d \n", 
                counter,
                params["objective"][counter],
                params["number_of_subpaths"][counter],
                params["number_of_charging_arcs"][counter],
                params["lp_relaxation_solution_time_taken"][counter],
                params["sp_base_time_taken"][counter],
                params["sp_full_time_taken"][counter],
                params["lp_relaxation_constraint_time_taken"][counter],
                params["number_of_new_subpaths"][counter],
                params["number_of_new_charging_arcs"][counter],
            ),
            verbose,
        )
    end

    for key in keys(some_subpaths)
        for p in 1:length(some_subpaths[key])
            JuMP.set_integer(z[key,p])
        end
    end
    for key in keys(some_charging_arcs)
        for p in 1:length(some_charging_arcs[key])
            JuMP.set_integer(w[key,p])
        end
    end
    cgip_solution_start_time = time()
    @suppress optimize!(mp_model)
    cgip_solution_end_time = time()

    CGLP_results = Dict(
        "objective" => mp_results["objective"],
        "z" => mp_results["z"],
        "w" => mp_results["w"],
        "κ" => mp_results["κ"],
        "μ" => mp_results["μ"],
        "ν" => mp_results["ν"],
    )
    CGIP_results = Dict(
        "objective" => objective_value(mp_model),
        "z" => Dict(
            (key, p) => value.(z[(key, p)])
            for (key, p) in keys(z)
        ),
        "w" => Dict(
            (key, p) => value.(w[(key, p)])
            for (key, p) in keys(w)
        ),
    )
    params["CGIP_time_taken"] = round(cgip_solution_end_time - cgip_solution_start_time, digits = 3)
    params["converged"] = converged
    params["counter"] = counter
    end_time = time() 
    time_taken = round(end_time - start_time, digits = 3)
    params["time_taken"] = time_taken
    params["time_limit_reached"] = (time_taken > time_limit)
    params["lp_relaxation_time_taken"] = params["lp_relaxation_constraint_time_taken"] .+ params["lp_relaxation_solution_time_taken"]

    for message in [
        @sprintf("\n"),
        @sprintf("(CG) Objective:               %.4e\n", CGLP_results["objective"]),
        @sprintf("Total time (LP):              %10.3f s\n", sum(params["lp_relaxation_solution_time_taken"])),
        @sprintf("Total time (SP):              %10.3f s\n", sum(params["sp_total_time_taken"])),
        @sprintf("Total time (LP construction): %10.3f s\n", sum(params["lp_relaxation_constraint_time_taken"])),
        @sprintf("Total time:                   %10.3f s\n", time_taken),
        @sprintf("\n"),
        @sprintf("(CGIP) Objective:             %.4e\n", CGIP_results["objective"]),
        @sprintf("(CGIP) Total time:            %10.3f s\n", params["CGIP_time_taken"]),
    ]
        add_message!(printlist, message, verbose)
    end
    return CGLP_results, CGIP_results, params, printlist, some_subpaths, some_charging_arcs
end

function construct_paths_from_subpath_solution(
    results, 
    data, 
    subpaths::Dict{
        Tuple{
            Tuple{Int64, Float64, Float64}, 
            Tuple{Int64, Float64, Float64}
        }, 
        Vector{Subpath}
    },
    charging_arcs::Dict{
        Tuple{
            Tuple{Int64, Float64, Float64}, 
            Tuple{Int64, Float64, Float64}
        }, 
        Vector{ChargingArc}
    },
    ;
)
    results_subpaths = []
    for key in keys(subpaths)
        for p in 1:length(subpaths[key])
            val = results["z"][key,p]
            if val != 0
                push!(results_subpaths, (val, subpaths[key][p]))
            end
        end
    end
    results_charging_arcs = []
    for key in keys(charging_arcs)
        for p in 1:length(charging_arcs[key])
            val = results["w"][key,p]
            if val != 0
                push!(results_charging_arcs, (val, charging_arcs[key][p]))
            end
        end
    end

    all_paths = []

    while sum(x[1] for x in results_subpaths) + sum(x[1] for x in results_charging_arcs) > 0
        pathlist = []
        current_states = [
            (depot, 0.0, data["B"]) for depot in data["N_depots"]
        ]
        while true
            i = findfirst(
                s -> (
                    (s[2].starting_node, s[2].starting_time, s[2].starting_charge) 
                    in current_states
                    && s[1] > 0
                ),
                results_subpaths
            )
            (val, s) = results_subpaths[i]
            push!(pathlist, (val, i, s))
            if s.current_node in data["N_depots"]
                break
            end
            current_states = [(s.current_node, s.current_time, s.current_charge)]
            i = findfirst(
                s -> (
                    (s[2].starting_node, s[2].starting_time, s[2].starting_charge) 
                    in current_states
                    && s[1] > 0
                ),
                results_charging_arcs
            )
            (val, a) = results_charging_arcs[i]
            push!(pathlist, (val, i, a))
            current_states = [(a.starting_node, a.current_time, a.current_charge)]
        end
        minval = minimum(x[1] for x in pathlist)
        for (val, i, s) in pathlist
            if typeof(s) == Subpath
                newval = results_subpaths[i][1] - minval
                if abs.(newval) < 1e-6
                    newval = 0.0
                end
                results_subpaths[i] = (
                    newval,
                    results_subpaths[i][2],
                )
            else
                newval = results_charging_arcs[i][1] - minval
                if abs.(newval) < 1e-6
                    newval = 0.0
                end
                results_charging_arcs[i] = (
                    newval,
                    results_charging_arcs[i][2],
                )
            end
        end
        path = Path(
            subpaths = [x[3] for (i, x) in enumerate(pathlist) if i % 2 == 1],
            charging_arcs = [x[3] for (i, x) in enumerate(pathlist) if i % 2 == 0],
        )
        push!(all_paths, (minval, path))
    end
    
    return all_paths
end

function subpath_results_printout(
    results,
    params,
    data,
    subpaths::Dict{
        Tuple{
            Tuple{Int64, Float64, Float64}, 
            Tuple{Int64, Float64, Float64}
        }, 
        Vector{Subpath}
    },
    charging_arcs::Dict{
        Tuple{
            Tuple{Int64, Float64, Float64}, 
            Tuple{Int64, Float64, Float64}
        }, 
        Vector{ChargingArc}
    },
    ;
)

    function print_subpath(
        s::Subpath,
        data,
    )
        subpath_state_list = [(s.starting_node, s.starting_time, s.starting_charge)]
        for arc in s.arcs
            prev_time = subpath_state_list[end][2]
            prev_charge = subpath_state_list[end][3]
            push!(subpath_state_list, (arc[2], prev_time + data["t"][arc...], prev_charge - data["q"][arc...]))
        end
        for (state1, state2) in zip(subpath_state_list[1:end-2], subpath_state_list[2:end-1])
            @printf(
                "               %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) |           | \n", 
                data["node_labels"][state1[1]], state1[2], state1[3],
                data["node_labels"][state2[1]], state2[2], state2[3],
            )
        end
        if length(subpath_state_list) > 1
            (state1, state2) = subpath_state_list[end-1:end]
        else
            (state1, state2) = (
                (s.starting_node, s.starting_time, s.starting_charge),
                (s.current_node, s.current_time, s.current_charge),
            )
        end
        @printf(
            "               %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | +  %6.1f | -    %6.1f \n", 
            data["node_labels"][state1[1]], state1[2], state1[3],
            data["node_labels"][state2[1]], state2[2], state2[3],
            s.current_time - s.starting_time,
            abs(s.current_charge - s.starting_charge),
        )
    end
    function print_charging_arc(
        a::ChargingArc,
        data,
    )
        @printf(
            "               %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | +  %6.1f | +    %6.1f \n", 
            data["node_labels"][a.starting_node],
            a.starting_time,
            a.starting_charge,
            data["node_labels"][a.starting_node],
            a.current_time,
            a.current_charge,
            a.delta_time,
            a.delta_charge,
        )
    end

    @printf("Objective:                        %.4e\n", results["objective"])
    @printf("Total time (LP):                  %10.3f s\n", sum(params["lp_relaxation_solution_time_taken"]))
    @printf("Total time (SP):                  %10.3f s\n", sum(params["sp_total_time_taken"]))
    @printf("Total time (LP construction):     %10.3f s\n", sum(params["lp_relaxation_constraint_time_taken"]))
    @printf("Total time:                       %10.3f s\n", params["time_taken"])

    println("Vehicles             From      time   charge             To      time   charge  |   sp_time |   sp_charge ")
    println("----------------------------------------------------------------------------------------------------------")

    all_paths = construct_paths_from_subpath_solution(
        results, data, subpaths, charging_arcs,
    )

    for (i, (val, path)) in enumerate(all_paths)
        println("----------------------------------------------------------------------------------------------------------")
        @printf("Vehicle %3d\n", i)
        @printf("Weight:   %5.3f\n", val)
        if length(path.subpaths) == 0
            continue
        end
        for (s, a) in zip(path.subpaths, path.charging_arcs)
            print_subpath(s, data)
            println("             -------------------------------------------------------------------|-----------|-------------")
            print_charging_arc(a, data)
            println("             -------------------------------------------------------------------|-----------|-------------")
        end
        print_subpath(path.subpaths[end], data)
    end
    println("----------------------------------------------------------------------------------------------------------")
end