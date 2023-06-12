using JuMP
using Gurobi
using Suppressor
using Printf

using Plots
using ColorSchemes

include("utils.jl")
include("desaulniers_benchmark.jl")
include("subpath_stitching.jl")

function generate_artificial_subpaths(data)
    artificial_subpaths = Dict{
        Tuple{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}},
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
            data["N_depots"], 
            outer = Int(ceil((data["n_vehicles"] - sum(values(data["v_end"]))) / data["n_depots"]))
        )
    )
    end_depots = end_depots[1:data["n_vehicles"]]

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
        served = zeros(Int, data["n_customers"])
        for i in 1:length(served)
            if mod1(i, data["n_vehicles"]) == v
                served[i] = 1
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

function add_subpath_to_generated_subpaths!(
    generated_subpaths::Dict{Tuple{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}, Vector{Subpath}},
    subpath::Subpath,
)
    state_pair = (
        (subpath.starting_node, subpath.starting_time, subpath.starting_charge),
        (subpath.current_node, subpath.current_time, subpath.current_charge)
    )
    if state_pair in keys(generated_subpaths)
        if !any(isequal(subpath, s) for s in generated_subpaths[state_pair])
            push!(generated_subpaths[state_pair], subpath)
        end
    else
        generated_subpaths[state_pair] = [subpath]
    end
    return
end

function add_charging_arc_to_generated_charging_arcs!(
    generated_charging_arcs::Dict{Tuple{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}, Vector{ChargingArc}},
    charging_arc::ChargingArc,
)
    state_pair = (
        (charging_arc.starting_node, charging_arc.starting_time, charging_arc.starting_charge),
        (charging_arc.starting_node, charging_arc.current_time, charging_arc.current_charge)
    )
    if state_pair in keys(generated_charging_arcs)
        if !any(isequal(charging_arc, a) for a in generated_charging_arcs[state_pair])
            push!(generated_charging_arcs[state_pair], charging_arc)
        end
    else
        generated_charging_arcs[state_pair] = [charging_arc]
    end
    return
end

function get_subpaths_charging_arcs_from_negative_path_labels(
    data, 
    full_labels = Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PathLabel,
        Base.Order.ForwardOrdering
    }}},
)
    generated_subpaths = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Subpath},
    }()
    generated_charging_arcs = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{ChargingArc},
    }()
    for starting_node in data["N_depots"], end_node in data["N_depots"]
        for path in values(full_labels[starting_node][end_node])
            current_time, current_charge = (0.0, data["B"])
            prev_time, prev_charge = current_time, current_charge
            s_labels = copy(path.subpath_labels)
            deltas = copy(path.charging_actions)
            while true
                s_label = popfirst!(s_labels)
                prev_time = current_time
                prev_charge = current_charge
                current_time = current_time + s_label.time_taken
                current_charge = current_charge - s_label.charge_taken
                s = Subpath(
                    n_customers = data["n_customers"],
                    starting_node = s_label.nodes[1],
                    starting_time = prev_time,
                    starting_charge = prev_charge,
                    current_node = s_label.nodes[end],
                    arcs = collect(zip(s_label.nodes[1:end-1], s_label.nodes[2:end])),
                    current_time = current_time,
                    current_charge = current_charge,
                    served = s_label.served,
                )
                add_subpath_to_generated_subpaths!(generated_subpaths, s)
                if length(deltas) == 0 
                    break
                end
                delta = popfirst!(deltas)
                prev_time = current_time
                prev_charge = current_charge
                current_time = current_time + delta
                current_charge = current_charge + delta
                a = ChargingArc(
                    starting_node = s_label.nodes[end], 
                    starting_time = prev_time, 
                    starting_charge = prev_charge, 
                    delta = delta,
                    current_time = current_time, 
                    current_charge = current_charge,
                )
                add_charging_arc_to_generated_charging_arcs!(generated_charging_arcs, a)
            end
        end
    end
    return generated_subpaths, generated_charging_arcs
end

function get_subpaths_charging_arcs_from_negative_pure_path_labels(
    data, 
    pure_labels = Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PurePathLabel,
        Base.Order.ForwardOrdering
    }}},
)
    generated_subpaths = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Subpath},
    }()
    generated_charging_arcs = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{ChargingArc},
    }()
    for starting_node in data["N_depots"], end_node in data["N_depots"]
        for path in values(pure_labels[starting_node][end_node])
            states = Tuple{Int, Int, Int}[]
            current_subpath = Subpath(
                n_customers = data["n_customers"],
                starting_node = path.nodes[1],
                starting_time = 0.0, 
                starting_charge = data["B"],
            )
            i = path.nodes[1]
            for (j, e, s) in zip(path.nodes[2:end], path.excesses, path.slacks)
                current_subpath.current_node = j
                push!(current_subpath.arcs, (i, j))
                current_subpath.starting_time += (e + s)
                current_subpath.starting_charge += (e + s) 
                current_subpath.current_time += (data["t"][i,j] + e + s)
                current_subpath.current_charge += (- data["q"][i,j] + e + s)
                if j in data["N_charging"]
                    push!(
                        states, 
                        (current_subpath.starting_node, current_subpath.starting_time, current_subpath.starting_charge), 
                        (current_subpath.current_node, current_subpath.current_time, current_subpath.current_charge), 
                    )
                    add_subpath_to_generated_subpaths!(generated_subpaths, current_subpath)
                    current_subpath = Subpath(
                        n_customers = data["n_customers"],
                        starting_node = j,
                        starting_time = current_subpath.current_time, 
                        starting_charge = current_subpath.current_charge,
                    )
                elseif j in data["N_customers"]
                    current_subpath.served[j] += 1
                end
                i = j
            end
            push!(
                states, 
                (current_subpath.starting_node, current_subpath.starting_time, current_subpath.starting_charge), 
                (current_subpath.current_node, current_subpath.current_time, current_subpath.current_charge), 
            )
            add_subpath_to_generated_subpaths!(generated_subpaths, current_subpath)
            for i in 1:(length(states)÷2)-1
                add_charging_arc_to_generated_charging_arcs!(
                    generated_charging_arcs, 
                    ChargingArc(
                        states[2*i][1],
                        states[2*i][2],
                        states[2*i][3],
                        states[2*i+1][2] - states[2*i][2],
                        states[2*i+1][2],
                        states[2*i+1][3],
                    )
                )
            end
        end
    end
    return generated_subpaths, generated_charging_arcs
end

function subpath_formulation_column_generation_integrated_from_paths(
    G,
    data, 
    ;
    Env = nothing,
    method::String = "ours",
    time_windows::Bool = false,
    subpath_single_service::Bool = false,
    subpath_check_customers::Bool = false,
    path_single_service::Bool = false,
    path_check_customers::Bool = false,
    check_customers_accelerated::Bool = false,
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

    compute_minimum_time_to_nearest_depot!(data, G)
    compute_minimum_charge_to_nearest_depot_charging_station!(data, G)

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
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{ChargingArc}
    }()
    charging_arc_costs = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Int}
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
            # customers:                    %2d
            # depots:                       %2d
            # charging stations:            %2d
            # vehicles:                     %2d
            method:                         %s
            subpath_single_service:         %s
            subpath_check_customers:        %s
            path_single_service:            %s
            path_check_customers:           %s
            check_customers_accelerated:    %s

            """,
            data["n_customers"],
            data["n_depots"],
            data["n_charging"],
            data["n_vehicles"],
            method,
            subpath_single_service,
            subpath_check_customers,
            path_single_service,
            path_check_customers,
            check_customers_accelerated,
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

    if isnothing(Env)
        mp_model = @suppress Model(Gurobi.Optimizer)
    else
        mp_model = @suppress Model(() -> Gurobi.Optimizer(Env))
    end
    JuMP.set_attribute(mp_model, "MIPGapAbs", 1e-3)
    JuMP.set_string_names_on_creation(mp_model, false)
    z = Dict{
        Tuple{
            Tuple{
                Tuple{Int, Int, Int}, 
                Tuple{Int, Int, Int}
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
                Tuple{Int, Int, Int}, 
                Tuple{Int, Int, Int}
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
    
    flow_conservation_exprs_s_out = Dict{Tuple{Int, Int, Int}, AffExpr}()
    flow_conservation_exprs_s_in = Dict{Tuple{Int, Int, Int}, AffExpr}()
    flow_conservation_exprs_a_out = Dict{Tuple{Int, Int, Int}, AffExpr}()
    flow_conservation_exprs_a_in = Dict{Tuple{Int, Int, Int}, AffExpr}()
    flow_conservation_constrs_a_s = Dict{Tuple{Int, Int, Int}, ConstraintRef}()
    flow_conservation_constrs_s_a = Dict{Tuple{Int, Int, Int}, ConstraintRef}()

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

    checkpoint_reached = false

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

        if method == "ours"
            if check_customers_accelerated && !checkpoint_reached
                (negative_full_labels, negative_full_labels_count, base_labels_time, full_labels_time) = subproblem_iteration_ours(
                    G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                    ;
                    subpath_single_service = subpath_single_service,        
                    subpath_check_customers = false,
                    path_single_service = path_single_service,
                    path_check_customers = false,
                )
                if negative_full_labels_count == 0
                    checkpoint_reached = true
                    (negative_full_labels, _, base_labels_time_new, full_labels_time_new) = subproblem_iteration_ours(
                        G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                        ;
                        subpath_single_service = subpath_single_service,        
                        subpath_check_customers = subpath_check_customers,
                        path_single_service = path_single_service,
                        path_check_customers = path_check_customers,
                    )
                    base_labels_time += base_labels_time_new
                    full_labels_time += full_labels_time_new
                end
            elseif check_customers_accelerated && checkpoint_reached
                (negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
                    G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                    ;
                    subpath_single_service = subpath_single_service,        
                    subpath_check_customers = subpath_check_customers,
                    path_single_service = path_single_service,
                    path_check_customers = path_check_customers,
                )
            else
                (negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
                    G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                    ;
                    subpath_single_service = subpath_single_service,        
                    subpath_check_customers = subpath_check_customers,
                    path_single_service = path_single_service,
                    path_check_customers = path_check_customers,
                )
            end
            (generated_subpaths, generated_charging_arcs) = get_subpaths_charging_arcs_from_negative_path_labels(
                data, negative_full_labels,
            )
            push!(
                params["sp_base_time_taken"],
                round(base_labels_time, digits=3)
            )
            push!(
                params["sp_full_time_taken"],
                round(full_labels_time, digits=3)
            )
            push!(
                params["sp_total_time_taken"],
                round(base_labels_time + full_labels_time, digits=3)
            )
        elseif method == "benchmark"
            if check_customers_accelerated && !checkpoint_reached
                (negative_pure_path_labels, negative_pure_path_labels_count, pure_path_labels_time) = subproblem_iteration_benchmark(
                    G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                    ;
                    time_windows = time_windows,
                    path_single_service = true,
                    path_check_customers = false,
                )
                if negative_pure_path_labels_count == 0
                    checkpoint_reached = true
                    (negative_pure_path_labels, _, pure_path_labels_time_new) = subproblem_iteration_benchmark(
                        G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                        ;
                        time_windows = time_windows,
                        path_single_service = true,
                        path_check_customers = true,
                    )
                    pure_path_labels_time += pure_path_labels_time_new
                end
            elseif check_customers_accelerated && checkpoint_reached
                (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
                    G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                    ;
                    time_windows = time_windows,
                    path_single_service = true,
                    path_check_customers = true,
                )
            else
                (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
                    G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                    ;
                    time_windows = time_windows,
                    path_single_service = path_single_service,
                    path_check_customers = path_check_customers,
                )
            end
            (generated_subpaths, generated_charging_arcs) = get_subpaths_charging_arcs_from_negative_pure_path_labels(
                data, negative_pure_path_labels,
            )
            push!(
                params["sp_base_time_taken"],
                0.0
            )
            push!(
                params["sp_full_time_taken"],
                round(pure_path_labels_time, digits=3)
            )
            push!(
                params["sp_total_time_taken"],
                round(pure_path_labels_time, digits=3)
            )
        end

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
                        set_normalized_coefficient(ν[l], z[state_pair, count], s_new.served[l])
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
        if length(some_charging_arcs) == 0
            push!(params["number_of_charging_arcs"], 0)
        else
            push!(
                params["number_of_charging_arcs"],
                sum(length(v) for v in values(some_charging_arcs))
            )
        end
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
    params["lp_relaxation_time_taken_total"] = sum(params["lp_relaxation_time_taken"])
    params["sp_base_time_taken_total"] = sum(params["sp_base_time_taken"])
    params["sp_full_time_taken_total"] = sum(params["sp_full_time_taken"])
    params["sp_time_taken_total"] = params["sp_base_time_taken_total"] + params["sp_full_time_taken_total"]
    params["lp_relaxation_time_taken_mean"] = params["lp_relaxation_time_taken_total"] / length(params["lp_relaxation_time_taken"])
    params["sp_base_time_taken_mean"] = params["sp_base_time_taken_total"] / length(params["sp_base_time_taken"])
    params["sp_full_time_taken_mean"] = params["sp_full_time_taken_total"] / length(params["sp_full_time_taken"])
    params["sp_time_taken_mean"] = params["sp_base_time_taken_mean"] + params["sp_full_time_taken_mean"]

    params["LP_IP_gap"] = CGIP_results["objective"] / CGLP_results["objective"] - 1.0

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

function collect_subpath_solution_support(
    results, 
    subpaths::Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Subpath}
    },
    charging_arcs::Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{ChargingArc}
    },
    ;
)
    results_subpaths = Tuple{Float64, Subpath}[]
    for key in keys(subpaths)
        for p in 1:length(subpaths[key])
            val = results["z"][key,p]
            if val > 1e-5
                push!(results_subpaths, (val, subpaths[key][p]))
            end
        end
    end
    results_charging_arcs = Tuple{Float64, ChargingArc}[]
    for key in keys(charging_arcs)
        for p in 1:length(charging_arcs[key])
            val = results["w"][key,p]
            if val > 1e-5
                push!(results_charging_arcs, (val, charging_arcs[key][p]))
            end
        end
    end
    return results_subpaths, results_charging_arcs
end

function collect_subpath_solution_metrics!(
    results,
    data, 
    subpaths, 
    charging_arcs,
)

    results["subpaths"], results["charging_arcs"] = collect_subpath_solution_support(results, subpaths, charging_arcs)
    results["paths"] = construct_paths_from_subpath_solution(results, data, subpaths, charging_arcs)

    results["mean_subpath_length"] = sum(
        sum(s.served) + 1 for (val, s) in results["subpaths"]
    ) / length(results["subpaths"])
    results["weighted_mean_subpath_length"] = sum(
        val * (sum(s.served) + 1) for (val, s) in results["subpaths"]
    ) / sum(
        val for (val, _) in results["subpaths"]
    )
    results["mean_path_length"] = sum(
        sum(p.served) + 1 for (val, p) in results["paths"]
    ) / length(results["paths"])

    results["weighted_mean_path_length"] = sum(
        val * (sum(p.served) + 1) for (val, p) in results["paths"]
    ) / sum(
        val for (val, _) in results["paths"]
    )
    results["mean_ps_length"] = sum(
        length(p.subpaths) for (val, p) in results["paths"]
    ) / length(results["paths"])
    results["weighted_mean_ps_length"] = sum(
        val * length(p.subpaths) for (val, p) in results["paths"]
    ) / sum(
        val for (val, _) in results["paths"]
    )
    return results

end

function plot_subpath_solution(
    results,
    data,
    subpaths,
    charging_arcs,
)
    p = plot_instance(data)
    results["subpaths"], results["charging_arcs"] = collect_subpath_solution_support(results, subpaths, charging_arcs)
    results["paths"] = construct_paths_from_subpath_solution(results, data, subpaths, charging_arcs)

    n_paths = length(results["paths"]) 
    colors = get(ColorSchemes.tol_bright, collect(0:1:(n_paths-1)) ./(n_paths-1))

    all_plots = []
    for (i, (val, path)) in enumerate(results["paths"])
        p = plot_instance(data)
        arcs = vcat(collect(s.arcs for s in path.subpaths)...)
        for (j, arc) in enumerate(arcs)
            plot!(
                data["coords"][1,collect(arc)],
                data["coords"][2,collect(arc)],
                arrow = true,
                color = colors[i],
                alpha = 0.5,
                lw = 1,
            )
        end
        plot!(title = "Vehicle $i: $val")
        push!(all_plots, p)
    end

    P = plot(
        all_plots..., 
        layout = (n_paths, 1), 
        size = (500, 400 * n_paths), 
        legend = false,
        fmt = :png,
    )
    return P
end



function construct_paths_from_subpath_solution(
    results, 
    data, 
    subpaths::Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Subpath}
    },
    charging_arcs::Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{ChargingArc}
    },
    ;
)
    results_subpaths, results_charging_arcs = collect_subpath_solution_support(results, subpaths, charging_arcs)

    all_paths = Tuple{Float64, Path}[]

    while true
        remaining_vals = 0.0
        for x in results_subpaths
            remaining_vals += x[1]
        end
        for x in results_charging_arcs
            remaining_vals += x[1]
        end
        if remaining_vals < 1e-3
            break
        end
        pathlist = []
        current_states = [
            (depot, 0.0, data["B"]) for depot in data["N_depots"]
        ]
        while true
            i = findfirst(
                s -> (
                    (s[2].starting_node, s[2].starting_time, s[2].starting_charge) 
                    in current_states
                    && s[1] > 1e-4
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
                    && s[1] > 1e-4
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
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Subpath}
    },
    charging_arcs::Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
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
            a.delta,
            a.delta,
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