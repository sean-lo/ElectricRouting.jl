using JuMP
using Gurobi
using Suppressor
using Printf

include("utils.jl")
include("desaulniers_benchmark.jl")
include("subpath_stitching.jl")

function generate_artificial_paths(data)
    artificial_paths = Dict{
        Tuple{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}},
        Vector{Path},
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
        key = (
            (starting_node, starting_time, starting_charge),  
            (current_node, starting_time, starting_charge)
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
            arcs = [(starting_node, current_node)],
            current_time = starting_time,
            current_charge = starting_charge,
            served = served,
            artificial = true,
        )
        p = Path(
            subpaths = [s],
            charging_arcs = ChargingArc[],
            served = served,
        )
        if !(key in keys(artificial_paths))
            artificial_paths[key] = []
        end
        push!(artificial_paths[key], p)
    end
    return artificial_paths
end

function add_path_to_generated_paths!(
    generated_paths::Dict{Tuple{Tuple{Int, Int, Int}, Tuple{Int, Int, Int}}, Vector{Path}},
    p::Path,
)
    state_pair = (
        (p.subpaths[1].starting_node, p.subpaths[1].starting_time, p.subpaths[1].starting_charge),
        (p.subpaths[end].current_node, p.subpaths[end].current_time, p.subpaths[end].current_charge)
    )
    if state_pair in keys(generated_paths)
        if !any(isequal(p, p2) for p2 in generated_paths[state_pair])
            push!(generated_paths[state_pair], p)
        end
    else
        generated_paths[state_pair] = [p]
    end
    return
end

function get_paths_from_negative_path_labels(
    data,
    path_labels::Vector{PathLabel},
)
    generated_paths = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Path},
    }()
    for path_label in path_labels
        current_time, current_charge = (0.0, data["B"])
        prev_time, prev_charge = current_time, current_charge
        s_labels = copy(path_label.subpath_labels)
        deltas = copy(path_label.charging_actions)
        p = Path(
            subpaths = Subpath[],
            charging_arcs = ChargingArc[],
            served = zeros(Int, data["n_customers"]),
        )
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
            push!(p.subpaths, s)
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
            push!(p.charging_arcs, a)
        end
        p.served = sum(s.served for s in p.subpaths)
        add_path_to_generated_paths!(generated_paths, p)
    end
    return generated_paths
end

function get_paths_from_negative_pure_path_labels(
    data, 
    pure_path_labels::Vector{PurePathLabel},
)
    generated_paths = Dict{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Vector{Path},
    }()
    for path_label in pure_path_labels
        p = Path(
            subpaths = Subpath[],
            charging_arcs = ChargingArc[],
            served = zeros(Int, data["n_customers"]),
        )
        states = Tuple{Int, Int, Int}[]
        current_subpath = Subpath(
            n_customers = data["n_customers"],
            starting_node = path_label.nodes[1],
            starting_time = 0.0, 
            starting_charge = data["B"],
        )
        i = path_label.nodes[1]
        for (j, e, s) in zip(path_label.nodes[2:end], path_label.excesses, path_label.slacks)
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
                push!(
                    p.subpaths,
                    current_subpath,
                )
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
        push!(
            p.subpaths,
            current_subpath,
        )
        for i in 1:(length(states)÷2)-1
            push!(
                p.charging_arcs, 
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
        p.served = sum(s.served for s in p.subpaths)
        add_path_to_generated_paths!(generated_paths, p)
    end
    return generated_paths
end

function path_formulation_column_generation(
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
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(data["n_customers"]))),
    ngroute_neighborhood_charging_depots_size::String = "small",
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
    if ngroute
        compute_ngroute_neighborhoods!(
            data, 
            ngroute_neighborhood_size; 
            charging_depots_size = ngroute_neighborhood_charging_depots_size,
        )
    end

    some_paths = generate_artificial_paths(data)
    path_costs = compute_path_costs(
        data, 
        some_paths,
    )
    path_service = compute_path_service(
        data, 
        some_paths,
    )
    mp_results = Dict()
    params = Dict()
    params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
    params["objective"] = Float64[]
    params["κ"] = Dict{Int, Float64}[]
    params["μ"] = Dict{Int, Float64}[]
    params["ν"] = Vector{Float64}[]
    params["lp_relaxation_solution_time_taken"] = Float64[]
    params["sp_base_time_taken"] = Float64[]
    params["sp_full_time_taken"] = Float64[]
    params["sp_total_time_taken"] = Float64[]
    params["lp_relaxation_constraint_time_taken"] = Float64[]
    params["number_of_new_paths"] = Int[]

    printlist = String[]
    counter = 0
    converged = false

    add_message!(
        printlist,
        @sprintf(
            """
            Starting column generation on the path formulation.
            # customers:                    %2d
            # depots:                       %2d
            # charging stations:            %2d
            # vehicles:                     %2d
            time windows?:                  %s

            method:                         %s
            subpath_single_service:         %s
            subpath_check_customers:        %s
            path_single_service:            %s
            path_check_customers:           %s
            christofides:                   %s
            ngroute:                        %s
            ngroute neighborhood size:
                customers                   %2d
                charging / depots           %s

            """,
            data["n_customers"],
            data["n_depots"],
            data["n_charging"],
            data["n_vehicles"],
            time_windows,
            method,
            subpath_single_service,
            subpath_check_customers,
            path_single_service,
            path_check_customers,
            christofides,
            ngroute,
            ngroute_neighborhood_size,
            ngroute_neighborhood_charging_depots_size,
        ),
        verbose,
    )
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
        for key in keys(some_paths)
            for p in 1:length(some_paths[key])
    )
    @constraint(
        mp_model,
        κ[i in data["N_depots"]],
        sum(
            sum(
                z[((i,0,data["B"]),state2),p]
                for p in 1:length(some_paths[((i,0,data["B"]),state2)])
            )        
            for (state1, state2) in keys(some_paths)
                if state1[1] == i && state1[2] == 0 && state1[3] == data["B"]
        )
        == data["v_start"][findfirst(x -> (x == i), data["N_depots"])]
    )
    @constraint(
        mp_model,
        μ[n2 in data["N_depots"]],
        sum(
            sum(
                z[(state1, state2),p]
                for p in 1:length(some_paths[(state1, state2)])
            )
            for (state1, state2) in keys(some_paths)
                if state2[1] == n2
        ) ≥ data["v_end"][n2]
    )
    @constraint(
        mp_model,
        ν[j in data["N_customers"]],
        sum(
            sum(
                path_service[((state1, state2),j)][p] * z[(state1, state2),p]
                for p in 1:length(some_paths[(state1, state2)])
            )
            for (state1, state2) in keys(some_paths)
        ) == 1
    )
    @expression(
        mp_model,
        path_costs_expr,
        sum(
            sum(
                path_costs[state_pair][p] * z[state_pair,p]
                for p in 1:length(some_paths[state_pair])
            )
            for state_pair in keys(some_paths)
        )
    )
    @objective(mp_model, Min, path_costs_expr)

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

        if method == "ours"
            (negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
                G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                ;
                ngroute = ngroute,
                subpath_single_service = subpath_single_service,
                subpath_check_customers = subpath_check_customers,
                path_single_service = path_single_service,
                path_check_customers = path_check_customers,
                christofides = christofides,
            )
            (generated_paths) = get_paths_from_negative_path_labels(
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
            (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
                G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                ;
                time_windows = time_windows,
                path_single_service = path_single_service,
                path_check_customers = path_check_customers,
                christofides = christofides,
            )
            generated_paths = get_paths_from_negative_pure_path_labels(
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
        
        if length(generated_paths) == 0
            push!(params["number_of_new_paths"], 0)
            converged = true
        else
            push!(
                params["number_of_new_paths"],
                sum(length(v) for v in values(generated_paths))
            )
        end

        mp_constraint_start_time = time()
        for state_pair in keys(generated_paths)
            if !(state_pair in keys(some_paths))
                some_paths[state_pair] = []
                path_costs[state_pair] = []
                for i in 1:data["n_customers"]
                    path_service[(state_pair, i)] = []
                end
                count = 0
            else
                count = length(some_paths[state_pair])
            end
            for p_new in generated_paths[state_pair]
                if state_pair in keys(some_paths)
                    add = !any(isequal(p_new, s) for s in some_paths[state_pair])
                else
                    add = true
                end
                if add
                    # 1: include in some_paths
                    push!(some_paths[state_pair], p_new)
                    # 2: add path cost
                    push!(
                        path_costs[state_pair], 
                        compute_path_cost(data, p_new)
                    )
                    # 3: add path service
                    for i in 1:data["n_customers"]
                        push!(path_service[(state_pair, i)], p_new.served[i])
                    end
                    # 4: create variable
                    count += 1
                    z[(state_pair, count)] = @variable(mp_model, lower_bound = 0)
                    (state1, state2) = state_pair
                    # 5: modify constraints starting from depot, ending at depot
                    set_normalized_coefficient(κ[state1[1]], z[state_pair,count], 1)
                    set_normalized_coefficient(μ[state2[1]], z[state_pair,count], 1)
                    # 6: modify customer service constraints
                    for l in data["N_customers"]
                        set_normalized_coefficient(ν[l], z[state_pair, count], p_new.served[l])
                    end
                    # 7: modify objective
                    set_objective_coefficient(mp_model, z[state_pair, count], path_costs[state_pair][count])
                end
            end
        end
        mp_constraint_end_time = time()

        push!(
            params["number_of_paths"], 
            sum(length(v) for v in values(some_paths))
        )
        push!(
            params["lp_relaxation_constraint_time_taken"],
            round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %9.3f | %9.3f | %9.3f | %9.3f | %10d \n", 
                counter,
                params["objective"][counter],
                params["number_of_paths"][counter],
                params["lp_relaxation_solution_time_taken"][counter],
                params["sp_base_time_taken"][counter],
                params["sp_full_time_taken"][counter],                
                params["lp_relaxation_constraint_time_taken"][counter],
                params["number_of_new_paths"][counter],
            ),
            verbose,
        )
    end

    for key in keys(some_paths)
        for p in 1:length(some_paths[key])
            JuMP.set_integer(z[key,p])
        end
    end
    cgip_solution_start_time = time()
    @suppress optimize!(mp_model)
    cgip_solution_end_time = time()

    CGLP_results = Dict(
        "objective" => mp_results["objective"],
        "z" => mp_results["z"],
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
        @sprintf("Total time (SP base):         %10.3f s\n", sum(params["sp_base_time_taken"])),
        @sprintf("Total time (SP full):         %10.3f s\n", sum(params["sp_full_time_taken"])),
        @sprintf("Total time (LP construction): %10.3f s\n", sum(params["lp_relaxation_constraint_time_taken"])),
        @sprintf("Total time:                   %10.3f s\n", time_taken),
        @sprintf("\n"),
        @sprintf("(CGIP) Objective:             %.4e\n", CGIP_results["objective"]),
        @sprintf("(CGIP) Total time:            %10.3f s\n", params["CGIP_time_taken"]),
    ]
        add_message!(printlist, message, verbose)
    end
    return CGLP_results, CGIP_results, params, printlist, some_paths
end

function collect_path_solution_support(
    results, 
    paths::Dict{
        Tuple{
            Tuple{Int, Int, Int},
            Tuple{Int, Int, Int},
        }
    },
    ;
)
    results_paths = Tuple{Float64, Path}[]
    for key in keys(paths)
        for p in 1:length(paths[key])
            val = results["z"][key,p]
            if val > 1e-5
                push!(results_paths, (val, paths[key][p]))
            end
        end
    end
    return results_paths
end

function collect_path_solution_metrics!(
    results,
    data, 
    paths,
)
    results["paths"] = collect_path_solution_support(results, paths)

    results["mean_subpath_length"] = sum(
        sum(
            length(s.arcs) for s in p.subpaths
        )
        for (val, p) in results["paths"]
    ) / sum(
        length(p.subpaths)
        for (val, p) in results["paths"]
    )
    results["weighted_mean_subpath_length"] = sum(
        val * sum(
            length(s.arcs) for s in p.subpaths
        )
        for (val, p) in results["paths"]
    ) / sum(
        val * length(p.subpaths)
        for (val, p) in results["paths"]
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