using Suppressor
using Printf

using JuMP
using Gurobi

include("subpath_formulation.jl")
include("utils.jl")

function generate_artificial_paths(data)
    artificial_subpaths = generate_artificial_subpaths(data)
    artificial_paths = Dict(
        key => [
            Path(
                subpaths = [s],
            )
            for s in subpath_list
        ]
        for (key, subpath_list) in pairs(artificial_subpaths)
    )
end

function compute_path_reduced_cost(
    p::Path,
    data, 
    κ,
    μ, 
    ν,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
)
    return sum(
        compute_subpath_reduced_cost(
            s, data, κ, nothing, μ, ν; 
            with_lambda = false,
            with_charging_cost = with_charging_cost,
            time_windows = time_windows,
        )
        for s in p.subpaths
    )
end

function compute_path_cost(
    data,
    p::Path,
    M::Float64 = 1e6,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
)
    return sum(
        compute_subpath_cost(
            data, s,
            ;
            with_charging_cost = with_charging_cost,
            time_windows = time_windows,
        )
        for s in p.subpaths
    )
end

function compute_path_costs(
    data,
    all_paths,
    M::Float64 = 1e6,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
)
    path_costs = Dict(
        key => [
            compute_path_cost(
                data, p, M,
                ;
                with_charging_cost = with_charging_cost,
                time_windows = time_windows,
            )
            for p in all_paths[key]
        ]
        for key in keys(all_paths)
    )
    return path_costs
end

function compute_path_service(
    data,
    p::Path,
    i::Int,
)
    return sum(s.served[i] for s in p.subpaths)
end

function compute_path_services(
    data, 
    all_paths,
)
    path_service = Dict(
        (key, i) => [
            compute_path_service(data, p, i)
            for p in all_paths[key]
        ]
        for key in keys(all_paths), i in 1:data["n_customers"]
    )
    return path_service
end

function path_formulation(
    data,
    all_paths,
    path_costs,
    path_service,
    T_range,
    B_range,
    ;
    integral::Bool = true,
)
    start_time = time()

    model = Model(Gurobi.Optimizer)

    depot_start_states = Set((i,0,data["B"]) for i in data["N_depots"])
    depot_end_states = Set(k[2] for k in keys(all_paths))

    m = maximum(length(x) for x in values(all_paths))

    if integral
        @variable(
            model,
            y[
                k = keys(all_paths),
                p = 1:m; p ≤ length(all_paths[k])
            ],
            Bin
        )
    else
        @variable(
            model,
            y[
                k = keys(all_paths),
                p = 1:m; p ≤ length(all_paths[k])
            ]
            ≥ 0
        )
    end

    @constraint(
        model,
        κ[n1 in data["N_depots"]],
        sum(
            sum(
                y[((n1, 0, data["B"]), end_state), p]
                for p in 1:length(all_paths[((n1, 0, data["B"]), end_state)])
            )   
            for end_state in depot_end_states
                if ((n1, 0, data["B"]), end_state) in keys(all_paths)
        )
        == data["v_start"][findfirst(x -> (x == n1), data["N_depots"])]
    )

    @constraint(
        model,
        μ[n2 in data["N_depots"]],
        sum(
            sum(
                sum(
                    y[(start_state, end_state), p]
                    for p in 1:length(all_paths[(start_state, end_state)])
                )
                for start_state in depot_start_states
                    if (start_state, end_state) in keys(all_paths)
            )
            for end_state in depot_end_states
                if end_state[1] == n2
        ) ≥ data["v_end"][n2]
    )

    @constraint(
        model,
        ν[l in data["N_pickups"]],
        sum(
            sum(
                path_service[(state_pair, l)][p] * y[state_pair, p]
                for p in 1:length(all_paths[state_pair])
            )
            for state_pair in keys(all_paths)
        ) == 1
    )

    @expression(
        model,
        path_costs_expr,
        sum(
            sum(
                path_costs[state_pair][p] * y[state_pair,p]
                for p in 1:length(all_paths[state_pair])
            )
            for state_pair in keys(all_paths)
        )
    )
    @objective(model, Min, path_costs_expr)

    constraint_end_time = time() 

    optimize!(model)

    end_time = time()

    time_taken = end_time - start_time
    constraint_time_taken = constraint_end_time - start_time
    solution_time_taken = end_time - start_time

    params = Dict(
        "time_taken" => round(time_taken, digits=3),
        "constraint_time_taken" => round(constraint_time_taken, digits=3),
        "solution_time_taken" => round(solution_time_taken, digits=3),
    )
    results = Dict(
        "objective" => objective_value(model),
        "y" => value.(y),
    )

    if !integral
        results["κ"] = Dict(zip(data["N_depots"], dual.(model[:κ]).data))
        results["μ"] = Dict(zip(data["N_depots"], dual.(model[:μ]).data))
        results["ν"] = dual.(model[:ν]).data
    end
    return results, params
end

function generate_paths(
    G, 
    data, 
    T_range,
    B_range,
    κ, 
    μ, 
    ν,
    ;
    time_windows::Bool = true,
    with_charging_cost::Bool = false,
)
    generated_paths = Dict{
        Tuple{Tuple{Int, Float64, Float64}, Tuple{Int, Float64, Float64}}, 
        Vector{Path},
    }()
    sp_max_time_taken = 0.0
    for starting_node in data["N_depots"]
        r = @timed find_smallest_reduced_cost_paths(
            starting_node, 
            G, data, T_range, B_range, 
            κ, μ, ν, 
            ;
            time_windows = time_windows,
            with_charging_cost = with_charging_cost,
        )
        labels = r.value
        if sp_max_time_taken < r.time
            sp_max_time_taken == r.time
        end
        state1 = (starting_node, 0.0, data["B"])
        for (end_node, path_dict) in pairs(labels)
            for ((time, charge), path) in pairs(path_dict)
                if path.cost ≥ -1e-6
                    continue
                end
                round_time = dceil(time, T_range)
                round_charge = dfloor(charge, B_range)
                state2 = (end_node, round_time, round_charge)
                state_pair = (state1, state2)
                if state_pair in keys(generated_paths)
                    push!(generated_paths[state_pair], Path(path))
                else
                    generated_paths[state_pair] = [Path(path)]
                end
            end
        end
    end
    return generated_paths, sp_max_time_taken
end

function path_formulation_column_generation(
    G,
    data,
    T_range,
    B_range,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
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

    artificial_paths = generate_artificial_paths(data)
    if with_heuristic
        some_paths = construct_heuristic_paths(data, T_range, B_range)
        for (key, artificial_path_list) in pairs(artificial_paths)
            if !(key in keys(some_paths))
                some_paths[key] = []
            end
            append!(some_paths[key], artificial_path_list)
        end
    else
        some_paths = deepcopy(artificial_paths)
    end
    mp_results = Dict()
    params = Dict()
    params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
    params["number_of_keys"] = [length(some_paths)]
    params["objective"] = Float64[]
    params["κ"] = Dict{Int, Float64}[]
    params["μ"] = Dict{Int, Float64}[]
    params["ν"] = Vector{Float64}[]
    params["lp_relaxation_time_taken"] = Float64[]
    params["lp_relaxation_constraint_time_taken"] = Float64[]
    params["lp_relaxation_solution_time_taken"] = Float64[]
    params["sp_total_time_taken"] = Float64[]
    params["sp_max_time_taken"] = Float64[]
    params["number_of_current_paths"] = Int[]
    
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
        @sprintf("              |  Objective | # paths | Time (LP) | Time (SP) | # new paths \n"),
        verbose,
    )
    while (
        !converged
        && time_limit > time() - start_time
    )
        counter += 1
        path_costs = compute_path_costs(
            data, 
            some_paths,
            ;
            with_charging_cost = with_charging_cost,
            time_windows = time_windows,
        )
        path_service = compute_path_services(
            data, 
            some_paths,
        )
        mp_results, mp_params = @suppress path_formulation(
            data, 
            some_paths,
            path_costs, 
            path_service,
            T_range,
            B_range,
            ;
            integral = false,
        )

        push!(params["objective"], mp_results["objective"])
        push!(params["κ"], mp_results["κ"])
        push!(params["μ"], mp_results["μ"])
        push!(params["ν"], mp_results["ν"])
        push!(params["lp_relaxation_time_taken"], mp_params["time_taken"])
        push!(params["lp_relaxation_constraint_time_taken"], mp_params["constraint_time_taken"])
        push!(params["lp_relaxation_solution_time_taken"], mp_params["solution_time_taken"])

        # generate paths
        generate_paths_result = @timed generate_paths(
            G, data, T_range, B_range, 
            mp_results["κ"], mp_results["μ"], mp_results["ν"]
            ;
            with_charging_cost = with_charging_cost,
            time_windows = time_windows,
        )
        (current_paths, sp_max_time_taken) = generate_paths_result.value

        push!(
            params["sp_total_time_taken"],
            round(generate_paths_result.time, digits=3)
        )
        push!(
            params["sp_max_time_taken"],
            round(sp_max_time_taken, digits=3)
        )
        if length(current_paths) == 0
            push!(params["number_of_current_paths"], 0)
            converged = true
        else
            push!(
                params["number_of_current_paths"],
                sum(length(v) for v in values(current_paths))
            )
        end
        for key in keys(current_paths)
            if !(key in keys(some_paths))
                some_paths[key] = current_paths[key]
            else
                for p_new in current_paths[key]
                    if !any(isequal(p_new, p) for p in some_paths[key])
                        push!(some_paths[key], p_new)
                    end
                end
            end
        end
        push!(
            params["number_of_keys"],
            length(some_paths)
        )
        push!(
            params["number_of_paths"], 
            sum(length(v) for v in values(some_paths))
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %7d | %9.3f | %9.3f | %11d \n", 
                counter,
                params["objective"][counter],
                params["number_of_paths"][counter],
                params["lp_relaxation_time_taken"][counter],
                params["sp_total_time_taken"][counter],
                params["number_of_current_paths"][counter],
            ),
            verbose,
        )
    end
    results = Dict(
        "objective" => mp_results["objective"],
        "y" => mp_results["y"],
    )
    params["counter"] = counter
    end_time = time()
    time_taken = round(end_time - start_time, digits = 3)
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
    return results, params, printlist, some_paths
end