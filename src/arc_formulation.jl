using JuMP
using Gurobi
using Printf

function arc_formulation(
    data,
    with_charging::Bool = false,
    with_charging_separate::Bool = false,
    ;
    with_charging_cost::Bool = false,
    with_customer_delay_cost::Bool = false,
    integral::Bool = true,
    paths::Union{Dict, Nothing} = nothing,
    exclude_arcs::Vector = [],
    time_limit::Union{Float64, Int} = 60.0,
    formulate_only::Bool = false,
)

    if with_charging_separate
        if !with_charging
            error("""
            Cannot impose charging solution without enabling charging 
            (`with_charging == false`)
            """)
        elseif isnothing(paths)
            error("""
            For `with_charging_separate == true`, must provide paths 
            of solution obtained without charging!
            """)
        end
    end

    n_nodes = data["n_nodes"]
    n_customers = data["n_customers"]
    n_depots = data["n_depots"]
    n_vehicles = data["n_vehicles"]

    N_pickups = data["N_pickups"]
    N_dropoffs = data["N_dropoffs"]
    N_depots = data["N_depots"]
    N_vehicles = data["N_vehicles"]
    N_nodes = data["N_nodes"]
    if with_charging
        n_charging = data["n_charging"]
        N_charging = data["N_charging"]
    end

    V = data["V"]
    v_start = data["v_start"]
    v_end = data["v_end"]

    c = data["c"]
    t = data["t"]
    C = data["C"]
    d = data["d"]

    T = data["T"]
    A = data["A"]

    α = data["α"]
    β = data["β"]
    
    if with_charging
        q = data["q"]
        B = data["B"]
        μ = data["μ"]
    end

    start_time = time()

    model = Model(Gurobi.Optimizer)
    set_time_limit_sec(model, time_limit)
    if integral
        @variable(model, x[keys(A), N_vehicles], Bin)
    else
        @variable(model, 0 ≤ x[keys(A), N_vehicles] ≤ 1)
    end
    @variable(model, τ_reach[N_nodes, N_vehicles] ≥ 0)
    @variable(model, τ_leave[N_nodes, N_vehicles] ≥ 0)
    @variable(model, l[N_nodes, N_vehicles] ≥ 0)
    if with_charging
        @variable(model, b_start[N_nodes, N_vehicles] ≥ 0)
        @variable(model, b_end[N_nodes, N_vehicles] ≥ 0)
        @variable(model, δ[N_charging, N_vehicles] ≥ 0)
    end

    @constraint(
        model,
        [k in N_vehicles, (i,j) in exclude_arcs],
        x[(i,j),k] == 0,
    )

    if with_charging_separate
        for k in N_vehicles
            for a in paths[k]
                if (a[1] in N_depots && a[2] in N_depots)
                    @constraint(
                        model, 
                        [(i,j) in setdiff(keys(A), [a])], 
                        x[(i,j),k] == 0
                    )
                    @constraint(
                        model, 
                        x[a,k] == 1
                    )
                else
                    # This includes (for now) (N_pickups, N_pickups) 
                    # and (N_dropoffs, N_dropoffs) pairs,
                    # which are not possible when C = 1.

                    # cannot go anywhere else from a[1], aside from a[2] and CS
                    for j in setdiff(N_nodes, [a[2]], N_charging)
                        if (a[1],j) in keys(A)
                            @constraint(model, x[(a[1],j),k] == 0)
                        end
                    end
                    # cannot arrive at a[2] from anywhere else, aside from a[1] and CS
                    for i in setdiff(N_nodes, [a[1]], N_charging)
                        if (i,a[2]) in keys(A)
                            @constraint(model, x[(i,a[2]),k] == 0)
                        end
                    end
                end
            end
        end    
    end

    @constraint(
        model, 
        [i ∈ N_depots, k ∈ V[i]],
        sum(x[(i,j),k] for j in N_nodes if (i,j) in keys(A)) == 1
    ); # (1c): leave from depot
    @constraint(
        model, 
        [i ∈ N_depots, j ∈ N_nodes, k ∈ setdiff(N_vehicles, V[i]); (i,j) in keys(A)],
        x[(i,j),k] == 0
    ); # (1d): leave from depot (others)
    @constraint(
        model,
        [i ∈ setdiff(N_nodes, N_depots), k ∈ N_vehicles],
        sum(x[(i,j),k] for j in N_nodes if (i,j) in keys(A))
        == sum(x[(j,i),k] for j in N_nodes if (j,i) in keys(A))
    ); # (1e): flow conservation (outside of depots)
    @constraint(
        model,
        [i ∈ N_pickups],
        sum(x[(i,j),k] for k in N_vehicles, j in N_nodes if (i,j) in keys(A)) == 1
    ); # (1f): all pickups served
    @constraint(
        model,
        [i ∈ N_pickups, k ∈ N_vehicles],
        sum(x[(i,j),k] for j in N_nodes if (i,j) in keys(A))
        == sum(x[(j,i+n_customers),k] for j in N_nodes if (j, i+n_customers) in keys(A))
    ); # (1g): same vehicle for a (pickup, dropoff) pair
    @constraint(
        model, 
        [k ∈ N_vehicles],
        sum(x[(j,i),k] for i in N_depots, j in N_nodes if (j,i) in keys(A)) == 1
    ); # (1h): all vehicles reach a depot
    @constraint(
        model, 
        [i ∈ N_depots],
        sum(x[(j,i),k] for k in N_vehicles, j in N_nodes if (j, i) in keys(A)) ≥ v_end[i]
    ); # (1i): each depot receives at least its required quota
    @constraint(
        model, 
        [i ∈ N_depots, k ∈ N_vehicles],
        τ_leave[i, k] == 0
    ); # (1j): time leaving the depot
    @constraint(
        model, 
        [i ∈ N_pickups, k ∈ N_vehicles],
        τ_leave[i, k] + t[i, i+n_customers] ≤ τ_reach[i+n_customers, k]
    ); # (1k): pickups reached before its dropoff
    if with_charging
        @constraint(
            model,
            [i ∈ setdiff(N_nodes, N_charging, N_depots), k ∈ N_vehicles],
            τ_leave[i,k] ≥ τ_reach[i,k],
        ); # logical: linking time reaching a node and time leaving it (for non-depots / non-charging stations)
        @constraint(
            model,
            [i ∈ N_charging, k ∈ N_vehicles],
            τ_leave[i,k] ≥ τ_reach[i,k] + δ[i,k]
        ); # logical: linking time reaching a node and time leaving it (for charging stations)    
    else
        @constraint(
            model,
            [i ∈ setdiff(N_nodes, N_depots), k ∈ N_vehicles],
            τ_leave[i,k] ≥ τ_reach[i,k],
        ); # logical: linking time reaching a node and time leaving it (for non-depots)
    end

    @constraint(
        model,
        [(i,j) in keys(A), k ∈ N_vehicles],
        τ_leave[i,k] + t[i,j] ≤ τ_reach[j,k] + (1 - x[(i,j),k]) * 2 * T
    ); # (1l, 1m): precedence conditions
    
    @constraint(
        model, 
        [i ∈ N_nodes, k ∈ N_vehicles],
        α[i] ≤ τ_leave[i,k] ≤ β[i]
    ); # (1n): time window constraints
    @constraint(
        model, 
        [i ∈ N_nodes, k ∈ N_vehicles],
        α[i] ≤ τ_reach[i,k] ≤ β[i]
    ); # (1n): time window constraints
    
    @constraint(
        model,
        [i ∈ N_depots, k ∈ N_vehicles],
        l[i,k] == 0
    ); # (1o): load leaving depot
    @constraint(
        model,
        [(i,j) in keys(A), k ∈ N_vehicles],
        l[j,k] ≤ l[i,k] + d[j] + (1 - x[(i,j),k]) * 2,
    ); # (1p, 1q, 1r): load after serving a node
    @constraint(
        model,
        [i ∈ N_nodes, k ∈ N_vehicles],
        0 ≤ l[i,k] ≤ C,
    ); # (1s): load within capacity
    
    if with_charging
        @constraint(
            model,
            [i ∈ N_depots, k ∈ N_vehicles],
            b_end[i,k] == B
        ); # (1t): charge leaving depot
        @constraint(
            model,
            [(i,j) in keys(A), k ∈ N_vehicles],
            b_start[j,k] ≤ b_end[i,k] - q[i,j] + (1 - x[(i,j),k]) * 2 * B
        ); # (1u): charge change based on travel cost
        @constraint( 
            model,
            [i ∈ setdiff(N_nodes, N_depots, N_charging), k ∈ N_vehicles],
            b_end[i,k] == b_start[i,k]
        ); # (1v): charge leaving not more than charge reaching a node (non-charging stations)
        @constraint(
            model,
            [i ∈ N_charging, k ∈ N_vehicles],
            b_end[i,k] ≤ b_start[i,k] + μ * δ[i,k]
        ); # (1w): charge leaving not more than charge reaching a node + amount charged (charging stations)
        @constraint(
            model,
            [i ∈ N_nodes, k ∈ N_vehicles],
            0 ≤ b_start[i,k] ≤ B
        ); # (1x): charge within capacity
        @constraint(
            model,
            [i ∈ N_nodes, k ∈ N_vehicles],
            0 ≤ b_end[i,k] ≤ B
        ); # (1x): charge within capacity
        @constraint(
            model,
            δ0[i ∈ N_charging, k ∈ N_vehicles],
            δ[i,k] ≤ sum(x[(j,i),k] for j in N_nodes if (j,i) in keys(A)) * 2 * (B / μ)
        )
    end
    
    @expression(
        model,
        total_travel_cost,
        sum(c[i,j] * x[(i,j),k] for (i,j) in keys(A), k ∈ N_vehicles)
    )
    @expression(
        model,
        total_vehicle_time, # Total journey time of vehicles
        sum(τ_reach[i,k] for i in N_depots, k in N_vehicles)
    )
    @expression(
        model,
        total_customer_delay_cost,
        sum(τ_reach[i,k] - α[i] for i in N_dropoffs, k in N_vehicles)
    )
    if with_charging
        @expression(
            model,
            total_charge_cost, 
            sum(δ[i,k] for i in N_charging, k in N_vehicles)
        )
    end

    if (with_charging && with_charging_cost) && with_customer_delay_cost
        @objective(
            model,
            Min,
            data["travel_cost_coeff"] * total_travel_cost 
            + data["charge_cost_coeff"] * total_charge_cost
            + data["customer_delay_cost_coeff"] * total_customer_delay_cost
        ); # (1a): objective
    elseif (with_charging && with_charging_cost)
        @objective(
            model,
            Min,
            data["travel_cost_coeff"] * total_travel_cost 
            + data["charge_cost_coeff"] * total_charge_cost
        ); # (1a): objective
    elseif with_customer_delay_cost
        @objective(
            model,
            Min,
            data["travel_cost_coeff"] * total_travel_cost 
            + data["customer_delay_cost_coeff"] * total_customer_delay_cost
        ); # (1a): objective
    else
        @objective(
            model,
            Min,
            data["travel_cost_coeff"] * total_travel_cost
        ); # (1a): objective
    end

    if formulate_only
        return
    end

    constraint_end_time = time()

    optimize!(model)
    end_time = time()
    time_taken = end_time - start_time
    constraint_time_taken = constraint_end_time - start_time
    solution_time_taken = end_time - constraint_end_time

    params = Dict(
        "time_taken" => round(time_taken, digits = 3),
        "constraint_time_taken" => round(constraint_time_taken, digits = 3),
        "solution_time_taken" => round(solution_time_taken, digits = 3),
    )
    results = Dict{Any, Any}(
        "model" => model,
    )
    if (
        termination_status(model) == OPTIMAL
        || (
            termination_status(model) == TIME_LIMIT 
            && has_values(model)
        )
    )
        if termination_status(model) == OPTIMAL
            results["status"] = "optimal"
        else
            results["status"] = "time_limit_with_values"
        end
        results["objective"] = objective_value(model)
        results["total_travel_cost"] = value(total_travel_cost)
        results["total_vehicle_time"] = value(total_vehicle_time)
        results["total_customer_delay_cost"] = value(total_customer_delay_cost)
        results["x"] = value.(x[:,:])
        results["τ_reach"] = value.(τ_reach[:,:])
        results["τ_leave"] = value.(τ_leave[:,:])
        results["l"] = value.(l[:,:])
        if with_charging
            results["total_charge_cost"] = value.(total_charge_cost)
            results["b_start"] = value.(b_start)
            results["b_end"] = value.(b_end)
            results["δ"] = value.(δ)
        end
    elseif (
        termination_status(model) == TIME_LIMIT
        && !has_values(model)
    )
        results["status"] = "time_limit_no_values"
    elseif termination_status(model) == INFEASIBLE
        results["status"] = "infeasible"
    else
        results["status"] = "unknown"
    end
    return results, params
end

function construct_paths_from_arc_solution(
    results,
    data,
)
    paths = Dict()
    for k in data["N_vehicles"] 
        arcs = [
            (i,j) for (i,j) in keys(data["A"]) 
            if results["x"][(i,j),k] > 0.5
        ]
        path = []
        i = findfirst(x -> (x[1] ∈ data["N_depots"]), arcs)
        while true
            a = popat!(arcs, i)
            push!(path, a)
            current_node = a[2]
            if length(arcs) == 0
                break
            end
            i = findfirst(x -> (x[1] == current_node), arcs)
        end
        paths[k] = path
    end
    return paths
end

function arc_results_printout(
    results, 
    params,
    data,
    with_charging::Bool = false,
)

    printlist = [
        @sprintf("Objective:                 %7.1f\n", results["objective"]),
        @sprintf("Total travel cost:         %7.1f\n", results["total_travel_cost"]),
        @sprintf("Total vehicle time:        %7.1f\n", results["total_vehicle_time"]),
        @sprintf("Total customer delay cost: %7.1f\n", results["total_customer_delay_cost"]),
    ]
    if with_charging
        push!(
            printlist, 
            @sprintf("Total charge cost:     %7.1f\n", results["total_charge_cost"])
        )
    end
    append!(
        printlist, 
        [
            @sprintf("Time taken:                %7.1f s\n", params["time_taken"]),
            @sprintf("Time taken (formulation):  %7.1f s\n", params["constraint_time_taken"]),
            @sprintf("Time taken (solving):      %7.1f s\n", params["solution_time_taken"]),
            "\n",
            @sprintf("Time taken:                %7.1f s\n", params["time_taken"]),
            @sprintf("Time taken (formulation):  %7.1f s\n", params["constraint_time_taken"]),
            @sprintf("Time taken (solving):      %7.1f s\n", params["solution_time_taken"]),
            "\n",
        ]
    )

    if with_charging
        append!(
            printlist, 
            [
                "Vehicles         From      time   charge             To      time   charge  |  arc_time |  arc_charge \n",
                "----------------------------------------------------------------------------|-----------|-------------\n",
            ]
        )
        else
        append!(
            printlist, 
            [
                "Vehicles         From      time             To      time  |  arc_time \n",
                "----------------------------------------------------------|-----------\n",
            ]
        )
    end

    paths = construct_paths_from_arc_solution(results, data)

    for k in data["N_vehicles"]
        for a in paths[k]
            if with_charging
                if a[2] in data["N_charging"]
                    push!(
                        printlist, 
                        @sprintf(
                            "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | -  %6.1f | -    %6.1f \n", 
                            k, 
                            data["node_labels"][a[1]], 
                            results["τ_leave"][a[1],k],
                            results["b_end"][a[1],k],
                            data["node_labels"][a[2]],
                            results["τ_reach"][a[2],k],
                            results["b_start"][a[2],k],
                            data["t"][a[1],a[2]],
                            data["q"][a[1],a[2]],
                        )
                    )
                    push!(
                        printlist, 
                        @sprintf(
                            "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | -  %6.1f | +    %6.1f \n", 
                            k, 
                            data["node_labels"][a[2]], 
                            results["τ_reach"][a[2],k],
                            results["b_start"][a[2],k],
                            data["node_labels"][a[2]], 
                            results["τ_leave"][a[2],k],
                            results["b_end"][a[2],k],
                            results["δ"][a[2],k],
                            results["δ"][a[2],k] * data["μ"],
                        )
                    )
                else
                    push!(
                        printlist, 
                        @sprintf(
                            "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | -  %6.1f | -    %6.1f \n", 
                            k, 
                            data["node_labels"][a[1]], 
                            results["τ_leave"][a[1],k],
                            results["b_end"][a[1],k],
                            data["node_labels"][a[2]],
                            results["τ_reach"][a[2],k],
                            results["b_start"][a[2],k],
                            data["t"][a[1],a[2]],
                            data["q"][a[1],a[2]],
                        )
                    )
                end
            else
                push!(
                    printlist,
                    @printf(
                        "Vehicle %s: %12s (%6.1f) -> %12s (%6.1f) | -  %6.1f\n", 
                        k, 
                        data["node_labels"][a[1]], 
                        results["τ_leave"][a[1],k],
                        data["node_labels"][a[2]],
                        results["τ_reach"][a[2],k],
                        data["t"][a[1],a[2]],
                    )
                )
            end
        end
        push!(printlist, "\n")
    end

    for message in printlist
        print(message)
    end
    return printlist
end
