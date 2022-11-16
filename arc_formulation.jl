using JuMP
using Gurobi

function arc_formulation(
    data,
    with_charging::Bool = false,
    with_charging_separate::Bool = false,
    ;
    paths::Union{Dict, Nothing} = nothing,
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
    @variable(model, x[keys(A), N_vehicles], Bin)
    @variable(model, τ_reach[N_nodes, N_vehicles] ≥ 0)
    @variable(model, τ_leave[N_nodes, N_vehicles] ≥ 0)
    @variable(model, l[N_nodes, N_vehicles] ≥ 0)
    if with_charging
        @variable(model, b_start[N_nodes, N_vehicles] ≥ 0)
        @variable(model, b_end[N_nodes, N_vehicles] ≥ 0)
        @variable(model, δ[N_nodes, N_vehicles] ≥ 0)
    end

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
    ); # (1g): save vehicle for a (pickup, dropoff) pair
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
        τ_leave[i,k] + t[i,j] ≤ τ_reach[j,k] + (1 - x[(i,j),k]) * T
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
        l[j,k] ≤ l[i,k] + d[j] + (1 - x[(i,j),k]),
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
            b_start[j,k] ≤ b_end[i,k] - q[i,j] + (1 - x[(i,j),k]) * B
        ); # (1u): charge change based on travel cost
        @constraint( 
            model,
            [i ∈ setdiff(N_nodes, N_depots, N_charging), k ∈ N_vehicles],
            b_end[i,k] ≤ b_start[i,k]
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
    if with_charging
        @expression(
            model,
            total_charge_required, # Total journey time of vehicles
            sum(b_end[i,k] - b_start[i,k] for i in N_depots, k in N_vehicles)
            + μ * sum(δ[i,k] for i in N_charging, k in N_vehicles)
        )
    end

    if with_charging
        @objective(
            model,
            Min,
            total_travel_cost 
            # + total_vehicle_time + 0.9 * total_charge_required
        ); # (1a): objective
    else
        @objective(
            model,
            Min,
            total_travel_cost 
            # + total_vehicle_time 
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
        "time_taken" => time_taken,
        "constraint_time_taken" => constraint_time_taken,
        "solution_time_taken" => solution_time_taken,
    )
    results = Dict(
        "model" => model,
        "objective" => objective_value(model),
        "total_travel_cost" => value(total_travel_cost),
        "total_vehicle_time" => value(total_vehicle_time),
        "x" => value.(x[:,:]),
        "τ_reach" => value.(τ_reach[:,:]),
        "τ_leave" => value.(τ_leave[:,:]),
        "l" => value.(l[:,:]),
    )
    if with_charging
        results["total_charge_required"] = value.(total_charge_required)
        results["b_start"] = value.(b_start)
        results["b_end"] = value.(b_end)
        results["δ"] = value.(δ)
    end
    return results, params
end

