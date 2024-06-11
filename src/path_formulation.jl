using JuMP
using Gurobi
using Suppressor

include("utils.jl")
include("desaulniers_benchmark.jl")
include("subpath_stitching.jl")

function generate_artificial_paths(
    data::EVRPData,
    graph::EVRPGraph,
)
    artificial_paths = Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    }()
    start_depots = zeros(Int, data.n_vehicles)
    for (k, v_list) in pairs(data.V)
        for v in v_list
            start_depots[v] = k
        end
    end
    end_depots = []
    for k in graph.N_depots
        append!(end_depots, repeat([k], data.v_end[k]))
    end
    append!(end_depots,
        repeat(
            graph.N_depots, 
            outer = Int(ceil((data.n_vehicles - sum(values(data.v_end))) / graph.n_depots))
        )
    )
    end_depots = end_depots[1:data.n_vehicles]
    
    for (v, (starting_node, current_node)) in enumerate(zip(start_depots, end_depots))
        starting_time = 0.0
        starting_charge = graph.B
        key = (
            (starting_node, starting_time, starting_charge),  
            (current_node, starting_time, starting_charge)
        )
        # initialise a proportion of the customers to be served
        served = zeros(Int, graph.n_customers)
        for i in 1:length(served)
            if mod1(i, data.n_vehicles) == v
                served[i] = 1
            end
        end
        s = Subpath(
            n_customers = graph.n_customers,
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
            arcs = [(starting_node, current_node)],
            customer_arcs = NTuple{2, Int}[],
            artificial = true,
        )
        if !(key in keys(artificial_paths))
            artificial_paths[key] = []
        end
        push!(artificial_paths[key], p)
    end
    return artificial_paths
end

function add_path_to_generated_paths!(
    generated_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Path},
    },
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

function convert_path_label_to_path(
    path_label::PathLabel,
    data::EVRPData,
    graph::EVRPGraph,
)
    current_time, current_charge = (0.0, graph.B)
    prev_time, prev_charge = current_time, current_charge
    s_labels = copy(path_label.subpath_labels)
    deltas = copy(path_label.charging_actions)
    p = Path(
        subpaths = Subpath[],
        charging_arcs = ChargingArc[],
        served = zeros(Int, graph.n_customers),
        arcs = NTuple{2, Int}[],
        customer_arcs = NTuple{2, Int}[],
    )
    while true
        s_label = popfirst!(s_labels)
        prev_time = current_time
        prev_charge = current_charge
        current_node = s_label.nodes[end]
        current_time = current_time + s_label.time_taken
        current_charge = current_charge - s_label.charge_taken
        s = Subpath(
            n_customers = graph.n_customers,
            starting_node = s_label.nodes[1],
            starting_time = prev_time,
            starting_charge = prev_charge,
            current_node = current_node,
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
            starting_node = current_node, 
            starting_time = prev_time, 
            starting_charge = prev_charge, 
            delta = delta,
            charge_cost_coeff = data.charge_cost_coeffs[current_node],
            current_time = current_time, 
            current_charge = current_charge,
        )
        push!(p.charging_arcs, a)
    end
    p.served = sum(s.served for s in p.subpaths)
    p.arcs = vcat([s.arcs for s in p.subpaths]...)
    customers = [a[1] for a in p.arcs if a[1] in graph.N_customers]
    p.customer_arcs = collect(zip(customers[1:end-1], customers[2:end]))
    return p
end

function convert_pure_path_label_to_path(
    pure_path_label::PurePathLabel,
    data::EVRPData,
    graph::EVRPGraph,
)
    p = Path(
        subpaths = Subpath[],
        charging_arcs = ChargingArc[],
        served = zeros(Int, graph.n_customers),
        arcs = NTuple{2, Int}[],
        customer_arcs = NTuple{2, Int}[],
    )
    states = NTuple{3, Int}[]
    current_subpath = Subpath(
        n_customers = graph.n_customers,
        starting_node = pure_path_label.nodes[1],
        starting_time = 0.0, 
        starting_charge = graph.B,
    )
    i = pure_path_label.nodes[1]
    for (j, e, s) in zip(pure_path_label.nodes[2:end], pure_path_label.excesses, pure_path_label.slacks)
        current_subpath.current_node = j
        push!(current_subpath.arcs, (i, j))
        current_subpath.starting_time += (e + s)
        current_subpath.starting_charge += (e + s) 
        current_subpath.current_time += (graph.t[i,j] + e + s)
        current_subpath.current_charge += (- graph.q[i,j] + e + s)
        if j in graph.N_charging
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
                n_customers = graph.n_customers,
                starting_node = j,
                starting_time = current_subpath.current_time, 
                starting_charge = current_subpath.current_charge,
            )
        elseif j in graph.N_customers
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
                data.charge_cost_coeffs[states[2*i][1]],
                states[2*i+1][2],
                states[2*i+1][3],
            )
        )
    end
    p.served = sum(s.served for s in p.subpaths)
    p.arcs = vcat([s.arcs for s in p.subpaths]...)
    customers = [a[1] for a in p.arcs if a[1] in graph.N_customers]
    p.customer_arcs = collect(zip(customers[1:end-1], customers[2:end]))
    return p
end

function get_paths_from_negative_path_labels(
    data::EVRPData,
    graph::EVRPGraph,
    path_labels::Vector{PathLabel},
)
    generated_paths = Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Path},
    }()
    for path_label in path_labels
        p = convert_path_label_to_path(path_label, data, graph)
        add_path_to_generated_paths!(generated_paths, p)
    end
    return generated_paths
end

function get_paths_from_negative_pure_path_labels(
    data::EVRPData,
    graph::EVRPGraph,
    pure_path_labels::Vector{PurePathLabel},
)
    generated_paths = Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Path},
    }()
    for pure_path_label in pure_path_labels
        p = convert_pure_path_label_to_path(pure_path_label, data, graph)
        add_path_to_generated_paths!(generated_paths, p)
    end
    return generated_paths
end

function path_formulation_build_model(
    data::EVRPData,
    graph::EVRPGraph,
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Path},
    },
    path_costs::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Int},
    },
    path_service::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
            Int,
        }, 
        Vector{Int}
    },
    ;
    Env = nothing,
)
    if isnothing(Env)
        model = @suppress Model(Gurobi.Optimizer)
    else
        model = @suppress Model(() -> Gurobi.Optimizer(Env))
    end
    JuMP.set_attribute(model, "MIPGapAbs", 1e-3)
    JuMP.set_string_names_on_creation(model, false)

    z = Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    }(
        (key, p) => @variable(model, lower_bound = 0)
        for key in keys(some_paths)
            for p in 1:length(some_paths[key])
    )
    @constraint(
        model, 
        κ[i in graph.N_depots],
        sum(
            sum(
                z[((i,0,graph.B),state2),p]
                for p in 1:length(some_paths[((i,0,graph.B),state2)])
            )
            for (state1, state2) in keys(some_paths)
                if state1[1] == i && state1[2] == 0 && state1[3] == graph.B
        )
        == data.v_start[findfirst(x -> (x == i), graph.N_depots)]
    )
    @constraint(
        model,
        μ[n2 in graph.N_depots],
        sum(
            sum(
                z[(state1, state2),p]
                for p in 1:length(some_paths[(state1, state2)])
            )
            for (state1, state2) in keys(some_paths)
                if state2[1] == n2
        )
        ≥ data.v_end[n2]
    )
    @constraint(
        model,
        ν[j in graph.N_customers],
        sum(
            sum(
                path_service[((state1, state2), j)][p] * z[(state1, state2), p]
                for p in 1:length(some_paths[(state1, state2)])
            )
            for (state1, state2) in keys(some_paths)
        )
        == 1
    )
    @expression(
        model, 
        path_costs_expr,
        sum(
            sum(
                path_costs[state_pair][p] * z[state_pair,p]
                for p in 1:length(some_paths[state_pair])
            )
            for state_pair in keys(some_paths)
        )
    )
    @objective(model, Min, path_costs_expr)

    return model, z
end

function check_path_in_SR3_constraint(
    path::Path,
    S::NTuple{3, Int},
)
    return Int(floor(sum(path.served[collect(S)]) / 2))
end

function compute_path_coefficient_in_lmSRnk_constraint(
    path::Path,
    S::Tuple{Vararg{Int}},
    M::Tuple{Vararg{Int}},
    k::Int,
)
    nodes = get_nodes(path)
    coeff = 0
    state = 0
    for node in nodes
        if !(node in M)
            state = 0
        elseif node in S
            state += 1
            if state ≥ k
                coeff += 1
                state -= k
            end
        end
    end
    return coeff
end

function add_paths_to_path_model!(
    model::Model,
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
            Int,
        }, 
        VariableRef,
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    path_costs::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Int},
    },
    path_service::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        Vector{Int},
    },
    generated_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    SR3_constraints::Dict{
        T, 
        ConstraintRef, 
    },
    data::EVRPData,
    graph::EVRPGraph,
) where {T}
    mp_constraint_start_time = time()
    for state_pair in keys(generated_paths)
        if !(state_pair in keys(some_paths))
            some_paths[state_pair] = Path[]
            path_costs[state_pair] = Int[]
            for i in 1:graph.n_customers
                path_service[(state_pair, i)] = Int[]
            end
            count = 0
        else
            count = length(some_paths[state_pair])
        end
        for p_new in generated_paths[state_pair]
            if state_pair in keys(some_paths)
                add = !any(isequal(p_new, p) for p in some_paths[state_pair])
            else
                add = true
            end
            if add
                count += 1
                # 1: include in some_paths
                push!(some_paths[state_pair], p_new)
                # 2: add path cost
                push!(
                    path_costs[state_pair], 
                    compute_path_cost(data, graph, p_new)
                )
                # 3: add path service
                for i in 1:graph.n_customers
                    push!(path_service[(state_pair, i)], p_new.served[i])
                end
                # 4: create variable
                z[(state_pair, count)] = @variable(model, lower_bound = 0)
                (state1, state2) = state_pair
                # 5: modify constraints starting from depot, ending at depot
                set_normalized_coefficient(model[:κ][state1[1]], z[state_pair,count], 1)
                set_normalized_coefficient(model[:μ][state2[1]], z[state_pair,count], 1)
                # 6: modify customer service constraints
                for l in graph.N_customers
                    set_normalized_coefficient(model[:ν][l], z[state_pair, count], p_new.served[l])
                end
                # 7: modify objective
                set_objective_coefficient(model, z[state_pair, count], path_costs[state_pair][count])
                # 8: add variable to violated SR3 constraints (if applicable)
                if keytype(SR3_constraints) == NTuple{3, Int}
                    for (SR3_key, SR3_con) in pairs(SR3_constraints)
                        val = check_path_in_SR3_constraint(p_new, SR3_key)
                        if val > 0
                            set_normalized_coefficient(SR3_con, z[state_pair,count], val)
                        end
                    end
                elseif keytype(SR3_constraints) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
                    for ((S, M), SR3_con) in pairs(SR3_constraints)
                        val = compute_path_coefficient_in_lmSRnk_constraint(p_new, S, M, 2)
                        if val > 0
                            set_normalized_coefficient(SR3_con, z[state_pair,count], val)
                        end
                    end
                end
            end
        end
    end
    mp_constraint_end_time = time()
    mp_constraint_time = round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
    return mp_constraint_time
end

function path_formulation_column_generation!(
    model::Model,
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    },
    SR3_constraints::Dict{
        T,
        ConstraintRef,
    },
    data::EVRPData,
    graph::EVRPGraph,
    artificial_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    path_costs::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Int},
    },
    path_service::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        Vector{Int},
    },
    printlist::Vector{String},
    ;
    method::String = "ours",
    charge_cost_heterogenous::Bool = false,
    time_windows::Bool = false,
    elementary::Bool = true,
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    verbose::Bool = true,
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
) where {T}

    start_time = time()
    counter = 0
    converged = false
    local CGLP_results = Dict{String, Any}()

    CG_params = Dict{String, Any}()
    CG_params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
    CG_params["objective"] = Float64[]
    CG_params["κ"] = Dict{Int, Float64}[]
    CG_params["μ"] = Dict{Int, Float64}[]
    CG_params["ν"] = Vector{Float64}[]
    CG_params["λ"] = Dict{keytype(SR3_constraints), Float64}[]
    CG_params["lp_relaxation_solution_time_taken"] = Float64[]
    CG_params["sp_base_time_taken"] = Float64[]
    CG_params["sp_full_time_taken"] = Float64[]
    CG_params["sp_total_time_taken"] = Float64[]
    CG_params["lp_relaxation_constraint_time_taken"] = Float64[]
    CG_params["number_of_new_paths"] = Int[]
    CG_params["converged"] = false
    CG_params["cycled"] = false

    while (
        !converged
        && time_limit ≥ (time() - start_time)
        && max_iters > counter
    )
        counter += 1
        mp_solution_start_time = time()
        @suppress optimize!(model)
        mp_solution_end_time = time()
        if JuMP.termination_status(model) in [
            JuMP.INFEASIBLE,
            JuMP.INFEASIBLE_OR_UNBOUNDED,
            JuMP.TIME_LIMIT,
            JuMP.MEMORY_LIMIT,
        ]
            CGLP_results = Dict(
                "errored" => true,
                "artificial" => false,
            )
            add_message!(printlist, "CGLP model errored.\n", verbose)
            break
        end
        CGLP_results = Dict(
            "errored" => false,
            "objective" => objective_value(model),
            "z" => Dict(
                (key, p) => value.(z[(key, p)])
                for (key, p) in keys(z)
            ),
            "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
            "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
            "ν" => dual.(model[:ν]).data,
            "λ" => Dict{keytype(SR3_constraints), Float64}(
                S => dual(SR3_constraints[S])
                for S in keys(SR3_constraints)
            ),
        )
        CGLP_results["artificial"] = any(
            value.(z[(key, p)]) > 1e-3
            for key in keys(artificial_paths)
                for p in 1:length(artificial_paths[key])
        )
        push!(CG_params["objective"], CGLP_results["objective"])
        push!(CG_params["κ"], CGLP_results["κ"])
        push!(CG_params["μ"], CGLP_results["μ"])
        push!(CG_params["ν"], CGLP_results["ν"])
        push!(CG_params["λ"], CGLP_results["λ"])
        push!(CG_params["lp_relaxation_solution_time_taken"], round(mp_solution_end_time - mp_solution_start_time, digits = 3))

        if method == "ours"
            local negative_full_labels  
            local base_labels_time
            local full_labels_time
            try
                (negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
                    data, graph, 
                    CGLP_results["κ"], 
                    CGLP_results["μ"], 
                    CGLP_results["ν"], 
                    CGLP_results["λ"], 
                    ;
                    charge_cost_heterogenous = charge_cost_heterogenous,
                    neighborhoods = neighborhoods,
                    ngroute = ngroute,
                    elementary = elementary,
                    time_limit = time_limit - (time() - start_time),
                )
            catch e
                if isa(e, TimeLimitException)
                    break
                else
                    throw(e)
                end
            end
            generated_paths = get_paths_from_negative_path_labels(
                data, graph, negative_full_labels,
            )
            push!(
                CG_params["sp_base_time_taken"],
                round(base_labels_time, digits=3)
            )
            push!(
                CG_params["sp_full_time_taken"],
                round(full_labels_time, digits=3)
            )
            push!(
                CG_params["sp_total_time_taken"],
                round(base_labels_time + full_labels_time, digits=3)
            )
        elseif method == "benchmark"
            local negative_pure_path_labels  
            local pure_path_labels_time
            try
                (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
                    data, graph, 
                    CGLP_results["κ"], 
                    CGLP_results["μ"], 
                    CGLP_results["ν"], 
                    CGLP_results["λ"], 
                    ;
                    neighborhoods = neighborhoods, 
                    ngroute = ngroute, 
                    time_windows = time_windows,
                    elementary = elementary,
                    time_limit = time_limit - (time() - start_time),
                )
            catch e
                if isa(e, TimeLimitException)
                    break
                else
                    throw(e)
                end
            end
            generated_paths = get_paths_from_negative_pure_path_labels(
                data, graph, negative_pure_path_labels,
            )
            push!(
                CG_params["sp_base_time_taken"],
                0.0
            )
            push!(
                CG_params["sp_full_time_taken"],
                round(pure_path_labels_time, digits=3)
            )
            push!(
                CG_params["sp_total_time_taken"],
                round(pure_path_labels_time, digits=3)
            )
        end
        
        if length(generated_paths) == 0
            push!(CG_params["number_of_new_paths"], 0)
            converged = true
        else
            push!(
                CG_params["number_of_new_paths"],
                sum(length(v) for v in values(generated_paths))
            )
        end

        mp_constraint_time = add_paths_to_path_model!(
            model,
            z,
            some_paths, 
            path_costs,
            path_service,
            generated_paths,
            SR3_constraints,
            data, graph,
        )

        push!(
            CG_params["number_of_paths"], 
            sum(length(v) for v in values(some_paths))
        )
        push!(
            CG_params["lp_relaxation_constraint_time_taken"],
            mp_constraint_time,
        )
        add_message!(
            printlist, 
            @sprintf(
                "Iteration %3d | %.4e | %10d | %9.3f | %9.3f | %9.3f | %9.3f | %10d \n", 
                counter,
                CG_params["objective"][counter],
                CG_params["number_of_paths"][counter],
                CG_params["lp_relaxation_solution_time_taken"][counter],
                CG_params["sp_base_time_taken"][counter],
                CG_params["sp_full_time_taken"][counter],                
                CG_params["lp_relaxation_constraint_time_taken"][counter],
                CG_params["number_of_new_paths"][counter],
            ),
            verbose,
        )

        # cycling detection
        if (
            length(CG_params["number_of_paths"]) >= 2
            && CG_params["number_of_paths"][end] == CG_params["number_of_paths"][end-1]
            && !converged
        )
            CG_params["cycled"] = true
            add_message!(printlist, "CG cycling error.\n", verbose)
            break
        end
    end

    CG_params["converged"] = converged
    CG_params["counter"] = counter
    end_time = time() 
    time_taken = round(end_time - start_time, digits = 3)
    CG_params["time_taken"] = time_taken
    CG_params["time_limit_reached"] = (time_taken > time_limit)
    CG_params["errored"] = CGLP_results["errored"]
    CG_params["artificial"] = CGLP_results["artificial"]
    CG_params["lp_relaxation_time_taken"] = sum.(zip(CG_params["lp_relaxation_constraint_time_taken"], CG_params["lp_relaxation_solution_time_taken"]))
    CG_params["lp_relaxation_time_taken_total"] = sum(CG_params["lp_relaxation_time_taken"])
    CG_params["sp_base_time_taken_total"] = sum(CG_params["sp_base_time_taken"])
    CG_params["sp_full_time_taken_total"] = sum(CG_params["sp_full_time_taken"])
    CG_params["sp_time_taken_total"] = CG_params["sp_base_time_taken_total"] + CG_params["sp_full_time_taken_total"]
    CG_params["lp_relaxation_time_taken_mean"] = CG_params["lp_relaxation_time_taken_total"] / length(CG_params["lp_relaxation_time_taken"])
    CG_params["sp_base_time_taken_mean"] = CG_params["sp_base_time_taken_total"] / length(CG_params["sp_base_time_taken"])
    CG_params["sp_full_time_taken_mean"] = CG_params["sp_full_time_taken_total"] / length(CG_params["sp_full_time_taken"])
    CG_params["sp_time_taken_mean"] = CG_params["sp_base_time_taken_mean"] + CG_params["sp_full_time_taken_mean"]

    add_message!(
        printlist, 
        @sprintf(
            "Total         |            | %10d | %9.3f | %9.3f | %9.3f | %9.3f | \n", 
            CG_params["number_of_paths"][end],
            sum(CG_params["lp_relaxation_solution_time_taken"]),
            sum(CG_params["sp_base_time_taken"]),
            sum(CG_params["sp_full_time_taken"]),
            sum(CG_params["lp_relaxation_constraint_time_taken"]),
        ),
        verbose,
    )

    return CGLP_results, CG_params
end

function path_formulation_solve_integer_model!(
    model::Model,
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    },
    artificial_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
)
    for (key, p) in keys(z)
        JuMP.set_integer(z[key, p])
    end

    CGIP_solution_start_time = time()
    @suppress optimize!(model)
    CGIP_solution_end_time = time()
    if JuMP.termination_status(model) in [
        JuMP.INFEASIBLE,
        JuMP.INFEASIBLE_OR_UNBOUNDED,
        JuMP.TIME_LIMIT,
        JuMP.MEMORY_LIMIT,
    ]
        CGIP_results = Dict(
            "errored" => true,
            "artificial" => false,
        )
        add_message!(printlist, "CGIP model errored.", verbose)
    else
        CGIP_results = Dict(
            "errored" => false,
            "objective" => objective_value(model),
            "z" => Dict(
                (key, p) => value.(z[key, p])
                for (key, p) in keys(z)
            ),
            "time_taken" => round(CGIP_solution_end_time - CGIP_solution_start_time, digits = 3),
        )
        CGIP_results["artificial"] = any(
            value.(z[(key, p)]) > 1e-3
            for key in keys(artificial_paths)
                for p in 1:length(artificial_paths[key])
        )
    end

    for (key, p) in keys(z)
        JuMP.unset_integer(z[key, p])
    end

    return CGIP_results
end

function enumerate_violated_path_WSR3_inequalities(
    paths::Vector{Tuple{Float64, Path}},
    graph::EVRPGraph,
)
    S_list = Tuple{Float64, Vararg{Int}}[]
    for S in combinations(graph.N_customers, 3)
        S_paths = [
            (val, p) for (val, p) in values(paths)
            if !isdisjoint(p.customer_arcs, Tuple.(permutations(S, 2)))
        ]
        violation = sum([x[1] for x in S_paths], init = 0.0) - 1.0
        if violation ≤ 1e-6
            continue
        end
        # Find included charging stations if they exist
        included_charging_stations = Set{Int}()
        for (val, p) in S_paths
            customer_arcs = intersect(Tuple.(permutations(S, 2)), p.customer_arcs)
            for ca in customer_arcs
                union!(included_charging_stations, 
                    Set{Int}(
                        a1[2]
                        for (a1, a2) in zip(p.arcs[1:end-1], p.arcs[2:end])
                            if (
                                a1[1] == ca[1] 
                                && a1[2] != ca[2] 
                                && a2[2] == ca[2]
                                && a1[2] in graph.N_charging
                            )
                    )
                )
            end
        end
        push!(S_list, (violation, S..., included_charging_stations...))
    end
    sort!(S_list, by = x -> x[1], rev = true)
    return S_list
end


function detect_cycle_in_path(
    p::Path,
    node::Int,
)
    # Assumes that `node` is a customer.

    # finds first and last occurrences of node
    nodes = get_nodes(p)
    start_ind = findfirst(x -> x == node, nodes)
    end_ind = findlast(x -> x == node, nodes)
    return nodes[start_ind:end_ind]
end


function detect_cycles_in_path_solution(
    path_results::Vector{Path},
    graph::EVRPGraph,
)
    cycles = Dict{Int, Set{Int}}()
    for p in path_results
        for node in findall(x -> x > 1, p.served)
            if !(node in keys(cycles))
                cycles[node] = Set{Int}()
            end
            union!(cycles[node], detect_cycle_in_path(p, node))
        end
    end
    return cycles
end


function modify_neighborhoods_with_found_cycles!(
    neighborhoods::BitMatrix,
    cycles_lookup::Dict{Int, Set{Int}},
)
    for (node, changed_nodes) in pairs(cycles_lookup)
        neighborhoods[node, collect(changed_nodes)] .= 1
    end
end


function delete_paths_with_found_cycles_from_model!(
    model::Model,
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
            Int,
        }, 
        VariableRef,
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    path_costs::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Int},
    },
    path_service::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        Vector{Int},
    },
    cycles_lookup::Dict{
        Int,
        Set{Int},
    },
    graph::EVRPGraph,
)
    for state_pair in keys(some_paths)
        delete_inds = Int[]
        for (count, path) in enumerate(some_paths[state_pair])
            # check if path needs to be removed
            removed = false
            for (node, changed_nodes) in pairs(cycles_lookup)
                if path.served[node] ≤ 1
                    continue
                end
                cycle = detect_cycle_in_path(path, node)
                if !isempty(setdiff(intersect(cycle, changed_nodes), [node]))
                    removed = true
                    # println("""
                    # Deleted path due to $node, $changed_nodes: 
                    # $cycle; \t $(path.arcs)
                    # """)
                    break
                end
            end
            if removed
                push!(delete_inds, count)
            end
        end
        # delete variables from model
        orig_inds = 1:length(some_paths[state_pair])
        for count in delete_inds
            delete(model, z[(state_pair, count)])
            pop!(z, (state_pair, count))
        end
        # rename variable bindings in dictionary `z`
        for (i, ind) in enumerate(sort(setdiff(orig_inds, delete_inds)))
            z[(state_pair, i)] = pop!(z, (state_pair, ind))
        end
        # delete paths, their costs, and service info via `delete_inds`
        deleteat!(some_paths[state_pair], delete_inds)
        deleteat!(path_costs[state_pair], delete_inds)
        for node in graph.N_customers
            deleteat!(path_service[(state_pair, node)], delete_inds)
        end
    end
end


function enumerate_violated_path_SR3_inequalities(
    paths::Vector{Tuple{Float64, Path}},
    graph::EVRPGraph,
)
    S_list = Tuple{Float64, NTuple{3, Int}}[]
    for S in combinations(graph.N_customers, 3)
        violation = sum(
            val * floor(sum(p.served[S]) / 2)
            for (val, p) in values(paths)
        ) - 1.0
        if violation > 1e-6
            push!(S_list, (violation, Tuple(S)))
        end
    end
    sort!(S_list, by = x -> x[1], rev = true)
    return S_list
end

function select_representative_violated_path_SR3_inequalities(
    SR3_list_in::Vector{Tuple{Float64, NTuple{3, Int}}},
    data::EVRPData,
)

    function find_closest_triplet(
        clique::Tuple{Float64, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}},
        data::EVRPData,
    )
        @assert length(clique[2]) ≥ 1
        @assert length(clique[3]) ≥ 1
        @assert length(clique[4]) ≥ 1
        (_, S) = minimum([
            (
                data.q[i,j] + data.q[j,k] + data.q[k,i],
                Tuple(sort([i,j,k]))
            )
            for i in clique[2] for j in clique[3] for k in clique[4]
        ])
        return S
    end

    found_cliques = Tuple{Float64, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}, Tuple{Vararg{Int}}}[]
    SR3_list = copy(SR3_list_in)
    while length(SR3_list) != 0
        f_found = 0.0
        local found_triples
        local clique
        for (i, (f, S)) in enumerate(SR3_list)
            c1 = Tuple(vcat([S[3]],
                [setdiff(S_, S)[1] for (f_, S_) in union(SR3_list[1:i-1], SR3_list[i+1:end])
                if isapprox(f, f_, atol=1e-12) && S[1] in S_ && S[2] in S_],
            ))
            c2 = Tuple(vcat([S[2]],
                [setdiff(S_, S)[1] for (f_, S_) in union(SR3_list[1:i-1], SR3_list[i+1:end])
                if isapprox(f, f_, atol=1e-12) && S[1] in S_ && S[3] in S_],
            ))
            c3 = Tuple(vcat([S[1]],
                [setdiff(S_, S)[1] for (f_, S_) in union(SR3_list[1:i-1], SR3_list[i+1:end])
                if isapprox(f, f_, atol=1e-12) && S[2] in S_ && S[3] in S_],
            ))
            clique = (f, c1, c2, c3)
            found_triples = [
                Tuple(sort([i, j, k]))
                for i in clique[2], j in clique[3], k in clique[4]
            ]
            matches_in_SR3_list = [
                (f_, S_) for (f_, S_) in SR3_list
                    if (S_ in found_triples && isapprox(f, f_, atol=1e-12))
            ]
            if length(found_triples) == length(matches_in_SR3_list)
                f_found = f
                break
            end
        end
        if f_found == 0.0
            break
        end
        SR3_list = [
            (f_, S_) for (f_, S_) in SR3_list
                if !(S_ in found_triples && f_found ≈ f_)
        ]
        push!(found_cliques, clique)
    end
    select_triples = [
        find_closest_triplet(clique, data)
        for clique in found_cliques
    ]
    SR3_list_final = vcat(
        Tuple{Float64, NTuple{3, Int}}[
            (f, S) for (f, S) in SR3_list_in
                if S in select_triples
        ],
        SR3_list,
    )
    return sort(SR3_list_final, by = x -> x[1], rev = true)
end

function enumerate_violated_path_SRnk_inequalities(
    paths::Vector{Tuple{Float64, Path}},
    graph::EVRPGraph,
    n::Int,
    k::Int,
)
    S_list = Tuple{Float64, NTuple{n, Int}}[]
    for S in combinations(graph.N_customers, n)
        violation = sum(
            val * floor(sum(p.served[S]) / k)
            for (val, p) in values(paths)
        ) - floor(n / k)
        if violation > 1e-6
            push!(S_list, (violation, Tuple(S)))
        end
    end
    sort!(S_list, by = x -> x[1], rev = true)
    return S_list
end

function compute_memory_set_of_lmSRnk_inequality(
    paths::Vector{Tuple{Float64, Path}},
    S::Tuple{Vararg{Int}},
    k::Int,
)
    M = Set(S)
    for (_, p) in paths
        M_add = Set{Int}()
        state = 0
        nodes = get_nodes(p)
        for node in nodes
            if node in S
                state += 1
                if state ≥ k
                    union!(M, M_add)
                    M_add = Set{Int}()
                    state -= k
                end
            elseif state > 0
                union!(M_add, node)
            end
        end
    end
    return Tuple{Vararg{Int}}(sort(collect(M)))
end


function add_SR3_constraints_to_path_model!(
    model::Model, 
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    SR3_constraints::Dict{
        NTuple{3, Int},
        ConstraintRef
    }, 
    generated_SR3_list::Vector{Tuple{Float64, NTuple{3, Int}}}, 
)
    for item in generated_SR3_list
        S = item[2]
        SR3_constraints[S] = @constraint(
            model,
            sum(
                sum(
                    floor(sum(p.served[collect(S)]) / 2) * z[state_pair, p_ind]
                    for (p_ind, p) in enumerate(some_paths[state_pair])
                )
                for state_pair in keys(some_paths)
            ) ≤ 1
        )
    end
end

function update_lmSR3_constraint_list!(
    lmSR3_list::Vector{Tuple{Float64, NTuple{3, Int}, T}},
    implemented_lmSR3_list::Vector{Tuple{Float64, NTuple{3, Int}, T}},
) where {T <: Tuple{Vararg{Int}}}
    lmSR3_list_todelete = Tuple{Float64, NTuple{3, Int}, Tuple{Vararg{Int}}}[]
    lmSR3_list_toadd = Tuple{Float64, NTuple{3, Int}, Tuple{Vararg{Int}}}[]
    for item in implemented_lmSR3_list
        S, M = item[2], item[3]
        big_M = collect(M)
        lmSR3_list_todeleteinds = Int[]
        for (i, item_) in enumerate(lmSR3_list)
            S_, M_ = item_[2], item_[3]
            if S == S_
                append!(big_M, collect(M_))
                push!(lmSR3_list_todeleteinds, i)
                push!(lmSR3_list_todelete, item_)
            end
        end
        deleteat!(lmSR3_list, lmSR3_list_todeleteinds)
        big_M = Tuple(sort(unique(big_M)))
        push!(lmSR3_list_toadd, (item[1], S, big_M))
    end
    append!(lmSR3_list, lmSR3_list_toadd)
    return lmSR3_list_todelete, lmSR3_list_toadd
end

function add_lmSR3_constraints_to_path_model!(
    model::Model, 
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    lmSR3_constraints::Dict{
        Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
        ConstraintRef,
    }, 
    lmSR3_list_todelete::Vector{
        Tuple{Float64, NTuple{3, Int}, T}
    }, 
    lmSR3_list_toadd::Vector{
        Tuple{Float64, NTuple{3, Int}, T}
    }, 
) where {T <: Tuple{Vararg{Int}}}
    for item in lmSR3_list_todelete
        S, M = item[2], item[3]
        delete(model, lmSR3_constraints[(S, M)])
        pop!(lmSR3_constraints, (S, M))
    end
    for item in lmSR3_list_toadd
        S, M = item[2], item[3]
        lmSR3_constraints[(S, M)] = @constraint(
            model,
            sum(
                sum(
                    compute_path_coefficient_in_lmSRnk_constraint(p, S, M, 2) * z[state_pair, p_ind]
                    for (p_ind, p) in enumerate(some_paths[state_pair])
                )
                for state_pair in keys(some_paths)
            ) ≤ 1
        )
    end
    return
end

function add_lmSR3_constraints_to_path_model!(
    model::Model, 
    z::Dict{
        Tuple{
            Tuple{NTuple{3, Int}, NTuple{3, Int}},
            Int,
        },
        VariableRef,
    },
    some_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    lmSR3_constraints::Dict{
        Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
        ConstraintRef,
    }, 
    lmSR3_list_toadd::Vector{
        Tuple{Float64, NTuple{3, Int}, T}
    }, 
) where {T <: Tuple{Vararg{Int}}}
    for item in lmSR3_list_toadd
        S, M = item[2], item[3]
        lmSR3_constraints[(S, M)] = @constraint(
            model,
            sum(
                sum(
                    compute_path_coefficient_in_lmSRnk_constraint(p, S, M, 2) * z[state_pair, p_ind]
                    for (p_ind, p) in enumerate(some_paths[state_pair])
                )
                for state_pair in keys(some_paths)
            ) ≤ 1
        )
    end
    return
end


function path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data::EVRPData, 
    graph::EVRPGraph,
    ;
    Env = nothing,
    method::String = "ours",
    charge_cost_heterogenous::Bool = false,
    time_windows::Bool = false,
    elementary::Bool = false,
    ngroute::Bool = true,
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size::String = "small",
    ngroute_neighborhood_charging_size::String = "small",
    verbose::Bool = true,
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
    use_adaptive_ngroute::Bool = true,
    use_SR3_cuts::Bool = true,
    use_lmSR3_cuts::Bool = true,
    max_SR3_cuts::Int = 5, 
    randomize_cuts::Bool = false,
)
    start_time = time()

    if ngroute && isnothing(neighborhoods)
        neighborhoods = compute_ngroute_neighborhoods(
            graph,
            ngroute_neighborhood_size; 
            depots_size = ngroute_neighborhood_depots_size,
            charging_size = ngroute_neighborhood_charging_size,
        )
    end

    artificial_paths = generate_artificial_paths(data, graph)
    some_paths = deepcopy(artificial_paths)
    path_costs = compute_path_costs(
        data, graph, 
        some_paths,
    )
    path_service = compute_path_service(
        graph,
        some_paths,
    )

    printlist = String[]

    start_printlist = String[]
    push!(start_printlist, 
        @sprintf(
            """
            Starting column generation on the path formulation.
            # customers:                    %3d
            # depots:                       %3d
            # charging stations:            %3d
            # vehicles:                     %3d
            time windows?:                  %s
            heterogenous charging costs?:   %s
            charge_cost_nlevels:            %3d

            method:                         %s
            """,
            graph.n_customers,
            graph.n_depots,
            graph.n_charging,
            data.n_vehicles,
            time_windows,
            charge_cost_heterogenous,
            data.charge_cost_nlevels,
            method,
        )
    )
    if ngroute
        push!(start_printlist, 
            @sprintf(
                """
                ngroute:                        %s
                ngroute neighborhood size:
                    customers                   %3d
                    depots                      %s
                    charging                    %s
                """,
                ngroute,
                ngroute_neighborhood_size,
                ngroute_neighborhood_depots_size,
                ngroute_neighborhood_charging_size,
            )
        )
        push!(start_printlist,
            @sprintf(
                "use_adaptive_ngroute:           %s\n",
                use_adaptive_ngroute
            )
        )
        if use_SR3_cuts
            push!(start_printlist,
                @sprintf(
                    """
                    use_SR3_cuts:                   %s
                        use_lmSR3_cuts:             %s
                        max_SR3_cuts:               %3d
                    """,
                    use_SR3_cuts,
                    use_lmSR3_cuts,
                    max_SR3_cuts,
                )
            )
        else
            push!(start_printlist,
                @sprintf(
                    """
                    use_SR3_cuts:                   %s
                    """,
                    use_SR3_cuts,
                )
            )
        end
    else
        push!(start_printlist,
            @sprintf(
                """
                ngroute:                        %s
                elementary:                     %s
                """,
                ngroute,
                elementary,
            )
        )
    end

    for message in start_printlist
        add_message!(printlist, message, verbose)
    end

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

    model, z = path_formulation_build_model(
        data, graph, some_paths, path_costs, path_service,
        ; 
        Env = Env,
    )

    local CGLP_results
    local CG_params
    local CGIP_results
    local converged = false

    if use_lmSR3_cuts
        SR3_constraints = Dict{
            Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
            ConstraintRef,
        }()
        SR3_list = Tuple{Float64, NTuple{3, Int}, Tuple{Vararg{Int}}}[]
    else
        SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
        SR3_list = Tuple{Float64, NTuple{3, Int}}[]
    end
    CGLP_all_results = Dict[]
    CGIP_all_results = Dict[]
    all_params = Dict[]
    CG_all_params = Dict[]
    if ngroute
        CG_all_neighborhoods = BitMatrix[]
    else
        CG_all_neighborhoods = nothing
    end

    while true
        continue_flag = false
        iteration_params = Dict{String, Any}()

        CGLP_results, CG_params = path_formulation_column_generation!(
            model, z, SR3_constraints,
            data, graph,
            artificial_paths,
            some_paths, path_costs, path_service,
            printlist,
            ;
            method = method,
            charge_cost_heterogenous = charge_cost_heterogenous,
            time_windows = time_windows,
            elementary = elementary,
            neighborhoods = neighborhoods,
            ngroute = ngroute,
            verbose = verbose,
            time_limit = time_limit - (time() - start_time),
            max_iters = max_iters,
        )
        CGIP_results = path_formulation_solve_integer_model!(
            model,
            z,
            artificial_paths,
        )
        push!(CGLP_all_results, CGLP_results)
        push!(CGIP_all_results, CGIP_results)
        push!(CG_all_params, CG_params)
        if ngroute
            push!(CG_all_neighborhoods, copy(neighborhoods))
        end

        if CGLP_results["errored"] || CGIP_results["errored"] || CG_params["cycled"]
            add_message!(printlist, "Terminating column generation...\n", verbose)
            iteration_params["errored"] = true
            break
        else
            iteration_params["errored"] = false
        end

        # Termination criteria
        CG_params["CGLP_objective"] = CGLP_results["objective"]
        CG_params["CGIP_objective"] = CGIP_results["objective"]
        CG_params["LP_IP_gap"] = 1.0 - CGLP_results["objective"] / CGIP_results["objective"]
        for message in [
            @sprintf("\n"),
            @sprintf("Time taken (s):       %9.3f s\n", CG_params["time_taken"]),
            @sprintf("(CGLP) Objective:         %.4e\n", CGLP_results["objective"]),
            @sprintf("(CGIP) Objective:         %.4e\n", CGIP_results["objective"]),
            @sprintf("%% gap:                %9.3f %%\n", CG_params["LP_IP_gap"] * 100.0),
        ]
            add_message!(printlist, message, verbose)
        end

        if CG_params["converged"]
            CGLP_results["paths"] = collect_path_solution_support(
                CGLP_results, some_paths, data, graph
            )
            CGIP_results["paths"] = collect_path_solution_support(
                CGIP_results, some_paths, data, graph
            )
        end
        
        iteration_params["CGLP_objective"] = CG_params["CGLP_objective"]
        iteration_params["CGIP_objective"] = CG_params["CGIP_objective"]
        iteration_params["CGLP_artificial"] = CGLP_results["artificial"]
        iteration_params["CGIP_artificial"] = CGIP_results["artificial"]
        iteration_params["CG_LP_IP_gap"] = CG_params["LP_IP_gap"]
        iteration_params["CG_time_taken"] = CG_params["time_taken"]
        iteration_params["CG_sp_time_taken_mean"] = CG_params["sp_time_taken_mean"]
        iteration_params["method"] = "none"
        iteration_params["ngroute_neighborhood_size"] = ngroute ? (graph.n_customers + graph.n_charging) * mean(neighborhoods[graph.N_customers, vcat(graph.N_customers, graph.N_charging)]) : 0
        iteration_params["cycles_lookup_length"] = 0
        iteration_params["implemented_SR3_cuts_count"] = 0
        iteration_params["converged"] = false
        iteration_params["time_limit_reached"] = false

        # check if converged
        if CGIP_results["objective"] ≈ CGLP_results["objective"]
            converged = true
            iteration_params["converged"] = true
        end

        if CG_params["converged"] && !continue_flag && !converged && ngroute && use_adaptive_ngroute
            # see if any path in solution was non-elementary
            cycles_lookup = detect_cycles_in_path_solution([p for (val, p) in CGLP_results["paths"]], graph)
            if length(cycles_lookup) > 0 
                delete_paths_with_found_cycles_from_model!(model, z, some_paths, path_costs, path_service, cycles_lookup, graph)
                modify_neighborhoods_with_found_cycles!(neighborhoods, cycles_lookup)
                continue_flag = true
                iteration_params["method"] = "use_adaptive_ngroute"
                iteration_params["cycles_lookup_length"] = length(cycles_lookup)
                add_message!(printlist, "Expanded ng-route neighborhoods by $(length(cycles_lookup))\n", verbose)
                add_message!(printlist, "\n", verbose)
            end
        end

        if CG_params["converged"] && !continue_flag && !converged && ngroute && use_SR3_cuts
            # if no path in solution was non-elementary, 
            # generate violated SR3 inequalities
            generated_SR3_list = enumerate_violated_path_SR3_inequalities(
                CGLP_results["paths"],
                graph,
            )
            add_message!(printlist, "Found SR3 cuts:\t\t\t$(length(generated_SR3_list))\n", verbose)
            if length(generated_SR3_list) != 0
                generated_SR3_list = select_representative_violated_path_SR3_inequalities(
                    generated_SR3_list, data,
                )
                add_message!(printlist, "Selected SR3 cuts:\t\t$(length(generated_SR3_list))\n", verbose)
            end

            # sample cuts if too many
            if length(generated_SR3_list) ≤ max_SR3_cuts
                implemented_SR3_list = generated_SR3_list
            elseif randomize_cuts
                implemented_SR3_list = sample(
                    generated_SR3_list, 
                    Weights([val for (val, _) in generated_SR3_list]),
                    max_SR3_cuts, 
                    replace = false,
                )
            else
                implemented_SR3_list = generated_SR3_list[1:max_SR3_cuts]
            end
            add_message!(printlist, "Sampled SR3 cuts:\t\t$(length(implemented_SR3_list))\n", verbose)
            
            if length(implemented_SR3_list) != 0
                if use_lmSR3_cuts
                    implemented_SR3_list = Tuple{Float64, Tuple{Int64, Int64, Int64}, Tuple{Vararg{Int64}}}[
                        (val, S, compute_memory_set_of_lmSRnk_inequality(CGLP_results["paths"], S, 2))
                        for (val, S) in implemented_SR3_list
                    ]
                    SR3_list_todelete, SR3_list_toadd = update_lmSR3_constraint_list!(
                        SR3_list, implemented_SR3_list,
                    )
                    add_lmSR3_constraints_to_path_model!(
                        model, z, some_paths, 
                        SR3_constraints, SR3_list_todelete, SR3_list_toadd,
                    )
                    # append!(SR3_list, implemented_SR3_list)
                    # add_lmSR3_constraints_to_path_model!(
                    #     model, z, some_paths, 
                    #     SR3_constraints, implemented_SR3_list,
                    # )
                    iteration_params["method"] = "use_lmSR3_cuts"
                    iteration_params["implemented_SR3_cuts_count"] = length(implemented_SR3_list)
                    continue_flag = true
                    add_message!(printlist, "Imposed lm-SR3 cuts:\t\t$(length(implemented_SR3_list))\n", verbose)
                    for (val, S, M) in implemented_SR3_list
                        add_message!(printlist, "$S, $M: $val\n", verbose)
                    end
                    add_message!(printlist, "\n", verbose)
                else
                    append!(SR3_list, implemented_SR3_list)
                    # Add violated inequalities to master problem
                    add_SR3_constraints_to_path_model!(
                        model, z, some_paths, 
                        SR3_constraints, implemented_SR3_list, 
                    )
                    iteration_params["method"] = "use_SR3_cuts"
                    iteration_params["implemented_SR3_cuts_count"] = length(implemented_SR3_list)
                    continue_flag = true
                    add_message!(printlist, "Imposed SR3 cuts:\t\t$(length(implemented_SR3_list))\n", verbose)
                    for (val, S) in implemented_SR3_list
                        add_message!(printlist, "$S: $val\n", verbose)
                    end
                    add_message!(printlist, "\n", verbose)
                end
            end
        end

        push!(all_params, iteration_params)

        add_message!(
            printlist,
            @sprintf("Total time taken (s): %9.3f s\n\n", time() - start_time),
            verbose,
        )

        if !(time_limit > time() - start_time)
            iteration_params["time_limit_reached"] = true
            break
        end
        
        if !continue_flag
            break
        end
    end

    return (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    )
end


function collect_path_solution_support(
    results, 
    paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    },
    data::EVRPData,
    graph::EVRPGraph,
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

    sort!(
        results_paths,
        by = x -> (-compute_path_cost(data, graph, x[2]), -x[1]),
    )
    return results_paths
end

function collect_solution_metrics!(
    results, 
    data::EVRPData,
)
    if !("paths" in keys(results))
        error()
    end

    total_subpath_length = 0.0
    weighted_total_subpath_length = 0.0
    total_subpath_ncust = 0.0
    weighted_total_subpath_ncust = 0.0
    num_subpaths = 0.0
    weighted_num_subpaths = 0.0
    total_path_length = 0.0
    weighted_total_path_length = 0.0
    total_path_ncust = 0.0
    weighted_total_path_ncust = 0.0
    num_paths = 0.0
    weighted_num_paths = 0.0
    total_ps_length = 0.0
    weighted_total_ps_length = 0.0
    utilization_total = 0.0
    driving_time_total = 0.0
    charging_time_total = 0.0
    for (val, p) in results["paths"]
        if (
            length(p.subpaths) == 1 
            && (
                p.subpaths[1].artificial # artificial path
                || length(p.subpaths[1].arcs) == 1 # path from depot to depot
            )
        )
            continue
        end
        total_subpath_length += sum(sum(s.served) + 1 for s in p.subpaths) 
        weighted_total_subpath_length += val * sum(sum(s.served) + 1 for s in p.subpaths)
        total_subpath_ncust += sum(sum(s.served) for s in p.subpaths)
        weighted_total_subpath_ncust += val * sum(sum(s.served) for s in p.subpaths)

        num_subpaths += length(p.subpaths)
        weighted_num_subpaths += val * length(p.subpaths)
        
        total_path_length += sum(p.served) + length(p.subpaths)
        weighted_total_path_length += val * (sum(p.served) + length(p.subpaths))
        total_path_ncust += sum(p.served)
        weighted_total_path_ncust += val * sum(p.served)

        num_paths += 1
        weighted_num_paths += val

        total_ps_length += length(p.subpaths)
        weighted_total_ps_length += val * length(p.subpaths)

        utilization_total += val * p.subpaths[end].current_time
        driving_time_total += val * sum(s.current_time - s.starting_time for s in p.subpaths)
        charging_time_total += val * sum([a.delta for a in p.charging_arcs], init = 0.0)
    end

    results["mean_subpath_length"] = total_subpath_length / num_subpaths
    results["weighted_mean_subpath_length"] = weighted_total_subpath_length / weighted_num_subpaths
    results["mean_subpath_ncust"] = total_subpath_ncust / num_subpaths
    results["weighted_mean_subpath_ncust"] = weighted_total_subpath_ncust / weighted_num_subpaths
    
    results["mean_path_length"] = total_path_length / num_paths
    results["weighted_mean_path_length"] = weighted_total_path_length / weighted_num_paths
    results["mean_path_ncust"] = total_path_ncust / num_paths
    results["weighted_mean_path_ncust"] = weighted_total_path_ncust / weighted_num_paths

    results["mean_ps_length"] = total_ps_length / num_paths
    results["weighted_mean_ps_length"] = weighted_total_ps_length / weighted_num_paths

    results["utilization"] = utilization_total / (weighted_num_paths * data.T)
    results["driving_time_proportion"] = driving_time_total / (weighted_num_paths * data.T)
    results["charging_time_proportion"] = charging_time_total / (weighted_num_paths * data.T)

    return results

end

function collect_path_solution_metrics!(
    results,
    data::EVRPData, 
    graph::EVRPGraph,
    paths,
)
    results["paths"] = collect_path_solution_support(results, paths, data, graph)
    collect_solution_metrics!(results, data)
    return results
end

function compute_objective_from_path_solution(
    results_paths::Vector{Tuple{Float64, Path}},
    data::EVRPData,
    graph::EVRPGraph,
)
    return sum(
        [val * compute_path_cost(data, graph, p)
        for (val, p) in results_paths],
        init = 0.0,
    )
end

function plot_path_solution(
    results,
    data::EVRPData,
    graph::EVRPGraph,
    paths,
)
    results["paths"] = collect_path_solution_support(results, paths, data, graph)
    return plot_solution(results["paths"], data)
end


function plot_solution(
    results_paths::Vector{Tuple{Float64, Path}}, 
    data::EVRPData,
)
    p = plot_instance(data)
    
    n_paths = length(results_paths) 
    colors = get(ColorSchemes.tol_bright, collect(0:1:(n_paths-1)) ./(n_paths-1))

    all_plots = []
    for (i, (val, path)) in enumerate(results_paths)
        p = plot_instance(data)
        arcs = vcat(collect(s.arcs for s in path.subpaths)...)
        for (j, arc) in enumerate(arcs)
            Plots.plot!(
                data.coords[1,collect(arc)],
                data.coords[2,collect(arc)],
                color = colors[i],
                alpha = 0.5,
                lw = 1,
            )
        end
        Plots.plot!(title = "Vehicle $i: $val")
        push!(all_plots, p)
    end

    P = Plots.plot(
        all_plots..., 
        layout = (n_paths, 1), 
        size = (500, 400 * n_paths), 
        legend = false,
        fmt = :png,
    )
    return P
end

function plot_solution(
    results::Dict{String, Any},
    data::EVRPData,
)
    return plot_solution(results["paths"], data)
end