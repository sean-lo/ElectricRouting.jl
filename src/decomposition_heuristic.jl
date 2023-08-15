using JuMP
using Gurobi
using Suppressor
using DataStructures
using Printf

include("subpath_stitching.jl")
include("path_formulation.jl")
include("utils.jl")

function find_nondominated_paths_nocharge(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_windows::Bool = false,
    christofides::Bool = false,
    single_service::Bool = false,
    check_customers::Bool = false,
)

    modified_costs = compute_arc_modified_costs(graph, data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in union(graph.N_depots, graph.N_customers)
        )
        for starting_node in graph.N_depots
    )
    unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
    for node in graph.N_depots
        served = zeros(Int, graph.n_customers)
        if check_customers
            key = (0, served...)
        else
            key = (0,)
        end
        base_labels[node][node][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(unexplored_states, (key..., node, node))
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        if !(current_key in keys(base_labels[starting_node][current_node]))
            continue
        end
        current_subpath = base_labels[starting_node][current_node][current_key]
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging_extra,
        )
            if single_service
                if next_node in graph.N_customers && current_subpath.served[next_node] > 0
                    continue
                end
            end
            # Preventing customer 2-cycles (Christofides)
            if christofides && next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            # time feasibility
            new_time_taken = current_subpath.time_taken + graph.t[current_node, next_node]
            if time_windows
                new_time_taken = max(
                    new_time_taken,
                    graph.α[next_node]
                )
                if new_time_taken > graph.β[next_node]
                    continue
                end
            else
                if new_time_taken + graph.min_t[next_node] > T_heuristic
                    continue
                end
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken = new_time_taken
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in graph.N_customers
                new_subpath.served[next_node] += 1
            end

            if check_customers
                new_key = (new_subpath.time_taken, new_subpath.served...)
            else
                new_key = (new_subpath.time_taken,)
            end
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for depot in graph.N_depots
        for path in values(base_labels[depot][depot])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_windows::Bool = false,
    christofides::Bool = false,
)

    modified_costs = compute_arc_modified_costs(graph, data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                NTuple{graph.n_customers + graph.n_depots, Int}, 
                SortedDict{
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            }()
            for current_node in union(graph.N_depots, graph.N_customers)
        )
        for starting_node in graph.N_depots
    )
    unexplored_states = SortedSet{NTuple{graph.n_customers + graph.n_depots + 3, Int}}()
    for node in graph.N_depots
        node_labels = zeros(Int, graph.n_customers + graph.n_depots)
        node_labels[node] = 1
        key = (0,)
        base_labels[node][node][(node_labels...)] = SortedDict{
            NTuple{1, Int}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            key => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
            )
        )
        push!(
            unexplored_states, 
            (
                key...,
                node_labels...,
                node, # starting_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[2:end-2]
        current_key = state[1:1]
        current_subpath = get(base_labels[starting_node][current_node][current_set], current_key, nothing)
        isnothing(current_subpath) && continue
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging_extra,
        )
            # Preventing customer 2-cycles (Christofides) # FIXME
            if christofides && next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            (feasible, new_set) = ngroute_check_create_set(
                graph.N_customers, neighborhoods, current_set, next_node,
            )
            !feasible && continue

            # time feasibility
            new_time_taken = current_subpath.time_taken + graph.t[current_node, next_node]
            if time_windows
                new_time_taken = max(
                    new_time_taken,
                    graph.α[next_node]
                )
                if new_time_taken > graph.β[next_node]
                    continue
                end
            else
                if new_time_taken + graph.min_t[next_node] > T_heuristic
                    continue
                end
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken = new_time_taken
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in graph.N_customers
                new_subpath.served[next_node] += 1
            end

            if !(new_set in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][new_set] = SortedDict{
                    NTuple{1, Int}, 
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][new_set],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for depot in graph.N_depots
        for set in keys(base_labels[depot][depot])
            for path in values(base_labels[depot][depot][set])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                end
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - μ[end_node]
                end
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute_alt(
    data::EVRPData,
    graph::EVRPGraph,
    neighborhoods::BitMatrix,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_windows::Bool = false,
    christofides::Bool = false,
)

    modified_costs = compute_arc_modified_costs(graph, data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                NTuple{graph.n_customers + graph.n_depots + 1, Int}, 
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in union(graph.N_depots, graph.N_customers)
        )
        for starting_node in graph.N_depots
    )
    unexplored_states = SortedSet{NTuple{graph.n_customers + graph.n_depots + 3, Int}}()
    for node in graph.N_depots
        node_labels = zeros(Int, graph.n_customers + graph.n_depots)
        node_labels[node] = 1
        key = (0,)
        base_labels[node][node][(key..., node_labels...)] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
        push!(
            unexplored_states, 
            (
                key..., 
                node_labels..., 
                node, # starting_node
                node, # current_node
            )
        )
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_set = state[2:end-2]
        current_key = state[1:end-2]
        current_subpath = get(base_labels[starting_node][current_node], current_key, nothing)
        isnothing(current_subpath) && continue
        for next_node in setdiff(
            outneighbors(graph.G, current_node), 
            current_node, 
            graph.N_charging_extra,
        )
            # Preventing customer 2-cycles (Christofides) # FIXME
            if christofides && next_node in graph.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            (feasible, new_set) = ngroute_check_create_set(
                graph.N_customers, neighborhoods, current_set, next_node,
            )
            !feasible && continue

            # time feasibility
            new_time_taken = current_subpath.time_taken + graph.t[current_node, next_node]
            if time_windows
                new_time_taken = max(
                    new_time_taken,
                    graph.α[next_node]
                )
                if new_time_taken > graph.β[next_node]
                    continue
                end
            else
                if new_time_taken + graph.min_t[next_node] > T_heuristic
                    continue
                end
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken = new_time_taken
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in graph.N_customers
                new_subpath.served[next_node] += 1
            end

            new_key = (new_subpath.time_taken,)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node],
                (new_key..., new_set...), new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in graph.N_customers
                new_state = (new_key..., new_set..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for depot in graph.N_depots
        for path in values(base_labels[depot][depot])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in graph.N_depots
        for starting_node in graph.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end


function get_negative_base_labels_from_base_labels(
    graph::EVRPGraph,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            SortedDict{
                Tuple{Vararg{Int}},
                BaseSubpathLabel,
                Base.Order.ForwardOrdering,
            },
        },
    },
)
    return BaseSubpathLabel[
        base_label
        for starting_node in graph.N_depots
            for end_node in graph.N_depots
                for (key, base_label) in base_labels[starting_node][end_node]
                    if base_label.cost < -1e-6
    ]
end

function get_negative_base_labels_from_base_labels_ngroute(
    graph::EVRPGraph,
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                Tuple{Vararg{Int}},
                SortedDict{
                    Tuple{Vararg{Int}},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                },
            },
        },
    },
)
    return BaseSubpathLabel[
        base_label
        for starting_node in graph.N_depots
            for end_node in graph.N_depots
                for set in keys(base_labels[starting_node][end_node])
                    for (key, base_label) in base_labels[starting_node][end_node][set]
                        if base_label.cost < -1e-6
    ]
end

function subproblem_iteration_nocharge(
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    time_windows::Bool = false,
    single_service::Bool = false,
    check_customers::Bool = false,
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
)
    if ngroute && !ngroute_alt
        base_labels_result = @timed find_nondominated_paths_nocharge_ngroute(
            graph, neighborhoods, κ, μ, ν, T_heuristic,
            ;
            time_windows = time_windows,
            christofides = christofides,
        )
    elseif ngroute && ngroute_alt
        base_labels_result = @timed find_nondominated_paths_nocharge_ngroute_alt(
            graph, neighborhoods, κ, μ, ν, T_heuristic,
            ;
            time_windows = time_windows,
            christofides = christofides,
        )
    else
        base_labels_result = @timed find_nondominated_paths_nocharge(
            graph, κ, μ, ν, T_heuristic,
            ;
            time_windows = time_windows,
            christofides = christofides,
            single_service = single_service,
            check_customers = check_customers,
        )
    end
    if ngroute && !ngroute_alt
        negative_base_labels = get_negative_base_labels_from_base_labels_ngroute(graph, base_labels_result.value)
    elseif (ngroute && ngroute_alt) || !ngroute
        negative_base_labels = get_negative_base_labels_from_base_labels(graph, base_labels_result.value)
    end
    return (negative_base_labels, length(negative_base_labels), base_labels_result.time)
end

function get_paths_from_negative_base_labels(
    graph::EVRPGraph,
    base_labels::Vector{BaseSubpathLabel},
)
    generated_paths = Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}},
        Vector{Path},
    }()
    for base_label in base_labels
        p = Path(
            subpaths = [
                Subpath(
                    n_customers = graph.n_customers,
                    starting_node = base_label.nodes[1],
                    starting_time = 0,
                    starting_charge = graph.B,
                    current_node = base_label.nodes[end],
                    arcs = collect(zip(base_label.nodes[1:end-1], base_label.nodes[2:end])),
                    current_time = base_label.time_taken,
                    current_charge = graph.B,
                    served = base_label.served,
                )
            ],
            charging_arcs = ChargingArc[],
            served = base_label.served,
        )
        add_path_to_generated_paths!(generated_paths, p)
    end
    return generated_paths
end

function path_formulation_column_generation_nocharge(
    data::EVRPData, 
    graph::EVRPGraph,
    ;
    time_heuristic_slack::Float64 = 0.9,
    Env = nothing,
    time_windows::Bool = false,
    path_single_service::Bool = false,
    path_check_customers::Bool = false,
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(graph.n_customers))),
    verbose::Bool = true,
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
)

    start_time = time()

    if ngroute
        neighborhoods = compute_ngroute_neighborhoods(
            graph,
            ngroute_neighborhood_size;
        )
    else
        neighborhoods = nothing
    end
    T_heuristic = Int(round(time_heuristic_slack * (graph.T + graph.B) * graph.μ / (1 + graph.μ)))

    some_paths = generate_artificial_paths(data, graph)
    path_costs = compute_path_costs(
        data, graph,
        some_paths,
    )
    path_service = compute_path_service(
        graph,
        some_paths,
    )
    CGLP_results = Dict{String, Any}()
    CG_params = Dict{String, Any}()
    CG_params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
    CG_params["objective"] = Float64[]
    CG_params["κ"] = Dict{Int, Float64}[]
    CG_params["μ"] = Dict{Int, Float64}[]
    CG_params["ν"] = Vector{Float64}[]
    CG_params["σ"] = Dict{Tuple{Vararg{Int}}, Float64}[]
    CG_params["lp_relaxation_solution_time_taken"] = Float64[]
    CG_params["sp_base_time_taken"] = Float64[]
    CG_params["sp_full_time_taken"] = Float64[]
    CG_params["sp_total_time_taken"] = Float64[]
    CG_params["lp_relaxation_constraint_time_taken"] = Float64[]
    CG_params["number_of_new_paths"] = Int[]

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

            path_single_service:            %s
            path_check_customers:           %s
            christofides:                   %s
            ngroute:                        %s
            ngroute_alt:                    %s
            ngroute neighborhood size:
                customers                   %2d
                charging / depots           %s
            
            """,
            graph.n_customers,
            graph.n_depots,
            graph.n_charging,
            data.n_vehicles,
            time_windows,
            path_single_service,
            path_check_customers,
            christofides,
            ngroute,
            ngroute_alt,
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

    model, z = path_formulation_build_model(
        data, graph, some_paths, path_costs, path_service,
        ;
        Env = Env,
    )

    while (
        !converged
        && time_limit ≥ (time() - start_time)
        && max_iters > counter
    )
        counter += 1
        mp_solution_start_time = time()
        @suppress optimize!(model)
        mp_solution_end_time = time()
        CGLP_results = Dict(
            "objective" => objective_value(model),
            "z" => Dict(
                (key, p) => value.(z[(key, p)])
                for (key, p) in keys(z)
            ),
            "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
            "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
            "ν" => dual.(model[:ν]).data,
            "σ" => Dict{Tuple{Vararg{Int}}, Float64}(
                S => dual(WSR3_constraints[S])
                for S in keys(WSR3_constraints)
            )
        )
        push!(CG_params["objective"], CGLP_results["objective"])
        push!(CG_params["κ"], CGLP_results["κ"])
        push!(CG_params["μ"], CGLP_results["μ"])
        push!(CG_params["ν"], CGLP_results["ν"])
        push!(CG_params["σ"], CGLP_results["σ"])
        push!(CG_params["lp_relaxation_solution_time_taken"], round(mp_solution_end_time - mp_solution_start_time, digits = 3))


        (negative_base_labels, _, base_labels_time) = subproblem_iteration_nocharge(
            graph, CGLP_results["κ"], CGLP_results["μ"], CGLP_results["ν"], T_heuristic,
            ;
            neighborhoods = neighborhoods,
            time_windows = time_windows,
            single_service = path_single_service,
            check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
        )
        (generated_paths) = get_paths_from_negative_base_labels(
            graph, negative_base_labels,
        )
        push!(
            CG_params["sp_base_time_taken"],
            round(base_labels_time, digits=3)
        )
        push!(
            CG_params["sp_full_time_taken"],
            0.0
        )
        push!(
            CG_params["sp_total_time_taken"],
            round(base_labels_time, digits=3)
        )
        
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
    end

    CG_params["converged"] = converged
    CG_params["counter"] = counter
    end_time = time() 
    time_taken = round(end_time - start_time, digits = 3)
    CG_params["time_taken"] = time_taken
    CG_params["time_limit_reached"] = (time_taken > time_limit)
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

    CGIP_results = path_formulation_solve_integer_model!(
        model, z
    )

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
    return CGLP_results, CGIP_results, CG_params, printlist, some_paths, model, z
end


function get_postcharge_shortest_pure_path_label(
    data::EVRPData, 
    graph::EVRPGraph,
    nodelist::Vector{Int},
    ;
    time_windows::Bool = false,
)
    modified_costs = compute_arc_modified_costs(graph, data, zeros(Float64, graph.n_customers))
    t = graph.t
    B = graph.B
    q = graph.q
    if time_windows
        α = graph.α
        β = graph.β
    else
        α = zeros(Int, graph.n_nodes_extra)
        β = repeat([graph.T], graph.n_nodes_extra)
    end

    pure_path_labels = Dict(
        (nodelist[1:j]...,) => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                PurePathLabel,
                Base.Order.ForwardOrdering,
            }(Base.Order.ForwardOrdering())
            for current_node in graph.N_nodes_extra
        )
        for j in eachindex(nodelist)
    )
    key = (0, -graph.B, graph.B)
    starting_node = nodelist[1]
    nodeseq = (starting_node,)
    pure_path_labels[nodeseq][starting_node][key] = PurePathLabel(
        0.0, 
        [starting_node],
        Int[],
        Int[],
        0,
        0,
        graph.B,
        graph.B,
        false,
        zeros(Int, graph.n_customers),
        false,
    )
    unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
    push!(unexplored_states, (key..., starting_node, nodeseq...))
    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        current_node = state[4]
        current_nodeseq = state[5:end]
        current_key = state[1:3]
        if !(current_key in keys(pure_path_labels[current_nodeseq][current_node]))
            continue
        end
        current_path = pure_path_labels[current_nodeseq][current_node][current_key]
        # println("current_path: $(current_path.nodes)")
        eventual_next_node = nodelist[length(current_nodeseq) + 1]
        for next_node in setdiff(
            vcat(eventual_next_node, graph.N_charging_extra), 
            current_node,
        )
            if !(next_node in outneighbors(graph.G, current_node))
                continue
            end
            
            (feasible, new_path) = compute_new_pure_path(
                current_path,
                current_node,
                next_node, 
                data, graph,
                α, β, modified_costs,
            )
            if !feasible
                continue 
            end
            # add new_path to collection
            if next_node in union(graph.N_customers, graph.N_depots)
                new_nodeseq = (current_nodeseq..., next_node)
            else
                new_nodeseq = current_nodeseq
            end
            new_key = (
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge,
            )
            # println("$next_node, $new_nodeseq, $new_key")
            added = add_pure_path_label_to_collection!(
                pure_path_labels[new_nodeseq][next_node], 
                new_key, new_path, 
                ;
                verbose = false,
            )
            if added && !(next_node in graph.N_depots)
                new_state = (new_key..., next_node, new_nodeseq...,)
                push!(unexplored_states, new_state)
            end
        end
    end
    if length(pure_path_labels[(nodelist...)][nodelist[end]]) == 0
        return (false, nothing)
    else
        return (true, first(values(pure_path_labels[(nodelist...)][nodelist[end]])))
    end
end

function path_formulation_decomposition_heuristic(
    data::EVRPData,
    graph::EVRPGraph,
    ;
    Env = nothing,
    time_windows::Bool = false,
    path_single_service::Bool = false,
    path_check_customers::Bool = false,
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(graph.n_customers))),
    use_integer_paths::Bool = false,
    verbose::Bool = true,
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
)
    start_time = time()
    time_heuristic_slack = 1.0
    results_paths_withcharge = Tuple{Float64, Path}[]
    while true
        feasible = true
        CGLP_results, CGIP_results, CG_params, printlist, some_paths = path_formulation_column_generation_nocharge(
            data, graph,
            ;
            time_heuristic_slack = time_heuristic_slack,
            Env = Env,
            time_windows = time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            ngroute_neighborhood_size = ngroute_neighborhood_size,
            verbose = verbose,
            time_limit = time_limit - (time() - start_time),
            max_iters = max_iters,
        )
        if use_integer_paths
            results_paths = collect_path_solution_support(CGIP_results, some_paths, data, graph)
        else
            results_paths = collect_path_solution_support(CGLP_results, some_paths, data, graph)
        end
        results_paths_withcharge = Tuple{Float64, Path}[]
        for (val, p) in results_paths
            if p.subpaths[1].artificial 
                error("Infeasible solution.")
            end
            nodelist = vcat(p.subpaths[1].arcs[1][1], [a[2] for a in p.subpaths[1].arcs])
            println(nodelist)
            if length(nodelist) == 2 && p.subpaths[1].current_time == 0
                push!(results_paths_withcharge, (val, p))
            else
                (p_feasible, pure_path_label) = get_postcharge_shortest_pure_path_label(
                    data, graph, nodelist,
                    ;
                    time_windows = time_windows,
                )
                println(p_feasible)
                feasible = feasible && p_feasible
                if !p_feasible
                    break
                end
                p_ = convert_pure_path_label_to_path(pure_path_label, graph)
                push!(results_paths_withcharge, (val, p_))
            end
        end
        if feasible
            break
        end
        time_heuristic_slack = time_heuristic_slack - 0.05
        if time_heuristic_slack ≤ 0.0
            error()
        end
    end
    return results_paths_withcharge
end