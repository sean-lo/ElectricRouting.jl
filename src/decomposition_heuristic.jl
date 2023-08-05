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

    modified_costs = compute_arc_modified_costs(data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                BaseSubpathLabel,
            }()
            for current_node in union(data.N_depots, data.N_customers)
        )
        for starting_node in data.N_depots
    )
    unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
    for node in data.N_depots
        served = zeros(Int, data.n_customers)
        if check_customers
            key = (0, served...)
        else
            key = (0,)
        end
        base_labels[node][node][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, data.n_customers),
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
        for next_node in setdiff(outneighbors(data.G, current_node), current_node, data.N_charging)
            if single_service
                if next_node in data.N_customers && current_subpath.served[next_node] > 0
                    continue
                end
            end
            # Preventing customer 2-cycles (Christofides)
            if christofides && next_node in data.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            # time feasibility
            new_time_taken = current_subpath.time_taken + data.t[current_node, next_node]
            if time_windows
                new_time_taken = max(
                    new_time_taken,
                    data.α[next_node]
                )
                if new_time_taken > data.β[next_node]
                    continue
                end
            else
                if new_time_taken + data.min_t[next_node] > T_heuristic
                    continue
                end
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken = new_time_taken
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in data.N_customers
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
            if added && next_node in data.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for depot in data.N_depots
        for path in values(base_labels[depot][depot])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end

function find_nondominated_paths_nocharge_ngroute(
    data::EVRPData,
    neighborhoods::NGRouteNeighborhood,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_windows::Bool = false,
    christofides::Bool = false,
)

    modified_costs = compute_arc_modified_costs(data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => Dict{
                Tuple{Vararg{Int}}, 
                SortedDict{
                    Tuple{Vararg{Int}}, 
                    BaseSubpathLabel,
                },
            }()
            for current_node in union(data.N_depots, data.N_customers)
        )
        for starting_node in data.N_depots
    )
    unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
    for node in data.N_depots
        base_labels[node][node][(node,)] = SortedDict{
            Tuple{Vararg{Int}}, 
            BaseSubpathLabel,
        }(
            Base.Order.ForwardOrdering(),
            (0,) => BaseSubpathLabel(
                0, 0, 0.0, [node,], zeros(Int, data.n_customers),
            )
        )
        push!(unexplored_states, (0, node, node))
    end

    while length(unexplored_states) > 0
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        current_node = state[end]
        current_key = state[1:end-2]
        for current_set in keys(base_labels[starting_node][current_node])
            if !(current_key in keys(base_labels[starting_node][current_node][current_set]))
                continue
            end
            current_subpath = base_labels[starting_node][current_node][current_set][current_key]
            for next_node in setdiff(outneighbors(data.G, current_node), current_node, data.N_charging)
                if next_node in current_set
                    # if next_node is a customer not yet visited, proceed
                    # only if one can extend current_subpath along next_node according to ng-route rules
                    continue
                end
                # Preventing customer 2-cycles (Christofides)
                if christofides && next_node in data.N_customers
                    if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                        continue
                    end
                end
                # time feasibility
                new_time_taken = current_subpath.time_taken + data.t[current_node, next_node]
                if time_windows
                    new_time_taken = max(
                        new_time_taken,
                        data.α[next_node]
                    )
                    if new_time_taken > data.β[next_node]
                        continue
                    end
                else
                    if new_time_taken + data.min_t[next_node] > T_heuristic
                        continue
                    end
                end

                new_subpath = copy(current_subpath)
                new_subpath.time_taken = new_time_taken
                new_subpath.cost += modified_costs[current_node, next_node]
                push!(new_subpath.nodes, next_node)
                if next_node in data.N_customers
                    new_subpath.served[next_node] += 1
                end

                new_set = ngroute_create_set(neighborhoods, current_set, next_node)
                if !(new_set in keys(base_labels[starting_node][next_node]))
                    base_labels[starting_node][next_node][new_set] = SortedDict{
                        Tuple{Vararg{Int}}, 
                        BaseSubpathLabel,
                    }()
                end
                new_key = (new_subpath.time_taken,)
                added = add_subpath_longlabel_to_collection!(
                    base_labels[starting_node][next_node][new_set],
                    new_key, new_subpath,
                    ;
                    verbose = false,
                )
                if added && next_node in data.N_customers
                    new_state = (new_key..., starting_node, next_node)
                    push!(unexplored_states, new_state)
                end
            end
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for depot in data.N_depots
        for set in keys(base_labels[depot][depot])
            for path in values(base_labels[depot][depot][set])
                if length(path.nodes) == 1
                    path.nodes = [depot, depot]
                end
            end
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots
            for set in keys(base_labels[starting_node][end_node])
                for v in values(base_labels[starting_node][end_node][set])
                    v.cost = v.cost - κ[starting_node]
                end
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots
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
    neighborhoods::NGRouteNeighborhood,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    time_windows::Bool = false,
    christofides::Bool = false,
)

    modified_costs = compute_arc_modified_costs(data, ν)

    base_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                BaseSubpathLabel,
            }()
            for current_node in union(data.N_depots, data.N_customers)
        )
        for starting_node in data.N_depots
    )
    unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
    for node in data.N_depots
        node_labels = zeros(Int, data.n_customers + data.n_depots)
        node_labels[node] = 1
        key = (0, node_labels...)
        base_labels[node][node][key] = BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, data.n_customers),
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
        current_set = state[2:end-2]
        current_subpath = base_labels[starting_node][current_node][current_key]
        for next_node in setdiff(outneighbors(data.G, current_node), current_node, data.N_charging)
            if next_node in data.N_customers && current_set[next_node] == 1
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                continue
            end
            # Preventing customer 2-cycles (Christofides)
            if christofides && next_node in data.N_customers
                if length(current_subpath.nodes) ≥ 2 && current_subpath.nodes[end-1] == next_node
                    continue
                end
            end
            # time feasibility
            new_time_taken = current_subpath.time_taken + data.t[current_node, next_node]
            if time_windows
                new_time_taken = max(
                    new_time_taken,
                    data.α[next_node]
                )
                if new_time_taken > data.β[next_node]
                    continue
                end
            else
                if new_time_taken + data.min_t[next_node] > T_heuristic
                    continue
                end
            end

            new_subpath = copy(current_subpath)
            new_subpath.time_taken = new_time_taken
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in data.N_customers
                new_subpath.served[next_node] += 1
            end

            new_set = ngroute_create_set_alt(neighborhoods, collect(current_set), next_node)
            new_key = (new_subpath.time_taken, new_set...)
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node],
                new_key, new_subpath,
                ;
                verbose = false,
            )
            if added && next_node in data.N_customers
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_customers
            delete!(base_labels[starting_node], end_node)
        end
    end

    for depot in data.N_depots
        for path in values(base_labels[depot][depot])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
            end
        end
    end

    for starting_node in data.N_depots
        for end_node in data.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
    for end_node in data.N_depots
        for starting_node in data.N_depots
            for v in values(base_labels[starting_node][end_node])
                v.cost = v.cost - μ[end_node]
            end
        end
    end

    return base_labels

end


function get_negative_base_labels_from_base_labels(
    data::EVRPData, 
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
        for starting_node in data.N_depots
            for end_node in data.N_depots
                for (key, base_label) in base_labels[starting_node][end_node]
                    if base_label.cost < -1e-6
    ]
end

function get_negative_base_labels_from_base_labels_ngroute(
    data::EVRPData, 
    base_labels::Dict{
        Int, 
        Dict{
            Int, 
            Dict{
                Tuple{Vararg{Int}},
                SortedDict{
                    Tuple{Vararg{Int}},
                    BaseSubpathLabel,
                },
            },
        },
    },
)
    return BaseSubpathLabel[
        base_label
        for starting_node in data.N_depots
            for end_node in data.N_depots
                for set in keys(base_labels[starting_node][end_node])
                    for (key, base_label) in base_labels[starting_node][end_node][set]
                        if base_label.cost < -1e-6
    ]
end

function subproblem_iteration_nocharge(
    data::EVRPData,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    T_heuristic::Int,
    ;
    neighborhoods::Union{Nothing, NGRouteNeighborhood} = nothing,
    time_windows::Bool = false,
    single_service::Bool = false,
    check_customers::Bool = false,
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
)
    if ngroute && !ngroute_alt
        base_labels_result = @timed find_nondominated_paths_nocharge_ngroute(
            data, neighborhoods, κ, μ, ν, T_heuristic,
            ;
            time_windows = time_windows,
            christofides = christofides,
        )
    elseif ngroute && ngroute_alt
        base_labels_result = @timed find_nondominated_paths_nocharge_ngroute_alt(
            data, neighborhoods, κ, μ, ν, T_heuristic,
            ;
            time_windows = time_windows,
            christofides = christofides,
        )
    else
        base_labels_result = @timed find_nondominated_paths_nocharge(
            data, κ, μ, ν, T_heuristic,
            ;
            time_windows = time_windows,
            christofides = christofides,
            single_service = single_service,
            check_customers = check_customers,
        )
    end
    if ngroute && !ngroute_alt
        negative_base_labels = get_negative_base_labels_from_base_labels_ngroute(data, base_labels_result.value)
    elseif (ngroute && ngroute_alt) || !ngroute
        negative_base_labels = get_negative_base_labels_from_base_labels(data, base_labels_result.value)
    end
    return (negative_base_labels, length(negative_base_labels), base_labels_result.time)
end

function get_paths_from_negative_base_labels(
    data::EVRPData,
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
                    n_customers = data.n_customers,
                    starting_node = base_label.nodes[1],
                    starting_time = 0,
                    starting_charge = data.B,
                    current_node = base_label.nodes[end],
                    arcs = collect(zip(base_label.nodes[1:end-1], base_label.nodes[2:end])),
                    current_time = base_label.time_taken,
                    current_charge = data.B,
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
    ;
    time_heuristic_slack::Float64 = 0.9,
    Env = nothing,
    time_windows::Bool = false,
    path_single_service::Bool = false,
    path_check_customers::Bool = false,
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(data.n_customers))),
    verbose::Bool = true,
    time_limit::Float64 = Inf,
    max_iters::Float64 = Inf,
)

    start_time = time()

    if ngroute
        neighborhoods = compute_ngroute_neighborhoods(
            data,
            ngroute_neighborhood_size;
        )
    else
        neighborhoods = nothing
    end
    T_heuristic = Int(round(time_heuristic_slack * (data.T + data.B) * data.μ / (1 + data.μ)))

    some_paths = generate_artificial_paths(data)
    path_costs = compute_path_costs(
        data, 
        some_paths,
    )
    path_service = compute_path_service(
        data, 
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
            data.n_customers,
            data.n_depots,
            data.n_charging,
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
        data, some_paths, path_costs, path_service,
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
            "κ" => Dict(zip(data.N_depots, dual.(model[:κ]).data)),
            "μ" => Dict(zip(data.N_depots, dual.(model[:μ]).data)),
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
            data, CGLP_results["κ"], CGLP_results["μ"], CGLP_results["ν"], T_heuristic,
            ;
            neigborhoods = neighborhoods,
            time_windows = time_windows,
            single_service = path_single_service,
            check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
        )
        (generated_paths) = get_paths_from_negative_base_labels(
            data, negative_base_labels,
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
            data,
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
    nodelist::Vector{Int},
    ;
    time_windows::Bool = false,
)
    modified_costs = compute_arc_modified_costs(data, zeros(Float64, data.n_customers))
    t = data.t
    B = data.B
    q = data.q
    if time_windows
        α = data.α
        β = data.β
    else
        α = zeros(Int, data.n_nodes)
        β = repeat([data.T], data.n_nodes)
    end

    pure_path_labels = Dict(
        (nodelist[1:j]...,) => Dict(
            current_node => SortedDict{
                Tuple{Vararg{Int}}, 
                PurePathLabel,
            }()
            for current_node in data.N_nodes
        )
        for j in eachindex(nodelist)
    )
    key = (0, -data.B, data.B)
    starting_node = nodelist[1]
    nodeseq = (starting_node,)
    pure_path_labels[nodeseq][starting_node][key] = PurePathLabel(
        0.0, 
        [starting_node],
        Int[],
        Int[],
        0,
        0,
        data.B,
        data.B,
        false,
        zeros(Int, data.n_customers),
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
        for next_node in setdiff(vcat(eventual_next_node, data.N_charging), current_node)
            if !(next_node in outneighbors(data.G, current_node))
                continue
            end
            
            (feasible, new_path) = compute_new_pure_path(
                current_path,
                current_node,
                next_node, 
                data,
                α, β, modified_costs,
            )
            if !feasible
                continue 
            end
            # add new_path to collection
            if next_node in union(data.N_customers, data.N_depots)
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
            if added && !(next_node in data.N_depots)
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
    ;
    Env = nothing,
    time_windows::Bool = false,
    path_single_service::Bool = false,
    path_check_customers::Bool = false,
    christofides::Bool = false,
    ngroute::Bool = false,
    ngroute_alt::Bool = false,
    ngroute_neighborhood_size::Int = Int(ceil(sqrt(data.n_customers))),
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
            data,
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
            results_paths = collect_path_solution_support(CGIP_results, some_paths, data)
        else
            results_paths = collect_path_solution_support(CGLP_results, some_paths, data)
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
                    data, nodelist,
                    ;
                    time_windows = time_windows,
                )
                println(p_feasible)
                feasible = feasible && p_feasible
                if !p_feasible
                    break
                end
                p_ = convert_pure_path_label_to_path(pure_path_label, data)
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