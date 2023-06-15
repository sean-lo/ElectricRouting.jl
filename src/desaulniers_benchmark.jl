include("utils.jl")
using DataStructures
using Printf

Base.@kwdef mutable struct PurePathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    explored::Bool
    served::Vector{Int}
    artificial::Bool = false
end

Base.copy(p::PurePathLabel) = PurePathLabel(
    p.cost,
    copy(p.nodes),
    copy(p.excesses),
    copy(p.slacks),
    p.time_mincharge,
    p.time_maxcharge,
    p.charge_mincharge,
    p.charge_maxcharge,
    p.explored,
    copy(p.served),
    p.artificial,
)

function compute_path_label_cost(
    p::PurePathLabel,
    data,
    M::Float64 = 1e10,
    ;
    verbose = false,
)
    if p.artificial
        return M
    elseif length(p.nodes) ≤ 1
        return 0
    end

    arcs = collect(zip(p.nodes[1:end-1], p.nodes[2:end]))
    cost = data["travel_cost_coeff"] * sum(data["c"][a...] for a in arcs)
    verbose && @printf("Path cost: \t\t%11.3f\n", cost)

    charging_cost = data["charge_cost_coeff"] * (sum(p.slacks) + sum(p.excesses))
    verbose && @printf("Charging cost: \t\t%11.3f\n", charging_cost)
    cost += charging_cost

    return cost
end


function compute_path_label_modified_cost(
    p::PurePathLabel,
    data,
    κ,
    μ,
    ν,
    ;
    verbose = false,
)
    reduced_cost = compute_path_label_cost(p, data, verbose = verbose)

    service_cost = 0.0
    for (j, c) in enumerate(p.served)
        service_cost += (c * -ν[j])
    end
    verbose && @printf("Service cost: \t\t%11.3f\n", service_cost)
    reduced_cost += service_cost

    verbose && @printf("Starting depot cost: \t%11.3f\n", (- κ[p.nodes[1]]))
    reduced_cost = reduced_cost - κ[p.nodes[1]]

    verbose && @printf("Ending depot cost: \t%11.3f\n", (- μ[p.nodes[end]]))
    reduced_cost = reduced_cost - μ[p.nodes[end]]

    verbose && @printf("Total modified cost: \t%11.3f\n\n", reduced_cost)

    return reduced_cost
end

function find_nondominated_paths(
    G,
    data, 
    κ,
    μ,
    ν,
    ;
    single_service::Bool = false,
    check_customers::Bool = true, 
    time_windows::Bool = true,
)
    if single_service
        S = data["N_customers"]
    else
        S = Int[]
    end
    return find_nondominated_paths_S(
        G, data, κ, μ, ν,
        ;
        S = S,     
        check_customers = check_customers, 
        time_windows = time_windows,
    )
end


function find_nondominated_paths_S(
    G,
    data, 
    κ,
    μ,
    ν,
    ;
    S::Vector{Int} = [],
    check_customers::Bool = true, # indicates if labels for customers are used in label domination
    time_windows::Bool = true,
    initial_pure_path_labels::Union{
        Nothing,
        Dict{Int, Dict{Int, SortedDict{
            Tuple{Vararg{Int}},
            PurePathLabel,
            Base.Order.ForwardOrdering
        }}},
    } = nothing,
)
    # finds non-dominated paths where customers in S appear at most once.
    function add_pure_path_label_to_collection!(
        collection::SortedDict{
            Tuple{Vararg{Int}}, 
            PurePathLabel,
            Base.Order.ForwardOrdering
        },
        key::Tuple{Vararg{Int}},
        path::PurePathLabel,
        ;
        verbose::Bool = false,
    )
        added = true
        for (k, p) in pairs(collection)
            if p.cost ≤ path.cost
                if all(k .≤ key)
                    added = false
                    if verbose
                        println("$(key), $(path.cost) dominated by $(k), $(p.cost)")
                    end
                    break
                end
            end
            if path.cost ≤ p.cost
                if all(key .≤ k)
                    if verbose
                        println("$(key), $(path.cost) dominates $(k), $(p.cost)")
                    end
                    pop!(collection, k)
                end
            end
        end
        if added
            if verbose
                println("$(key), $(path.cost) added!")
            end
            insert!(collection, key, path)
        end
        return added
    end

    modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
    for j in data["N_customers"]
        for i in data["N_nodes"]
            modified_costs[i,j] -= ν[j]
        end
    end

    if isnothing(initial_pure_path_labels)
        pure_path_labels = Dict(
            starting_node => Dict(
                current_node => SortedDict{
                    Tuple{Vararg{Int}}, 
                    PurePathLabel,
                }()
                for current_node in data["N_nodes"]
            )
            for starting_node in data["N_depots"]
        )
        if check_customers
            # label key here has the following fields:
            # 1) current minimum time T_i(min)
            # 2) negative of current max charge -B_i(max)
            # 3) difference between min time and min charge, T_i(min) - B_i(min)
            # 4) if applicable, whether i-th customer served
            key = (0, -data["B"], -data["B"], zeros(Int, length(S))...)
        else
            key = (0, -data["B"], -data["B"])
        end
        for depot in data["N_depots"]
            pure_path_labels[depot][depot][key] = PurePathLabel(
                0.0,
                [depot],
                Int[],
                Int[],
                0,
                0,
                data["B"],
                data["B"],
                false,
                zeros(Int, data["n_customers"]),
                false,
            )
        end
    else
        pure_path_labels = deepcopy(initial_pure_path_labels)
        # reset pure_path_labels
        for starting_node in data["N_depots"]
            for ending_node in keys(pure_path_labels[starting_node])
                for path in values(pure_path_labels[starting_node][ending_node])
                    path.explored = false
                    path.cost = path.cost + κ[starting_node]
                    if ending_node in data["N_depots"]
                        path.cost = path.cost + μ[ending_node]
                    end
                end
            end
        end
    end

    unexplored_states = SortedSet(
        [
            (key..., starting_node, current_node)
            for starting_node in data["N_depots"]
                for current_node in keys(pure_path_labels[starting_node])
                    for key in keys(pure_path_labels[starting_node][current_node])
        ]
    )

    t = data["t"]
    B = data["B"]
    q = data["q"]
    if time_windows
        α = data["α"]
        β = data["β"]
    else
        α = zeros(Int, data["n_nodes"])
        β = repeat([data["T"]], data["n_nodes"])
    end

    while length(unexplored_states) > 0
        # println(length(unexplored_states))
        state = pop!(unexplored_states)
        starting_node = state[end-1]
        i = state[end]
        if !(state[1:end-2] in keys(pure_path_labels[starting_node][i]))
            continue
        end
        path = pure_path_labels[starting_node][i][state[1:end-2]]
        # println("$(path.time_mincharge), $(path.time_maxcharge), $(path.charge_mincharge), $(path.charge_maxcharge)")
        for j in setdiff(outneighbors(G, i), i)
            # println("$state -> $j")
            if j in data["N_customers"]
                if j in S && path.served[j] > 0
                    # println("already served $j")
                    continue
                end
                if length(path.nodes) ≥ 2 && path.nodes[end-1] == j
                    continue
                end
            end
            # feasibility checks
            # (1) battery
            excess = max(
                0, 
                q[i,j] - path.charge_mincharge 
            )
            # (2) time windows
            if path.time_mincharge + excess + t[i,j] > β[j]
                # println("$(path.time_mincharge), $excess, $(t[i,j]), $(β[j])")
                # println("not time windows feasible")
                continue
            end
            if path.time_mincharge + excess + t[i,j] + data["min_t"][j] > data["T"]
                continue
            end
            # (3) charge interval 
            if (
                (i in data["N_charging"] && excess > max(B - path.charge_mincharge, 0))
                || 
                (!(i in data["N_charging"]) && excess > max(path.charge_maxcharge - path.charge_mincharge, 0))
            )
                # if i in data["N_charging"]
                #     println("$excess, $(B), $(path.charge_mincharge)")
                # else
                #     println("$excess, $(path.charge_maxcharge), $(path.charge_mincharge)")
                # end
                # println("not charge feasible")
                continue
            end
            
            new_path = copy(path)
            push!(new_path.nodes, j)
            if j in data["N_customers"]
                new_path.served[j] += 1
            end

            push!(new_path.excesses, excess)
            new_path.time_mincharge = max(
                α[j],
                path.time_mincharge + t[i,j] + excess
            )
            if i in data["N_charging"]
                slack = max(
                    # floating point accuracy
                    0, 
                    min(
                        new_path.time_mincharge - (path.time_mincharge + t[i,j] + excess),
                        B - (path.charge_mincharge + excess),
                    )
                )
                push!(new_path.slacks, slack)
                new_path.time_maxcharge = min(
                    β[j],
                    max(
                        α[j],
                        path.time_mincharge + (B - path.charge_mincharge) + t[i,j],
                    )
                )
            else
                slack = max(
                    # floating point accuracy
                    0, 
                    min(
                        new_path.time_mincharge - (path.time_mincharge + t[i,j] + excess),
                        path.charge_maxcharge - (path.charge_mincharge + excess),
                    )
                )
                push!(new_path.slacks, slack)
                new_path.time_maxcharge = min(
                    β[j],
                    max(
                        α[j],
                        path.time_maxcharge + t[i,j],
                    )
                )
            end
            
            new_path.charge_mincharge = (
                path.charge_mincharge 
                + excess 
                + slack
                - q[i,j]
            )
            new_path.charge_maxcharge = (
                new_path.charge_mincharge 
                + new_path.time_maxcharge 
                - new_path.time_mincharge
            )

            new_path.cost += modified_costs[i,j]
            new_path.cost += data["charge_cost_coeff"] * (slack + excess)
            
            # println("$(new_path.time_mincharge), $(new_path.time_maxcharge), $(new_path.charge_mincharge), $(new_path.charge_maxcharge)")
            # add new_path to collection
            if check_customers
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge, 
                    new_path.served[S]...,
                )
            else
                new_key = (
                    new_path.time_mincharge, 
                    - new_path.charge_maxcharge, 
                    new_path.time_mincharge - new_path.charge_mincharge
                )
            end

            added = add_pure_path_label_to_collection!(pure_path_labels[starting_node][j], new_key, new_path, verbose = false)
            if added && !(j in data["N_depots"])
                new_state = (new_key..., starting_node, j)
                # println("adding state: $(new_state)")
                push!(unexplored_states, new_state)
            end
        end
    end
    
    for starting_node in data["N_depots"]
        for ending_node in keys(pure_path_labels[starting_node])
            for path in values(pure_path_labels[starting_node][ending_node])
                path.cost = path.cost - κ[starting_node]
            end
        end
    end

    for starting_node in keys(pure_path_labels)
        for ending_node in intersect(data["N_depots"], keys(pure_path_labels[starting_node]))
            for path in values(pure_path_labels[starting_node][ending_node])
                path.cost = path.cost - μ[ending_node]
            end
        end
    end

    return pure_path_labels
end

function get_negative_pure_path_labels_from_pure_path_labels(
    data, 
    pure_path_labels::Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PurePathLabel,
        Base.Order.ForwardOrdering
    }}}
)
    return Dict(
        starting_node => Dict(
            end_node => SortedDict{
                Tuple{Vararg{Int}},
                PurePathLabel,
            }(
                key => path_label
                for (key, path_label) in pure_path_labels[starting_node][end_node] 
                    if path_label.cost < -1e-6
            )
            for end_node in keys(pure_path_labels[starting_node])
        )
        for starting_node in data["N_depots"]
    )
end

function subproblem_iteration_benchmark(
    G, data, κ, μ, ν,
    ;
    time_windows::Bool = false,
    path_single_service::Bool = true,
    path_check_customers::Bool = true,
)
    pure_path_labels_result = @timed find_nondominated_paths(
        G, data, κ, μ, ν,
        ;
        time_windows = time_windows, 
        single_service = path_single_service, 
        check_customers = path_check_customers,
    )
    pure_path_labels_time = pure_path_labels_result.time
    negative_pure_path_labels = get_negative_pure_path_labels_from_pure_path_labels(data, pure_path_labels_result.value)
    negative_pure_path_labels_count = sum(
        length(negative_pure_path_labels[starting_node][end_node]) 
        for starting_node in data["N_depots"]
            for end_node in keys(negative_pure_path_labels[starting_node]) 
    )
    return (negative_pure_path_labels, negative_pure_path_labels_count, pure_path_labels_time)
end

function subproblem_iteration_benchmark_incremental_elementarity(
    G, data, κ, μ, ν,
    ;
    time_windows::Bool = false,
    path_check_customers::Bool = true,
    verbose::Bool = false,
    warm_start::Bool = false,
)
    start_time = time()
    S = Int[]
    initial_pure_path_labels = nothing
    while true
        pure_path_labels = find_nondominated_paths_S(
            G, data, κ, μ, ν,
            ;
            S = S,
            time_windows = time_windows, 
            check_customers = path_check_customers,
            initial_pure_path_labels = initial_pure_path_labels
        )
        # filter for path labels: (i) ending at depot, 
        # (ii) w/ negative reduced cost,
        # (iii) elementary
        d_pure_path_labels = get_depot_pure_path_labels(data, pure_path_labels);
        n_d_pure_path_labels = get_negative_pure_path_labels_from_pure_path_labels(data, d_pure_path_labels);
        n_d_pure_path_labels_count = sum(
            length(n_d_pure_path_labels[starting_node][end_node]) 
            for starting_node in data["N_depots"]
                for end_node in keys(n_d_pure_path_labels[starting_node]) 
        )
        if n_d_pure_path_labels_count == 0
            # this activates if the overall CG converges
            return (n_d_pure_path_labels, 0, time() - start_time)
        end
        (e_n_d_pure_path_labels, ne_n_d_pure_path_labels) = get_elementary_nonelementary_pure_path_labels(
            data, n_d_pure_path_labels; 
            S = data["N_customers"]
        );
        e_n_d_pure_path_labels_count = sum(
            length(e_n_d_pure_path_labels[starting_node][end_node]) 
            for starting_node in data["N_depots"]
                for end_node in keys(e_n_d_pure_path_labels[starting_node]) 
        )
        if e_n_d_pure_path_labels_count > 0
            # this activates if there is some negative path found (elementary w.r.t. all customers)
            return (e_n_d_pure_path_labels, e_n_d_pure_path_labels_count, time() - start_time)
        end

        # else, expand S
        # println(length(S))
        if length(S) == data["n_customers"]
            error()
        end
        # select a path which is nondominated of minimal cost 
        (min_cost, _, _, _, path_label) = find_minimum_reduced_cost_pure_path_label(data, d_pure_path_labels)
        verbose && println("$min_cost, $path_label")
        if min_cost >= -1e-6
            error()
        end
        (val, j) = findmax(path_label.served)
        # println("$val, $j")
        push!(S, j)
        sort!(S)
        verbose && println(S)
        
        # prepare warm start
        if warm_start
            initial_pure_path_labels, _ = get_elementary_nonelementary_pure_path_labels(
                data, pure_path_labels, 
                ; 
                S = S, rename_keys = true,
            )
        end
    end
end


function get_depot_pure_path_labels(
    data, 
    pure_path_labels::Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PurePathLabel,
        Base.Order.ForwardOrdering
    }}},
)
    return Dict(
        starting_node => Dict(
            end_node => pure_path_labels[starting_node][end_node]
            for end_node in data["N_depots"]
        )
        for starting_node in data["N_depots"]
    )
end

function get_elementary_nonelementary_pure_path_labels(
    data, 
    pure_path_labels::Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PurePathLabel,
        Base.Order.ForwardOrdering
    }}},
    ;
    S::Vector{Int} = data["N_customers"],
    rename_keys::Bool = false,
)
    # divides pure path labels in pure_path_labels into two:
    # those that are elementary (with respect to S) 
    # and those that visit a customer in S more than once.
    elementary_pure_path_labels = Dict(
        starting_node => Dict(
            end_node => SortedDict{
                Tuple{Vararg{Int}}, 
                PurePathLabel,
            }()
            for end_node in keys(pure_path_labels[starting_node])
        )
        for starting_node in data["N_depots"]
    )
    nonelementary_pure_path_labels = Dict(
        starting_node => Dict(
            end_node => SortedDict{
                Tuple{Vararg{Int}}, 
                PurePathLabel,
            }()
            for end_node in keys(pure_path_labels[starting_node])
        )
        for starting_node in data["N_depots"]
    )
    for starting_node in data["N_depots"]
        for end_node in keys(pure_path_labels[starting_node])
            for (key, path_label) in pairs(pure_path_labels[starting_node][end_node])
                if rename_keys
                    new_key = (key[1:3]..., path_label.served[S]...)
                else
                    new_key = key 
                end
                if all(path_label.served[S] .≤ 1)
                    push!(elementary_pure_path_labels[starting_node][end_node], new_key => path_label)
                else
                    push!(nonelementary_pure_path_labels[starting_node][end_node], new_key => path_label)
                end
            end
        end
    end
    return (elementary_pure_path_labels, nonelementary_pure_path_labels)
end

function find_minimum_reduced_cost_pure_path_label(
    data,
    pure_path_labels::Dict{Int, Dict{Int, SortedDict{
        Tuple{Vararg{Int}},
        PurePathLabel,
        Base.Order.ForwardOrdering
    }}},
)
    min_cost = Inf
    min_starting_node = nothing
    min_end_node = nothing
    min_key = nothing 
    min_path_label = nothing
    for starting_node in data["N_depots"]
        for end_node in keys(pure_path_labels[starting_node])
            for (key, path_label) in pairs(pure_path_labels[starting_node][end_node])
                if path_label.cost < min_cost
                    min_cost = path_label.cost
                    min_starting_node = starting_node
                    min_end_node = end_node
                    min_key = key
                    min_path_label = path_label
                end
            end
        end
    end
    return (min_cost, min_starting_node, min_end_node, min_key, min_path_label)
end