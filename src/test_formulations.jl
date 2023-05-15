include("arc_formulation.jl")
include("subpath_formulation.jl")
include("utils.jl")

using Distributions

data = generate_instance(
    ;
    n_depots = 2,
    n_customers = 10,
    n_charging = 2,
    n_vehicles = 4,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.6,
    T = 1200.0,
    seed = 0,
    B = 800.0,
    μ = 5.0,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.2,
)
plot_instance(data)

arc_results, arc_params = arc_formulation(data, with_charging = true, time_limit = 60)
arc_paths = construct_paths_from_arc_solution(arc_results, data)
arc_results_printout(
    arc_results, 
    arc_params,
    data,
    with_charging = true,
)

G = construct_graph(data)
CGLP_results, CGIP_results, params, printlist, some_subpaths, some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data);
print.(printlist);

subpath_results_printout(
    CGLP_results,
    params,
    data,
    some_subpaths,
    some_charging_arcs,
)

subpath_results_printout(
    CGIP_results,
    params,
    data,
    some_subpaths,
    some_charging_arcs,
)

G_sparse = construct_sparse_graph(data, 1.0)
(
    CGLP_results_sparse, 
    CGIP_results_sparse, 
    params_sparse, 
    printlist_sparse, 
    some_subpaths_sparse, 
    some_charging_arcs_sparse,
) = subpath_formulation_column_generation_integrated_from_paths(G_sparse, data);

subpath_results_printout(
    CGLP_results,
    params,
    data,
    some_subpaths,
    some_charging_arcs,
)

subpath_results_printout(
    CGIP_results,
    params,
    data,
    some_subpaths,
    some_charging_arcs,
)



data_large = generate_instance(
    ;
    n_depots = 2,
    n_customers = 25,
    n_charging = 3,
    n_vehicles = 3,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.8,
    T = 2000.0,
    seed = 6,
    B = 600.0,
    μ = 5.0,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
)
plot_instance(data_large)
G_large = construct_graph(data_large)
(
    CGLP_results_large, CGIP_results_large, params_large, printlist_large, some_subpaths_large, some_charging_arcs_large 
) = subpath_formulation_column_generation_integrated_from_paths(G_large, data_large);
print.(printlist_large);
subpath_results_printout(
    CGLP_results_large,
    params_large,
    data_large,
    some_subpaths_large,
    some_charging_arcs_large,
)
### Scratch work

some_subpaths

CGLP_results, CGIP_results, params, printlist, some_subpaths, some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data);

κ = params["κ"][6]
μ = params["μ"][6]
ν = params["ν"][6]

modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
for j in data["N_customers"]
    for i in data["N_nodes"]
        modified_costs[i,j] -= ν[j]
    end
end
modified_costs

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
        if v2.cost ≤ v1.cost
            if k2 ≤ k1
                added = false
                break
            end
        elseif v1.cost ≤ v2.cost
            if k1 ≤ k2
                pop!(collection, k2)
            end
        end
    end
    if added
        insert!(collection, k1, v1)
    end
    return added
end

for new_node in 1:5
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

function merge_collections(
    labels1::SortedDict{Float64, SubpathWithCost},
    labels2::SortedDict{Float64, SubpathWithCost},
)
    keys1 = collect(keys(labels1))
    keys2 = collect(keys(labels2))

    new = []
    for (t, cost, i, j) in sort([
        (k1 + k2, s1.cost + s2.cost, i, j)
        for (i, (k1, s1)) in enumerate(pairs(labels1)),
            (j, (k2, s2)) in enumerate(pairs(labels2))
            if k1 + k2 ≤ data["T"]
    ])
        if length(new) == 0
            push!(new, (t, cost, i, j))
        elseif new[end][1] < t
            if new[end][2] > cost
                push!(new, (t, cost, i, j))
            end
        else
            if new[end][2] > cost
                new = new[1:end-1]
                push!(new, (t, cost, i, j))
            end
        end
    end
    new_labels = SortedDict(
        t => SubpathWithCost(
            cost = cost,
            n_customers = data["n_customers"],
            starting_node = labels1[keys1[i]].starting_node,
            starting_time = 0.0,
            starting_charge = data["B"],
            current_node = labels2[keys2[j]].current_node,
            current_time = t,
            current_charge = data["B"] - t,
            arcs = vcat(labels1[keys1[i]].arcs, labels2[keys2[j]].arcs),
            served = labels1[keys1[i]].served .|| labels2[keys2[j]].served,
        )
        for (t, cost, i, j) in new
    )
    return new_labels
end

function direct_sum_of_collections(
    labels1::SortedDict{Float64, SubpathWithCost},
    labels2::SortedDict{Float64, SubpathWithCost},
)
    keys1 = collect(keys(labels1))
    keys2 = collect(keys(labels2))

    new = []
    for (t, cost, i, j) in sort([
        (k1 + k2, s1.cost + s2.cost, i, j)
        for (i, (k1, s1)) in enumerate(pairs(labels1)),
            (j, (k2, s2)) in enumerate(pairs(labels2))
            if k1 + k2 ≤ data["T"]
    ])
        if length(new) == 0
            push!(new, (t, cost, i, j))
        elseif new[end][1] < t
            if new[end][2] > cost
                push!(new, (t, cost, i, j))
            end
        else
            if new[end][2] > cost
                new = new[1:end-1]
                push!(new, (t, cost, i, j))
            end
        end
    end
    new_labels = SortedDict(
        t => SubpathWithCost(
            cost = cost,
            n_customers = data["n_customers"],
            starting_node = labels1[keys1[i]].starting_node,
            starting_time = 0.0,
            starting_charge = data["B"],
            current_node = labels2[keys2[j]].current_node,
            current_time = t,
            current_charge = data["B"] - t,
            arcs = vcat(labels1[keys1[i]].arcs, labels2[keys2[j]].arcs),
            served = labels1[keys1[i]].served .|| labels2[keys2[j]].served,
        )
        for (t, cost, i, j) in new
    )
    return new_labels
end

other_labels = direct_sum_of_collections(base_labels[1][6], base_labels[6][2])
old_labels = base_labels[1][2]
keys1 = collect(keys(old_labels))
left = popfirst!(keys1)
keys2 = collect(keys(other_labels))
right = popfirst!(keys2)

new_pairs = []
new_labels = []
while length(keys1) > 0 || length(keys2) > 0
    if left < right
        time = left 
        cost = old_labels[left].cost
        println("left: $left")
        if length(keys1) > 0
            left = popfirst!(keys1)
        else
            left = Inf
        end
    else
        time = right 
        cost = other_labels[right].cost
        println("right: $right")
        if length(keys2) > 0
            right = popfirst!(keys2)
        else
            right = Inf
        end
    end
    if length(new_pairs) == 0
        push!(new_pairs, (time, cost))
    else
        if new_pairs[end][2] > cost 
            if new_pairs[end][1] < time
                push!(new_pairs, (time, cost))
            elseif new_pairs[end][1] == time_limit_sec
                new_pairs = vcat(new_pairs[1:end-1], [(time, cost)])
            end
        end
    end
end

[
    (t, s.cost) 
    for (t, s) in pairs(new_labels)
]

some_subpaths = generate_artificial_subpaths(data)
subpath_costs = compute_subpath_costs(
    data, 
    some_subpaths,
)
subpath_service = compute_subpath_service(
    data, 
    some_subpaths,
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

flow_conservation_exprs_out = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
flow_conservation_exprs_in = Dict{Tuple{Int, Float64, Float64}, AffExpr}()
flow_conservation_constrs = Dict{Tuple{Int, Float64, Float64}, ConstraintRef}()

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
@objective(mp_model, Min, subpath_costs_expr)



optimize!(mp_model)
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
)




base_labels = @time generate_base_subpaths(
    G, 
    data,
    mp_results["κ"],
    mp_results["μ"],
    mp_results["ν"],
)

full_labels = @time find_nondominated_paths(
    data, base_labels,
    mp_results["κ"],
    mp_results["μ"],
)

(generated_subpaths, generated_charging_arcs) = @time get_subpaths_charging_arcs_from_negative_paths(
    data,
    full_labels,
)
sum(length(x) for x in generated_subpaths)
sum(length(x) for x in generated_charging_arcs)

some_subpaths
full_labels[11][11]
full_labels[11][12]
full_labels[12][11]
full_labels[12][12]

full_labels[11][11][(0.0, 800.0)]
full_labels[12][11][(292.0, 508.0)]