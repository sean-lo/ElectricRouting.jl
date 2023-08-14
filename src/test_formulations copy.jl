include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames

using Test

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 10,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 6,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)
plot_instance(data)

p_b_ngs_ch_LP_results, p_b_ngs_ch_IP_results, p_b_ngs_ch_params, p_b_ngs_ch_printlist, p_b_ngs_ch_some_paths, p_b_ngs_ch_model = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
p_o_ngs_ch_LP_results, p_o_ngs_ch_IP_results, p_o_ngs_ch_params, p_o_ngs_ch_printlist, p_o_ngs_ch_some_paths, p_o_ngs_ch_model = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
sp_b_ngs_ch_LP_results, sp_b_ngs_ch_IP_results, sp_b_ngs_ch_params, sp_b_ngs_ch_printlist, sp_b_ngs_ch_some_subpaths, sp_b_ngs_ch_some_charging_arcs, sp_b_ngs_ch_model = subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
sp_o_ngs_ch_LP_results, sp_o_ngs_ch_IP_results, sp_o_ngs_ch_params, sp_o_ngs_ch_printlist, sp_o_ngs_ch_some_subpaths, sp_o_ngs_ch_some_charging_arcs, sp_o_ngs_ch_model = subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", christofides = true, verbose = true);


collect_path_solution_metrics!(p_b_ngs_ch_LP_results, data, graph, p_b_ngs_ch_some_paths)
collect_path_solution_metrics!(p_o_ngs_ch_LP_results, data, graph, p_o_ngs_ch_some_paths)
collect_subpath_solution_metrics!(sp_b_ngs_ch_LP_results, data, sp_b_ngs_ch_some_subpaths, sp_b_ngs_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngs_ch_LP_results, data, sp_o_ngs_ch_some_subpaths, sp_o_ngs_ch_some_charging_arcs)

## Plotting

plot_path_solution(
    p_b_ngs_ch_LP_results, 
    data, graph, 
    p_b_ngs_ch_some_paths,
)
plot_path_solution(
    p_o_ngs_ch_LP_results, 
    data, graph, 
    p_o_ngs_ch_some_paths,
)
plot_subpath_solution(
    sp_b_ngs_ch_LP_results, 
    data, 
    sp_b_ngs_ch_some_subpaths, 
    sp_b_ngs_ch_some_charging_arcs,
)
plot_subpath_solution(
    sp_o_ngs_ch_LP_results, 
    data, 
    sp_o_ngs_ch_some_subpaths, 
    sp_o_ngs_ch_some_charging_arcs,
)

## Separation of SR3 inequalities (i)

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


## Separation of WSR3 inequalities (iii)

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
                            if a1[1] == ca[1] && a1[2] != ca[2] && a2[2] == ca[2]
                    )
                )
            end
        end
        push!(S_list, (violation, S..., included_charging_stations...))
    end
    sort!(S_list, by = x -> x[1], rev = true)
    return S_list
end


## Separation of WWSR3 inequalities (ii)

function enumerate_violated_path_WWSR3_inequalities(
    paths::Vector{Tuple{Float64, Path}},
    graph::EVRPGraph,
)
    S_list = Tuple{Float64, NTuple{3, Int}}[]
    for S in combinations(graph.N_customers, 3)
        violation = sum(
            [
                val 
                for (val, p) in values(paths)
                    if !isdisjoint(p.arcs, Tuple.(permutations(S, 2)))
            ],
            init = 0.0
        ) - 1.0
        if violation > 1e-6
            push!(S_list, (violation, Tuple(S)))
        end
    end
    sort!(S_list, by = x -> x[1], rev = true)
    return S_list
end

p_b_ngs_ch_SR3_list = enumerate_violated_path_SR3_inequalities(p_b_ngs_ch_LP_results["paths"], data)
p_b_ngs_ch_WSR3_list = enumerate_violated_path_WSR3_inequalities(p_b_ngs_ch_LP_results["paths"], graph)
p_b_ngs_ch_WWSR3_list = enumerate_violated_path_WWSR3_inequalities(p_b_ngs_ch_LP_results["paths"], data)

p_o_ngs_ch_SR3_list = enumerate_violated_path_SR3_inequalities(p_o_ngs_ch_LP_results["paths"], data)
p_o_ngs_ch_WSR3_list = enumerate_violated_path_WSR3_inequalities(p_o_ngs_ch_LP_results["paths"], graph)
p_o_ngs_ch_WWSR3_list = enumerate_violated_path_WWSR3_inequalities(p_o_ngs_ch_LP_results["paths"], data)

sp_b_ngs_ch_SR3_list = enumerate_violated_path_SR3_inequalities(sp_b_ngs_ch_LP_results["paths"], data)
sp_b_ngs_ch_WSR3_list = enumerate_violated_path_WSR3_inequalities(sp_b_ngs_ch_LP_results["paths"], graph)
sp_b_ngs_ch_WWSR3_list = enumerate_violated_path_WWSR3_inequalities(sp_b_ngs_ch_LP_results["paths"], data)

sp_o_ngs_ch_SR3_list = enumerate_violated_path_SR3_inequalities(sp_o_ngs_ch_LP_results["paths"], data)
sp_o_ngs_ch_WSR3_list = enumerate_violated_path_WSR3_inequalities(sp_o_ngs_ch_LP_results["paths"], graph)
sp_o_ngs_ch_WWSR3_list = enumerate_violated_path_WWSR3_inequalities(sp_o_ngs_ch_LP_results["paths"], data)




## Adding inequalities to models

for (method, ngroute_alt) in [
    ("benchmark", false),
    ("benchmark", true),
    ("ours", false),
    ("ours", true),
]
    data = generate_instance(
        ;
        n_depots = 4,
        n_customers = 16,
        n_charging = 7,
        n_vehicles = 6,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "circular_packing",
        shrinkage_depots = 1.0,
        shrinkage_charging = 1.0,
        T = 40000,
        seed = 5,
        B = 15000,
        μ = 5,
        travel_cost_coeff = 7,
        charge_cost_coeff = 3,
        load_scale = 5.0,
        load_shape = 20.0,
        load_tolerance = 1.3,
        batch = 1,
        permissiveness = 0.2,
    )
    graph = generate_graph_from_data(data)
    LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = method, ngroute = true, ngroute_alt = ngroute_alt, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
end


data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 16,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 5,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)
plot_path_solution(
    LP_results,
    data, graph,
    some_paths,
)


data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 10,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 6,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

plot_path_solution(
    LP_results,
    data,
    graph,
    some_paths,
)

LP_results["paths"][1]

# LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = false, subpath_single_service = false, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

p = LP_results["paths"][1][2]

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);


plot_path_solution(
    LP_results,
    data, graph,
    some_paths,
)
LP_results["paths"]





data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 10,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 6,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

plot_path_solution(
    LP_results,
    data, graph,
    some_paths,
)

T1 = Dict{Int64, Dict{Int64, Dict{Tuple{Vararg{Int64}}, SortedDict{Tuple{Int64, Int64}, PathLabel, Base.Order.ForwardOrdering}}}}

T2 = Dict{Int, Dict{Int, T3}} where {T3 <: AbstractDict}

T1 <: T2
















data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 10,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 6,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)

some_paths = generate_artificial_paths(data, graph)
path_costs = compute_path_costs(
    data, graph,
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)
neighborhoods = compute_ngroute_neighborhoods(
    graph,
    Int(ceil(sqrt(graph.n_customers))); 
    charging_depots_size = "small",
)
printlist = String[]
converged = false

model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
    ;
)

WSR3_constraints = Dict{Tuple{Vararg{Int}}, ConstraintRef}()
WSR3_list = Tuple{Float64, Vararg{Int}}[]

CGLP_results, CG_params = path_formulation_column_generation!(
    model, z, WSR3_constraints,
    data, graph, 
    some_paths, path_costs, path_service,
    printlist,
    ;
    method = "ours",
    christofides = true,
    neighborhoods = neighborhoods,
    ngroute = true,
    ngroute_alt = false,
    verbose = true,
)

CGIP_results = path_formulation_solve_integer_model!(
    model,
    z,
)

CGLP_all_results = []
CGIP_all_results = []
push!(CGLP_all_results, CGLP_results)
push!(CGIP_all_results, CGIP_results)

CGLP_results["objective"]
CGIP_results["objective"]

CGLP_results["paths"] = collect_path_solution_support(
    CGLP_results, some_paths, data, graph,
)
plot_path_solution(
    CGLP_results, 
    data, graph, some_paths,
)
generated_WSR3_list = enumerate_violated_path_WSR3_inequalities(CGLP_results["paths"], graph)


append!(WSR3_list, generated_WSR3_list)
# Add violated inequalities to master problem
add_WSR3_constraints_to_path_model!(
    model, z, some_paths, 
    WSR3_constraints, generated_WSR3_list, 
)

generated_WSR3_to_charging_extra_map = generate_new_WSR3_to_charging_extra_map(
    generated_WSR3_list,
    graph,
)
augment_graph_with_WSR3_duals!(
    graph,
    generated_WSR3_to_charging_extra_map,
)
neighborhoods = augment_neighborhoods_extra_with_WSR3_duals(
    neighborhoods,
    generated_WSR3_to_charging_extra_map,
)
neighborhoods

# debug 

p = CGLP_all_results[2]["paths"][1][2]
compute_path_modified_cost(
    data,
    graph,
    p,
    CGLP_all_results[1]["κ"],
    CGLP_all_results[1]["μ"],
    CGLP_all_results[1]["ν"],
)

b = generate_base_labels_ngroute_sigma(
    data,
    graph,
    neighborhoods,
    CGLP_all_results[1]["κ"],
    CGLP_all_results[1]["μ"],
    CGLP_all_results[1]["ν"],
    Dict{Tuple{Vararg{Int}}, Float64}(),
    ;
    christofides = true,
)
p.arcs
b[18][21][7]
b[18][18][15]
b[21][21][7]

s1 = p.subpaths[2]
s2 = b[21][21][7][(21,)][(44530,)]
compute_subpath_modified_cost(
    data,
    graph,
    s1,
    CGLP_all_results[1]["κ"],
    CGLP_all_results[1]["μ"],
    CGLP_all_results[1]["ν"],
)
s2.cost


function generate_new_WSR3_to_charging_extra_map(
    generated_WSR3_list::Vector{Tuple{Float64, Vararg{Int}}},
    graph::EVRPGraph,
)
    WSR3_to_charging_extra_map = Dict{NTuple{4, Int}, Int}()
    count = graph.n_nodes_extra
    for WSR3 in generated_WSR3_list
        S = Tuple(WSR3[2:4])
        for cs in WSR3[5:end]
            if (S..., cs) in values(graph.charging_extra_to_WSR3_map)
                continue
            end
            count += 1
            WSR3_to_charging_extra_map[(S..., cs)] = count
        end
    end
    return WSR3_to_charging_extra_map
end

function augment_neighborhoods_with_WSR3_duals(
    neighborhoods::BitMatrix,
    generated_WSR3_to_charging_extra_map::Dict{NTuple{4, Int}, Int},
)
    new_neighborhoods = copy(neighborhoods)

    for (key, _) in pairs(generated_WSR3_to_charging_extra_map)
        S = Tuple(key[1:3])
        cs = key[4]

        # add cs and S to the neighborhood of i (for i in S)
        for i in S
            new_neighborhoods[i, cs] = true
            new_neighborhoods[i, collect(S)] .= true
        end
    end
    return new_neighborhoods
end

function augment_neighborhoods_extra_with_WSR3_duals(
    neighborhoods::BitMatrix,
    generated_WSR3_to_charging_extra_map::Dict{NTuple{4, Int}, Int},
)
    n_new_nodes = length(generated_WSR3_to_charging_extra_map)
    n_total_nodes = size(neighborhoods, 1)
    n_total_nodes_new = size(neighborhoods, 1) + n_new_nodes
    new_neighborhoods = falses(n_total_nodes_new, n_total_nodes_new)
    new_neighborhoods[1:n_total_nodes, 1:n_total_nodes] .= neighborhoods

    for (key, cs_new) in pairs(generated_WSR3_to_charging_extra_map)
        S = Tuple(key[1:3])
        cs = key[4]

        # add cs and S to the neighborhood of i (for i in S)
        # for i in S
        #     new_neighborhoods[i, cs] = true
        #     new_neighborhoods[i, cs_new] = true
        #     new_neighborhoods[i, collect(S)] .= true
        # end
        new_neighborhoods[cs_new, cs_new] = true
    end
    return new_neighborhoods
end

function augment_graph_with_WSR3_duals!(
    graph::EVRPGraph,
    generated_WSR3_to_charging_extra_map::Dict{NTuple{4, Int}, Int},
)
    n_new_nodes = length(generated_WSR3_to_charging_extra_map)
    n_total_nodes = graph.n_nodes_extra
    n_total_nodes_new = graph.n_nodes_extra + n_new_nodes
    
    c = zeros(Float64, (n_total_nodes_new, n_total_nodes_new))
    c[1:n_total_nodes, 1:n_total_nodes] .= graph.c
    t = zeros(Float64, (n_total_nodes_new, n_total_nodes_new))
    t[1:n_total_nodes, 1:n_total_nodes] .= graph.t
    q = zeros(Float64, (n_total_nodes_new, n_total_nodes_new))
    q[1:n_total_nodes, 1:n_total_nodes] .= graph.q

    count = n_total_nodes
    add_vertices!(graph.G, n_new_nodes)

    for (key, cs_new) in pairs(generated_WSR3_to_charging_extra_map)
        S = Tuple(key[1:3])
        cs = key[4]

        # initialize new vertex
        graph.node_labels[cs_new] = graph.node_labels[cs]
        graph.charging_extra_to_WSR3_map[cs_new] = (S..., cs)
        graph.nodes_extra_to_nodes_map[cs_new] = cs
            
        # remove edges cs -> i (for i in S)
        for i in S
            rem_edge!(graph.G, cs, i)
            delete!(graph.A, (cs, i))
        end

        # add edges i -> cs_new (for i in all nodes except cs)
        for i in setdiff(graph.N_nodes, cs) # do not include extra CS
            add_edge!(graph.G, i, cs_new)
            push!(graph.A, (i, cs_new))
            c[i, cs_new] = graph.c[i, cs]
            t[i, cs_new] = graph.t[i, cs]
            q[i, cs_new] = graph.q[i, cs]
        end

        # add edges cs_new -> i (for i in S)
        for i in S
            add_edge!(graph.G, cs_new, i)
            push!(graph.A, (cs_new, i))
            c[cs_new, i] = graph.c[cs, i]
            t[cs_new, i] = graph.t[cs, i]
            q[cs_new, i] = graph.q[cs, i]
        end
    end

    merge!(graph.WSR3_to_charging_extra_map, generated_WSR3_to_charging_extra_map)

    append!(graph.N_charging_extra, collect(n_total_nodes+1:n_total_nodes_new))
    append!(graph.N_depots_charging_extra, collect(n_total_nodes+1:n_total_nodes_new))    
    append!(graph.N_nodes_extra, collect(n_total_nodes+1:n_total_nodes_new))
    append!(graph.α, zeros(Int, n_new_nodes))
    append!(graph.β, fill(graph.T, n_new_nodes))

    graph.n_charging_extra += n_new_nodes
    graph.n_depots_charging_extra += n_new_nodes
    graph.n_nodes_extra += n_new_nodes

    graph.min_t = dijkstra_shortest_paths(graph.G, graph.N_depots, t).dists
    graph.min_q = dijkstra_shortest_paths(graph.G, graph.N_depots_charging_extra, q).dists

    graph.c = c
    graph.t = t
    graph.q = q

    return graph
end

CGLP_results["objective"]
optimize!(model)
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


κ = CGLP_results["κ"]
μ = CGLP_results["μ"]
ν = CGLP_results["ν"]
σ = CGLP_results["σ"]

base_labels_result = @timed generate_base_labels_ngroute_sigma(
    data, graph, neighborhoods, κ, μ, ν, σ,
    ;
    christofides = true,
    time_limit = Inf,
);
base_labels_result.time

full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_sigma(
    data, graph, neighborhoods, base_labels_result.value, 
    κ, μ,
    ;
    christofides = true,
);
full_labels_result.time

negative_full_labels = get_negative_path_labels_from_path_labels(full_labels_result.value)
normalize_path_labels!(negative_full_labels, graph)

negative_full_labels[1]
generated_paths = get_paths_from_negative_path_labels(
    graph, negative_full_labels,
)

[
    compute_path_modified_cost(
        data, graph, p,
        κ, μ, ν,
    )
    for p_list in values(generated_paths)
        for p in p_list
]

mp_constraint_time = add_paths_to_path_model!(
    model,
    z,
    some_paths, 
    path_costs,
    path_service,
    generated_paths,
    data, graph,
)

CGLP_results["paths"] = collect_path_solution_support(
    CGLP_results, 
    some_paths, 
    data, graph
)
CGLP_results["paths"]

CGIP_results = path_formulation_solve_integer_model!(
    model,
    z,
)

push!(CGLP_all_results, CGLP_results)
push!(CGIP_all_results, CGIP_results)

[r["objective"] for r in CGLP_all_results]
[r["objective"] for r in CGIP_all_results]


compute_path_modified_cost(data, graph, p, CGLP_all_results[1]["κ"], CGLP_all_results[1]["μ"], CGLP_all_results[1]["ν"])

b = generate_base_labels_ngroute_sigma(
    data, 
    graph, 
    neighborhoods,
    CGLP_all_results[1]["κ"], CGLP_all_results[1]["μ"], CGLP_all_results[1]["ν"], 
    Dict{Tuple{Vararg{Int}}, Float64}(),
    ;
    christofides = true,
)

b
p # this is an optimal path from 2nd iteration
p.arcs

b[18][21][7]
b[21][21][7]