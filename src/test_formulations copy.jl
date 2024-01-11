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
sp_b_ngs_ch_LP_results, sp_b_ngs_ch_IP_results, sp_b_ngs_ch_params, sp_b_ngs_ch_printlist, sp_b_ngs_ch_some_subpaths, sp_b_ngs_ch_some_charging_arcs, sp_b_ngs_ch_model = subpath_formulation_column_generation_integrated_from_paths(data, graph; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
sp_o_ngs_ch_LP_results, sp_o_ngs_ch_IP_results, sp_o_ngs_ch_params, sp_o_ngs_ch_printlist, sp_o_ngs_ch_some_subpaths, sp_o_ngs_ch_some_charging_arcs, sp_o_ngs_ch_model = subpath_formulation_column_generation_integrated_from_paths(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", christofides = false, verbose = true);


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

for (method) in [
    ("benchmark"),
    ("ours"),
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
    LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = method, ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
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

for (method, ngroute, christofides) in [
    # ("benchmark", true, false),
    # ("benchmark", true, true),
    ("ours", true, false),
    # ("ours", true, true),
]
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
    LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(
        data, graph; 
        method = method, 
        ngroute = ngroute, 
        ngroute_neighborhood_charging_depots_size = "small", 
        verbose = true, 
        christofides = christofides,
    )
end

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);


# LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = false, subpath_single_service = false, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = false);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = false);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_cuts(data, graph; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);


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
    n_customers = 15,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 1,
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

LP_results, IP_results, CG_params, printlist, some_paths, model, z, WSR3_constraints = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(data, graph; method = "ours", ngroute_neighborhood_charging_depots_size = "medium", verbose = true);

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 15,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 1,
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
    charging_depots_size = "medium",
)
printlist = String[]
converged = false

model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
    ;
)

SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
SR3_list = Tuple{Float64, NTuple{3, Int}}[]
WSR3_constraints = Dict{Tuple{Vararg{Int}}, ConstraintRef}()
WSR3_list = Tuple{Float64, Vararg{Int}}[]
CGLP_all_results = []
CGIP_all_results = []
# Start of while loop

CGLP_results, CG_params = path_formulation_column_generation!(
    model, z, WSR3_constraints, SR3_constraints,
    data, graph, 
    some_paths, path_costs, path_service,
    printlist,
    ;
    method = "ours",
    christofides = false,
    neighborhoods = neighborhoods,
    ngroute = true,
    verbose = true,
)

CGIP_results = path_formulation_solve_integer_model!(
    model,
    z,
)


push!(CGLP_all_results, CGLP_results)
push!(CGIP_all_results, CGIP_results)

CGLP_results["objective"]
CGIP_results["objective"]
CG_params["LP_IP_gap"] = 1.0 - CGLP_results["objective"] / CGIP_results["objective"]


CGLP_results["paths"] = collect_path_solution_support(
    CGLP_results, some_paths, data, graph,
)
plot_path_solution(
    CGLP_results, 
    data, graph, some_paths,
)

cycles_lookup = detect_cycles_in_path_solution([p for (val, p) in CGLP_results["paths"]], graph)
delete_paths_with_found_cycles_from_model!(model, z, some_paths, path_costs, path_service, cycles_lookup, graph)
modify_neighborhoods_with_found_cycles!(neighborhoods, cycles_lookup)


generated_WSR3_list = enumerate_violated_path_WSR3_inequalities(
    CGLP_results["paths"], 
    graph,
)


SR3_constraints = Dict{
    NTuple{3, Int},
    ConstraintRef
}()
generated_SR3_list = Vector{Tuple{Float64, NTuple{3, Int}}}[]

function check_path_in_SR3_constraint(
    path::Path,
    S::NTuple{3, Int},
)
    return Int(floor(sum(path.served[collect(S)]) / 2))
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
    charging_depots_size = "medium",
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
    christofides = false,
    neighborhoods = neighborhoods,
    ngroute = true,
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
CG_params["LP_IP_gap"] = 1.0 - CGLP_results["objective"] / CGIP_results["objective"]


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


# debug 

# p = CGLP_all_results[2]["paths"][1][2]
# compute_path_modified_cost(
#     data,
#     graph,
#     p,
#     CGLP_all_results[1]["κ"],
#     CGLP_all_results[1]["μ"],
#     CGLP_all_results[1]["ν"],
# )

# p.arcs
# b[18][21][7]
# b[18][18][15]
# b[21][21][7]

# s1 = p.subpaths[2]
# s2 = b[21][21][7][(21,)][(44530,)]
# compute_subpath_modified_cost(
#     data,
#     graph,
#     s1,
#     CGLP_all_results[1]["κ"],
#     CGLP_all_results[1]["μ"],
#     CGLP_all_results[1]["ν"],
# )
# s2.cost

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
)


κ = CGLP_results["κ"]
μ = CGLP_results["μ"]
ν = CGLP_results["ν"]

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
    WSR3_constraints,
    data, graph,
)

WSR3_constraints[(4, 5, 6, 21)]

list_of_constraint_types(model, )
all_constraints(model, AffExpr, MOI.LessThan{Float64})

generated_paths

WSR3_constraints

p = generated_paths[(13, 0, 15000), (11, 136962, 0)][1]

S = (4, 5, 6)
intersect(p.arcs, Tuple.(permutations(S, 2)))

function check_path_in_WSR3_constraint(
    path::Path,
    S::NTuple{3, Int},
    included_charging_stations::Tuple{Vararg{Int}},
)
    if length(intersect(path.arcs, Tuple.(permutations(S, 2)))) ≥ 1
        return true
    end
    if length(intersect(path.customer_arcs, Tuple.(permutations(S, 2)))) == 0
        return false
    end
    for cs in included_charging_stations
        if length(intersect(path.arcs, Tuple.(permutations(vcat(S, cs), 2)))) ≥ 2 
            return true
        end
    end
    return false
end

for (key, p_list) in pairs(generated_paths)
    for p in p_list
        println(p.arcs)
        println(p.customer_arcs)
        println(check_path_in_WSR3_constraint(p, (4, 5, 6), (21,)))
    end
end

T

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
