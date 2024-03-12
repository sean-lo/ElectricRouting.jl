include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Cthulhu
using Test

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 20,
    n_charging = 9,
    n_vehicles = 6,
    depot_pattern = "grid",    
    customer_pattern = "random_box",
    charging_pattern = "grid",
    shrinkage_depots = 2.0,
    shrinkage_customers = 2.0,
    shrinkage_charging = 1.7,
    T = 80000,
    seed = 1,
    B = 30000,
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

for (method, elementary, ngroute) in [
    ("ours", true, false)
    ("ours", false, false)
    ("ours", false, true)
    ("benchmark", true, false)
    ("benchmark", false, false)
    ("benchmark", false, true)
]
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = @time @suppress path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        elementary = elementary,
        ngroute = ngroute,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        verbose = true,
        use_adaptive_ngroute = false,
        use_SR3_cuts = false,
        use_lmSR3_cuts = false,
    );
    println("$(CGLP_all_results[end]["objective"])\t$(CGIP_all_results[end]["objective"])")
end

(method, use_lmSR3_cuts) = ("ours", true)
Env = nothing
time_windows = false
ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers)))
ngroute_neighborhood_depots_size = "small"
ngroute_neighborhood_charging_size = "small"

SR3_constraints = Dict{
    Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
    ConstraintRef,
}()
SR3_list = Tuple{Float64, NTuple{3, Int}, Tuple{Vararg{Int}}}[]



neighborhoods = compute_ngroute_neighborhoods(
    graph,
    ngroute_neighborhood_size; 
    depots_size = ngroute_neighborhood_depots_size,
    charging_size = ngroute_neighborhood_charging_size,
)

some_paths = generate_artificial_paths(data, graph)
path_costs = compute_path_costs(
    data, graph, 
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)
printlist = String[]

model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
    ; 
    Env = Env,
)
CGLP_results, CG_params = path_formulation_column_generation!(
    model, z, SR3_constraints,
    data, graph,
    some_paths, path_costs, path_service,
    printlist,
    ;
    method = method,
    time_windows = time_windows,
    subpath_single_service = false,
    subpath_check_customers = false,
    path_single_service = false,
    path_check_customers = false,
    neighborhoods = neighborhoods,
    ngroute = true,
)

CGIP_results = path_formulation_solve_integer_model!(
    model,
    z,
)

path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = method,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = use_lmSR3_cuts,
);

for (method, use_lmSR3_cuts) in [
    # ("benchmark", false),
    # ("benchmark", true),
    # ("ours", false), 
    ("ours", true),
]
    @time path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        verbose = true,
        use_adaptive_ngroute = true,
        use_SR3_cuts = true,
        use_lmSR3_cuts = use_lmSR3_cuts,
        max_SR3_cuts = 20,
    );
end

path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = method,
    elementary = false,
    ngroute = false,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = use_lmSR3_cuts,
    max_SR3_cuts = 20,
);



