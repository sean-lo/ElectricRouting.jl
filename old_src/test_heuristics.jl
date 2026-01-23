using Pkg
Pkg.activate("$(@__DIR__)/..")

include("utils.jl")

using CSV, DataFrames

args_df = CSV.read("$(@__DIR__)/../experiments/ours_benchmark/01a/args.csv", DataFrame)
row_index = 5000

data = generate_instance(
    n_depots = 4,
    n_customers = 24,
    n_charging = 11,
    n_vehicles = 6,
    depot_pattern = "grid",
    customer_pattern = "random_box",
    charging_pattern = "grid_clipped",
    customer_spread = 0.1,
    xmin = 0.0,
    xmax = 4.0,
    ymin = 0.0,
    ymax = 2.0,
    T = 90000,
    seed = 20,
    B = 15000,
    inverse_refueling_rate = 5.0,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.4,
    ;
    data_dir = "../../../data/",
)
graph = generate_graph_from_data(data)
plot_instance(data, plot_edges = true, graph = graph)

include("path_formulation.jl")
include("subpath_stitching_v2.jl")

@timev ours_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    use_load = true,
    use_heuristic = false,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
    time_limit = 120.0,
);


@timev ours_results_h = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    use_load = true,
    use_heuristic = true,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
    time_limit = 120.0,
);



include("read_evrptw_data.jl")
fp = "data/evrptw/Instances/c101_21_25.txt"

# data_u = read_evrptw_instance(
#     fp,
#     5, # n_vehicles
#     7, # travel_cost_coeff
#     3, # charge_cost_coeff
#     # n_charging = 10,
#     # n_customers = 10,
# )
# graph_u = generate_graph_from_data(data_u)

data_e = read_evrptw_instance(
    fp,
    5, # n_vehicles
    7, # travel_cost_coeff
    3, # charge_cost_coeff
    # n_charging = 10,
    # n_customers = 10,
    scale_time_horizon = 0.6,
    scale_charge_capacity = 0.2,
    scale_load_capacity = 0.6,
)
graph_e = generate_graph_from_data(data_e)



@timev ours_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data_e, graph_e,
    ;
    method = "ours",
    use_load = true,
    use_heuristic = false,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph_e.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
    time_limit = 60.0,
);


@timev ours_results_h = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data_e, graph_e,
    ;
    method = "ours",
    use_load = true,
    use_heuristic = true,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph_e.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
    time_limit = 60.0,
);
