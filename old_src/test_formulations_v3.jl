using Pkg
Pkg.activate(".")
include("nonlinear_charging.jl")
include("path_formulation.jl")
include("utils.jl")

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 32,
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
    T = 75000,
    seed = 1,
    B = 15000,
    inverse_refueling_rate = 5.0,
    travel_cost_coeff = 7,
    charge_cost_coeff = 0,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)
plot_instance(data)
data.α .= 0
data.β .= data.T
graph.α .= 0
graph.β .= data.T

include("subpath_stitching_v3.jl")

ours_results = path_formulation_column_generation_with_adaptive_ngroute_SR3_cuts(
    data,
    graph,
    ;
    method = "ours",
    use_time_windows = false,
    use_load = true,
    use_heuristic = false,
    use_nonlinear_charging = false,
    use_ngroute = true,
    use_adaptive_ngroute = false,
    use_SR3_cuts = false,
    use_lmSR3_cuts = false,
    verbose = true,
);

ours_results_TW = path_formulation_column_generation_with_adaptive_ngroute_SR3_cuts(
    data,
    graph,
    ;
    method = "ours",
    use_time_windows = true,
    use_load = true,
    use_heuristic = false,
    use_nonlinear_charging = false,
    use_ngroute = true,
    use_adaptive_ngroute = false,
    use_SR3_cuts = false,
    use_lmSR3_cuts = false,
    verbose = true,
);

no_TW_paths = ours_results[1][1]["paths"]
TW_paths = ours_results_TW[1][1]["paths"]

print.(ours_results[6]);
print.(ours_results_TW[6]);

graph.min_t
