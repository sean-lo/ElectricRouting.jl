include("path_formulation.jl")
include("utils.jl")

using BenchmarkTools, Test

data = generate_instance(
    n_depots = 4,
    n_customers = 40,
    n_charging = 21,
    n_vehicles = 8,
    depot_pattern = "grid",
    customer_pattern = "random_box",
    charging_pattern = "grid_clipped",
    customer_spread = 1e-3,
    xmin = 0.0,
    xmax = 4.0,
    ymin = 0.0,
    ymax = 4.0,
    T = 100000,
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
    ;
    # data_dir = "../../../data/",
)
graph = generate_graph_from_data(data, sparse = true, sparse_linear_max_q = 2.0 * data.B)
plot_instance(data, plot_edges = true, graph = graph)
length(edges(graph.G))

graph = generate_graph_from_data(data, sparse = false)
length(edges(graph.G))

using Profile

@timev SR3_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(graph.n_customers)),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
);