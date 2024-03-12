include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Cthulhu
using Test

(
    n_vehicles,
    n_depots,
    depot_pattern,
    customer_pattern,
    charging_pattern,
    customer_spread,
    xmin,
    ymin,
    ymax,
    B,
    μ,
) = (
    6, 4, "grid", "random_box", "grid_clipped", 
    0.1, 0.0, 0.0, 2.0, 
    15000, 5,
)

setting_params = [
    # load
    # time windows
    (false, false),
]

method_params = [
    # method, ngroute_neighborhood_charging_size
    ("benchmark", "small"),
    ("benchmark", "medium"),
    ("ours", "small"),
    ("ours", "medium"),
]

data_params = [
    # xmax
    # k
    # density
    (4.0, 4.0, 2.5), # 20
    (4.0, 4.0, 3.0), # 24
    (4.0, 4.0, 3.5), # 28
    (4.0, 4.0, 4.0), # 32
    (4.0, 4.0, 4.5), # 36
]

instance_params = [
    # xmax
    # k
    # density
    # method
    # ngroute_neighborhood_charging_size
    # seed
    (4.0, 4.5, 4.5, "benchmark", "medium", 14)
    (4.0, 4.5, 4.5, "benchmark", "small", 14)
    (4.0, 5.0, 4.5, "benchmark", "medium", 14)
    (4.0, 5.0, 4.5, "benchmark", "small", 4)
    (4.0, 5.0, 4.5, "ours", "small", 18)
]

df = CSV.read("experiments/adaptive_ngroute/01/args.csv", DataFrame)
df[774, [:xmax, :T, :density, :method, :ngroute_neighborhood_charging_size, :seed]]
df[1174, [:xmax, :T, :density, :method, :ngroute_neighborhood_charging_size, :seed]]
df[734, [:xmax, :T, :density, :method, :ngroute_neighborhood_charging_size, :seed]]
df[1124, [:xmax, :T, :density, :method, :ngroute_neighborhood_charging_size, :seed]]

(xmax, T, density, method, ngroute_neighborhood_charging_size, seed) = df[774, [:xmax, :T, :density, :method, :ngroute_neighborhood_charging_size, :seed]]
method = String(method)
ngroute_neighborhood_charging_size = String(ngroute_neighborhood_charging_size)

data = generate_instance(
    ;
    n_depots = n_depots,
    n_customers = n_customers,
    n_charging = n_charging,
    depot_pattern = depot_pattern,    
    customer_pattern = customer_pattern,
    charging_pattern = charging_pattern,
    customer_spread = 0.2,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    n_vehicles = n_vehicles,
    T = T,
    B = B,
    μ = μ,
    seed = seed,
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
(
    CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
    some_paths, model, z, SR3_constraints
) = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = method,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size, 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
    use_lmSR3_cuts = false,
)




DataFrame(all_params)
CG_all_neighborhoods = nothing
for (xmax, k, density, method, ngroute_neighborhood_charging_size, seed) in instance_params
    n_customers = Int(density * (xmax - xmin) * (ymax - ymin))
    n_charging = Int((xmax - xmin + 1) * (ymax - ymin + 1) - 4)
    T = Int(B * k * (μ + 1) / μ)
    data = generate_instance(
        ;
        n_depots = n_depots,
        n_customers = n_customers,
        n_charging = n_charging,
        depot_pattern = depot_pattern,    
        customer_pattern = customer_pattern,
        charging_pattern = charging_pattern,
        customer_spread = customer_spread,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        n_vehicles = n_vehicles,
        T = T,
        B = B,
        μ = μ,
        seed = seed,
        travel_cost_coeff = 7,
        charge_cost_coeff = 3,
        load_scale = 5.0,
        load_shape = 20.0,
        load_tolerance = 1.3,
        batch = 1,
        permissiveness = 0.2,
    )
    graph = generate_graph_from_data(data)
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        elementary = false,
        ngroute = true,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size, 
        verbose = true,
        use_adaptive_ngroute = true,
        use_SR3_cuts = false,
        use_lmSR3_cuts = false,
    );
end

d = DataFrame(CG_all_params)
names(d)

select(
    d, 
    [
        # :lp_relaxation_time_taken_total,
        # :lp_relaxation_time_taken_mean,
        # :sp_base_time_taken_total,
        # :sp_base_time_taken_mean,
        # :sp_full_time_taken_total,
        # :sp_full_time_taken_mean,
        # :sp_time_taken_total,
        # :sp_time_taken_mean,
        :LP_IP_gap,
        :converged,
        :time_limit_reached,
        :time_taken,
    ]
)