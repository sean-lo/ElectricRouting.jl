include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Test

args_df = CSV.read("$(@__DIR__)/../experiments/adaptive_ngroute_lmsr3/01a/args.csv", DataFrame)
row_index = 3

n_depots = args_df[row_index, :n_depots]
n_customers = args_df[row_index, :n_customers]
n_charging = args_df[row_index, :n_charging]
depot_pattern = String(args_df[row_index, :depot_pattern])
customer_pattern = String(args_df[row_index, :customer_pattern])
charging_pattern = String(args_df[row_index, :charging_pattern])
customer_spread = args_df[row_index, :customer_spread]
xmin = args_df[row_index, :xmin]
xmax = args_df[row_index, :xmax]
ymin = args_df[row_index, :ymin]
ymax = args_df[row_index, :ymax]
n_vehicles = args_df[row_index, :n_vehicles]
T = args_df[row_index, :T]
B = args_df[row_index, :B]
μ = args_df[row_index, :μ]
seed = args_df[row_index, :seed]
travel_cost_coeff = args_df[row_index, :travel_cost_coeff]
# charge_cost_coeff = args_df[row_index, :charge_cost_coeff]
charge_cost_coeff = 3
load_scale = args_df[row_index, :load_scale]
load_shape = args_df[row_index, :load_shape]
load_tolerance = args_df[row_index, :load_tolerance]
batch = args_df[row_index, :batch]
permissiveness = args_df[row_index, :permissiveness]

use_load = args_df[row_index, :use_load]
use_time_windows = args_df[row_index, :use_time_windows]

method = String(args_df[row_index, :method])
ngroute_neighborhood_charging_size = String(args_df[row_index, :ngroute_neighborhood_charging_size])
use_lmSR3_cuts = args_df[row_index, :use_lmSR3_cuts]
max_SR3_cuts = args_df[row_index, :max_SR3_cuts]

data = generate_instance(
    n_depots = n_depots,
    n_customers = n_customers,
    n_charging = n_charging,
    n_vehicles = n_vehicles,
    depot_pattern = depot_pattern,
    customer_pattern = customer_pattern,
    charging_pattern = charging_pattern,
    customer_spread = customer_spread,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    T = T,
    seed = seed,
    B = B,
    μ = μ,
    travel_cost_coeff = travel_cost_coeff,
    charge_cost_coeff = charge_cost_coeff,
    load_scale = load_scale,
    load_shape = load_shape,
    load_tolerance = load_tolerance,
    batch = batch,
    permissiveness = permissiveness,
    ;
    data_dir = "../../../data/",
)
graph = generate_graph_from_data(data)

@time SR3_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
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
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
);

@time lmSR3_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
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
    use_SR3_cuts = true,
    use_lmSR3_cuts = true,
);

using Plots
Plots.scatter(
    cumsum([x["implemented_SR3_cuts_count"] for x in lmSR3_results[5]]),
    [x["CG_sp_time_taken_mean"] for x in lmSR3_results[5]],
    xlabel = "Number of cuts",
    ylabel = "Pricing problem time taken (s)",
)
Plots.scatter(
    cumsum([x["implemented_SR3_cuts_count"] for x in SR3_results[5]]),
    [x["CG_sp_time_taken_mean"] for x in SR3_results[5]],
    xlabel = "Number of cuts",
    ylabel = "Pricing problem time taken (s)",
)


(method, use_lmSR3_cuts) = ("ours", false)
Env = nothing
time_windows = false
verbose = true
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

artificial_paths = generate_artificial_paths(data, graph)
some_paths = deepcopy(artificial_paths)
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