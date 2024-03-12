include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Cthulhu
using Test

(
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
    4, "grid", "random_box", "grid_clipped", 
    0.1, 0.0, 0.0, 2.0,
    15000, 5
)

xmax_range = [2.0, 3.0, 4.0]
n_customers_range = [12, 15, 18, 21]
k_range = [
    2.5, 3.0, 3.5,
]

xmax_k_bad_ranges = [
    (1.0, 1.5),
    (2.0, 1.5),
    (2.0, 2.0),
    (3.0, 1.5), 
    (3.0, 2.0),
]

(
    xmax, n_customers, k, seed
) = (3.0, 12, 2.0, 4)

n_charging = Int((xmax - xmin + 1)*(ymax - ymin + 1) - 4)
T = Int(B * k * (μ + 1) / μ)
# n_vehicles = n_customers / ((k / (xmax * ymax)) * 2 + 0.2)
n_vehicles = 6


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
    seed = 25,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)
p = plot_instance(data)


CGLP_all_results = nothing
CGIP_all_results = nothing
some_paths = nothing
CG_all_params = nothing
for (method, elementary, ngroute) in [
    # ("ours", true, false)
    # ("ours", false, true)
    # ("ours", false, false)
    # ("benchmark", true, false)
    # ("benchmark", false, true)
    ("benchmark", false, false)
]
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
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
        max_iters = 20.0,
    );
    println("$(all_params[end]["CGLP_objective"])\t$(all_params[end]["CGIP_objective"])")
    collect_solution_metrics!(CGLP_all_results[end], data)
    collect_solution_metrics!(CGIP_all_results[end], data)

    @printf("%6.4f\t%6.4f\n", CGLP_all_results[end]["mean_path_length"], CGIP_all_results[end]["mean_path_length"])
    @printf("%6.4f\t%6.4f\n", CGLP_all_results[end]["mean_subpath_length"], CGIP_all_results[end]["mean_subpath_length"])
    @printf("%6.4f\t%6.4f\n", CGLP_all_results[end]["mean_ps_length"], CGIP_all_results[end]["mean_ps_length"])
    @printf("%6.4f\t%6.4f\n", CGLP_all_results[end]["weighted_mean_path_length"], CGIP_all_results[end]["weighted_mean_path_length"])
    @printf("%6.4f\t%6.4f\n", CGLP_all_results[end]["weighted_mean_subpath_length"], CGIP_all_results[end]["weighted_mean_subpath_length"])
    @printf("%6.4f\t%6.4f\n", CGLP_all_results[end]["weighted_mean_ps_length"], CGIP_all_results[end]["weighted_mean_ps_length"])
    println()
end 

CG_all_params[1]

CGIP_all_results

CGLP_all_results[end]["utilization"]
CGIP_all_results[end]["utilization"]

total_subpath_length = 0.0
total_subpath_ncust = 0.0
num_subpaths = 0.0
total_path_length = 0.0
total_path_ncust = 0.0
num_paths = 0.0
total_ps_length = 0.0
for path_l in values(some_paths)
    for p in path_l
        if (
            length(p.subpaths) == 1 
            && (
                p.subpaths[1].artificial # artificial path
                # || length(p.subpaths[1].arcs) == 1 # path from depot to depot
            )
        )
            continue
        end
        total_subpath_length += sum(sum(s.served) + 1 for s in p.subpaths) 
        total_subpath_ncust += sum(sum(s.served) for s in p.subpaths)

        num_subpaths += length(p.subpaths)
        
        total_path_length += sum(p.served) + length(p.subpaths)
        total_path_ncust += sum(p.served)

        num_paths += 1

        total_ps_length += length(p.subpaths)
    end
end
total_subpath_length / num_subpaths

total_path_length / num_paths

total_ps_length / num_paths

some_paths



include("../experiments/ours_benchmark/01/script.jl")
