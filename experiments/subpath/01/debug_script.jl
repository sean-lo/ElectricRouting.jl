include("../../../src/subpath_formulation.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

row_index = 683
# Get paramters from args_df at row row_index
n_depots = args_df[row_index, :n_depots]
n_customers = args_df[row_index, :n_customers]
n_charging = args_df[row_index, :n_charging]
n_vehicles = args_df[row_index, :n_vehicles]
depot_pattern = String(args_df[row_index, :depot_pattern])
customer_pattern = String(args_df[row_index, :customer_pattern])
charging_pattern = String(args_df[row_index, :charging_pattern])
shrinkage_depots = args_df[row_index, :shrinkage_depots]
shrinkage_charging = args_df[row_index, :shrinkage_charging]
T = args_df[row_index, :T]
B = args_df[row_index, :B]
seed = args_df[row_index, :seed]
μ = args_df[row_index, :μ]
travel_cost_coeff = args_df[row_index, :travel_cost_coeff]
charge_cost_coeff = args_df[row_index, :charge_cost_coeff]
load_scale = args_df[row_index, :load_scale]
load_shape = args_df[row_index, :load_shape]
load_tolerance = args_df[row_index, :load_tolerance]
batch = args_df[row_index, :batch]
permissiveness = args_df[row_index, :permissiveness]

use_load = args_df[row_index, :use_load]
use_time_windows = args_df[row_index, :use_time_windows]
method = String(args_df[row_index, :method])
subpath_single_service = args_df[row_index, :subpath_single_service]
subpath_check_customers = args_df[row_index, :subpath_check_customers]
path_single_service = args_df[row_index, :path_single_service]
path_check_customers = args_df[row_index, :path_check_customers]
check_customers_accelerated = args_df[row_index, :check_customers_accelerated]

data = generate_instance(
    n_depots = n_depots,
    n_customers = n_customers,
    n_charging = n_charging,
    n_vehicles = n_vehicles,
    depot_pattern = depot_pattern,    
    customer_pattern = customer_pattern,
    charging_pattern = charging_pattern,
    shrinkage_depots = shrinkage_depots,
    shrinkage_charging = shrinkage_charging,
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
)
G = construct_graph(data)

(
    r_LP_results, r_IP_results, r_params, r_printlist, r_subpaths, r_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, 
    method = method, 
    subpath_single_service = subpath_single_service,
    subpath_check_customers = subpath_check_customers,
    path_single_service = path_single_service,
    path_check_customers = path_check_customers,
    check_customers_accelerated = check_customers_accelerated,
)
r_LP_results["subpaths"], r_LP_results["charging_arcs"] = collect_solution_support(r_LP_results, r_subpaths, r_charging_arcs)

results_subpaths, results_charging_arcs = collect_solution_support(r_LP_results, r_subpaths, r_charging_arcs)

all_paths = Tuple{Float64, Path}[]

while true
# begin
    remaining_vals = 0.0
    for x in results_subpaths
        remaining_vals += x[1]
    end
    for x in results_charging_arcs
        remaining_vals += x[1]
    end
    remaining_vals
    if remaining_vals < 1e-3
        break
        # println("oh no")
    end
    pathlist = []
    current_states = [
        (depot, 0.0, data["B"]) for depot in data["N_depots"]
    ]
    while true
    # begin
        i = findfirst(
            s -> (
                (s[2].starting_node, s[2].starting_time, s[2].starting_charge) 
                in current_states
                && s[1] > 1e-4
            ),
            results_subpaths
        )
        (val, s) = results_subpaths[i]
        push!(pathlist, (val, i, s))
        if s.current_node in data["N_depots"]
            break
            # println("oh no")
        end
        current_states = [(s.current_node, s.current_time, s.current_charge)]
        i = findfirst(
            s -> (
                (s[2].starting_node, s[2].starting_time, s[2].starting_charge) 
                in current_states
                && s[1] > 1e-4
            ),
            results_charging_arcs
        )
        (val, a) = results_charging_arcs[i]
        push!(pathlist, (val, i, a))
        current_states = [(a.starting_node, a.current_time, a.current_charge)]
    end
    pathlist
    minval = minimum(x[1] for x in pathlist)
    for (val, i, s) in pathlist
        if typeof(s) == Subpath
            newval = results_subpaths[i][1] - minval
            if abs.(newval) < 1e-6
                newval = 0.0
            end
            results_subpaths[i] = (
                newval,
                results_subpaths[i][2],
            )
        else
            newval = results_charging_arcs[i][1] - minval
            if abs.(newval) < 1e-6
                newval = 0.0
            end
            results_charging_arcs[i] = (
                newval,
                results_charging_arcs[i][2],
            )
        end
    end
    path = Path(
        subpaths = [x[3] for (i, x) in enumerate(pathlist) if i % 2 == 1],
        charging_arcs = [x[3] for (i, x) in enumerate(pathlist) if i % 2 == 0],
    )
    push!(all_paths, (minval, path))
end

all_paths
results_subpaths
results_charging_arcs
sum(x[1] for x in results_subpaths)
sum(x[1] for x in results_charging_arcs)



return all_paths




    collect_solution_metrics!(r_LP_results, data, r_subpaths, r_charging_arcs)
    collect_solution_metrics!(r_IP_results, data, r_subpaths, r_charging_arcs)
    records = [
        (
            n_depots = n_depots,
            n_customers = n_customers,
            n_charging = n_charging,
            n_vehicles = n_vehicles,
            depot_pattern = depot_pattern,    
            customer_pattern = customer_pattern,
            charging_pattern = charging_pattern,
            shrinkage_depots = shrinkage_depots,
            shrinkage_charging = shrinkage_charging,
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
            use_load = use_load,
            use_time_windows = use_time_windows,
            method = method,
            subpath_single_service = subpath_single_service,
            subpath_check_customers = subpath_check_customers,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            check_customers_accelerated = check_customers_accelerated,
            LP_objective = r_LP_results["objective"],
            IP_objective = r_IP_results["objective"],
            LP_IP_gap = r_params["LP_IP_gap"],
            counter = r_params["counter"],
            converged = r_params["converged"],
            time_taken = r_params["time_taken"],
            sp_base_time_taken_total = r_params["sp_base_time_taken_total"],
            sp_full_time_taken_total = r_params["sp_full_time_taken_total"],
            lp_relaxation_time_taken_total = r_params["lp_relaxation_time_taken_total"],
            sp_base_time_taken_mean = r_params["sp_base_time_taken_mean"],
            sp_full_time_taken_mean = r_params["sp_full_time_taken_mean"],
            lp_relaxation_time_taken_mean = r_params["lp_relaxation_time_taken_mean"],
        )
    ]
    CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
    println("Row $row_index processed: $(args_df[row_index, :])")
end