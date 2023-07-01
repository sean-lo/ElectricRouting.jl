include("../../../src/subpath_formulation.jl")
include("../../../src/path_formulation.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

using JuMP, Gurobi

sleep(rand() * 30.0)
const GRB_ENV = Gurobi.Env()

TIME_LIMIT = 3600.0

# simple test case to quickly compile 
begin
    sample_data = generate_instance(
        n_depots = 4,
        n_customers = 10,
        n_charging = 9,
        n_vehicles = 7,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "circular_packing",
        shrinkage_depots = 1.0,
        shrinkage_charging = 0.7,
        T = 40000,
        seed = 3,
        B = 15000,
        μ = 5,
        travel_cost_coeff = 7,
        charge_cost_coeff = 3,
        load_scale = 5.0,
        load_shape = 20.0,
        load_tolerance = 1.3,
        batch = 5,
        permissiveness = 0.7,
        data_dir = "../../../data/",
    )
    method_params = [
        # formulation
        # method
        # subpath_single_service
        # subpath_check_customers
        # path_single_service
        # path_check_customers
        # check_customers_accelerated
        # christofides
        # ngroute
        # ngroute_alt
        # ngroute_neighborhood_charging_depots_size
        ("path", "benchmark", false, false, false, false, false,  true, false, false, "none"),
        ("path", "benchmark", false, false,  true,  true, false, false, false, false, "none"),
        ("path", "benchmark", false, false, false, false, false,  true,  true, false, "small"),
        ("path", "benchmark", false, false, false, false, false,  true,  true, false, "large"),
        ("path", "benchmark", false, false, false, false, false,  true,  true,  true, "small"),
        ("path", "benchmark", false, false, false, false, false,  true,  true,  true, "large"),
        ("path", "ours", false, false, false, false, false,  true, false, false, "none"),
        ("path", "ours",  true,  true, false, false, false,  true, false, false, "none"),
        ("path", "ours",  true,  true,  true,  true, false,  true, false, false, "none"),
        ("path", "ours", false, false, false, false, false,  true,  true, false, "small"),
        ("path", "ours", false, false, false, false, false,  true,  true, false, "large"),
        ("path", "ours", false, false, false, false, false,  true,  true,  true, "small"),
        ("path", "ours", false, false, false, false, false,  true,  true,  true, "large"),
        ("subpath", "benchmark", false, false, false, false, false,  true, false, false, "none"),
        ("subpath", "benchmark", false, false,  true,  true, false, false, false, false, "none"),
        ("subpath", "benchmark", false, false,  true,  true,  true, false, false, false, "none"),
        ("subpath", "benchmark", false, false, false, false, false,  true,  true, false, "small"),
        ("subpath", "benchmark", false, false, false, false, false,  true,  true, false, "large"),
        ("subpath", "benchmark", false, false, false, false, false,  true,  true,  true, "small"),
        ("subpath", "benchmark", false, false, false, false, false,  true,  true,  true, "large"),
        ("subpath", "ours", false, false, false, false, false,  true, false, false, "none"),
        ("subpath", "ours",  true,  true, false, false, false,  true, false, false, "none"),
        ("subpath", "ours",  true,  true, false, false,  true,  true, false, false, "none"),
        ("subpath", "ours",  true,  true,  true,  true, false,  true, false, false, "none"),
        ("subpath", "ours",  true,  true,  true,  true,  true,  true, false, false, "none"),
        ("subpath", "ours", false, false, false, false, false,  true,  true, false, "small"),
        ("subpath", "ours", false, false, false, false, false,  true,  true, false, "large"),
        ("subpath", "ours", false, false, false, false, false,  true,  true,  true, "small"),
        ("subpath", "ours", false, false, false, false, false,  true,  true,  true, "large"),
    ]
    for method_param in method_params
        if method_param[1] == "path"
            (
                LP_results, IP_results, cgparams, printlist, paths
            ) = path_formulation_column_generation(
                sample_data,
                ;
                Env = GRB_ENV, 
                method = method_param[2],
                subpath_single_service = method_param[3],
                subpath_check_customers = method_param[4],
                path_single_service = method_param[5],
                path_check_customers = method_param[6],
                christofides = method_param[8],
                ngroute = method_param[9],
                ngroute_alt = method_param[10],
                ngroute_neighborhood_charging_depots_size = method_param[11],
                time_limit = TIME_LIMIT,
            )
        elseif method_param[1] == "subpath"
            (
                LP_results, IP_results, cgparams, printlist, subpaths, charging_arcs
            ) = subpath_formulation_column_generation_integrated_from_paths(
                sample_data,
                ;
                Env = GRB_ENV, 
                method = method_param[2],
                subpath_single_service = method_param[3],
                subpath_check_customers = method_param[4],
                path_single_service = method_param[5],
                path_check_customers = method_param[6],
                christofides = method_param[8],
                ngroute = method_param[9],
                ngroute_alt = method_param[10],
                ngroute_neighborhood_charging_depots_size = method_param[11],
                time_limit = TIME_LIMIT,
            )
        end
    end
end

println("Compilation complete.")

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

task_index = parse(Int, ARGS[1]) + 1
n_tasks = parse(Int, ARGS[2])

println("Processing rows: $(collect(task_index:n_tasks:size(args_df, 1)))")

for row_index in task_index:n_tasks:size(args_df, 1)

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
    formulation = String(args_df[row_index, :formulation])
    method = String(args_df[row_index, :method])
    subpath_single_service = args_df[row_index, :subpath_single_service]
    subpath_check_customers = args_df[row_index, :subpath_check_customers]
    path_single_service = args_df[row_index, :path_single_service]
    path_check_customers = args_df[row_index, :path_check_customers]
    check_customers_accelerated = args_df[row_index, :check_customers_accelerated]
    christofides = args_df[row_index, :christofides]
    ngroute = args_df[row_index, :ngroute]
    ngroute_alt = args_df[row_index, :ngroute_alt]
    ngroute_neighborhood_charging_depots_size = String(args_df[row_index, :ngroute_neighborhood_charging_depots_size])

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
        data_dir = "../../../data/",
    )

    if formulation == "subpath"
        (
            r_LP_results, r_IP_results, r_params, r_printlist, r_subpaths, r_charging_arcs
        ) = subpath_formulation_column_generation_integrated_from_paths(
            data, 
            Env = GRB_ENV, 
            method = method, 
            time_windows = use_time_windows,
            subpath_single_service = subpath_single_service,
            subpath_check_customers = subpath_check_customers,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            check_customers_accelerated = check_customers_accelerated,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
            time_limit = TIME_LIMIT,
        )
        try
            collect_subpath_solution_metrics!(r_LP_results, data, r_subpaths, r_charging_arcs)
            collect_subpath_solution_metrics!(r_IP_results, data, r_subpaths, r_charging_arcs)
        catch
            nothing
        end
    elseif formulation == "path"
        (
            r_LP_results, r_IP_results, r_params, r_printlist, r_paths
        ) = path_formulation_column_generation(
            data, 
            Env = GRB_ENV, 
            method = method, 
            time_windows = use_time_windows,
            subpath_single_service = subpath_single_service,
            subpath_check_customers = subpath_check_customers,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            time_limit = TIME_LIMIT,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
        )
        try
            collect_path_solution_metrics!(r_LP_results, data, r_paths)
            collect_path_solution_metrics!(r_IP_results, data, r_paths)
        catch
            nothing
        end
    end
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
            formulation = formulation,
            method = method,
            subpath_single_service = subpath_single_service,
            subpath_check_customers = subpath_check_customers,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            check_customers_accelerated = check_customers_accelerated,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
            LP_objective = r_LP_results["objective"],
            IP_objective = r_IP_results["objective"],
            LP_IP_gap = r_params["LP_IP_gap"],
            counter = r_params["counter"],
            converged = r_params["converged"],
            time_taken = r_params["time_taken"],
            time_limit_reached = r_params["time_limit_reached"],
            sp_base_time_taken_total = r_params["sp_base_time_taken_total"],
            sp_full_time_taken_total = r_params["sp_full_time_taken_total"],
            lp_relaxation_time_taken_total = r_params["lp_relaxation_time_taken_total"],
            sp_base_time_taken_mean = r_params["sp_base_time_taken_mean"],
            sp_full_time_taken_mean = r_params["sp_full_time_taken_mean"],
            lp_relaxation_time_taken_mean = r_params["lp_relaxation_time_taken_mean"],
            lp_mean_subpath_length = get(r_LP_results, "mean_subpath_length", missing),
            lp_weighted_mean_subpath_length = get(r_LP_results, "weighted_mean_subpath_length", missing),
            lp_mean_subpath_ncust = get(r_LP_results, "mean_subpath_ncust", missing),
            lp_weighted_mean_subpath_ncust = get(r_LP_results, "weighted_mean_subpath_ncust", missing),
            lp_mean_path_length = get(r_LP_results, "mean_path_length", missing),
            lp_weighted_mean_path_length = get(r_LP_results, "weighted_mean_path_length", missing),
            lp_mean_path_ncust = get(r_LP_results, "mean_path_ncust", missing),
            lp_weighted_mean_path_ncust = get(r_LP_results, "weighted_mean_path_ncust", missing),
            lp_mean_ps_length = get(r_LP_results, "mean_ps_length", missing),
            lp_weighted_mean_ps_length = get(r_LP_results, "weighted_mean_ps_length", missing),
            lp_utilization = get(r_LP_results, "utilization", missing),
            lp_driving_time_proportion = get(r_LP_results, "driving_time_proportion", missing),
            lp_charging_time_proportion = get(r_LP_results, "charging_time_proportion", missing),
            ip_mean_subpath_length = get(r_IP_results, "mean_subpath_length", missing),
            ip_weighted_mean_subpath_length = get(r_IP_results, "weighted_mean_subpath_length", missing),
            ip_mean_subpath_ncust = get(r_IP_results, "mean_subpath_ncust", missing),
            ip_weighted_mean_subpath_ncust = get(r_IP_results, "weighted_mean_subpath_ncust", missing),
            ip_mean_path_length = get(r_IP_results, "mean_path_length", missing),
            ip_weighted_mean_path_length = get(r_IP_results, "weighted_mean_path_length", missing),
            ip_mean_path_ncust = get(r_IP_results, "mean_path_ncust", missing),
            ip_weighted_mean_path_ncust = get(r_IP_results, "weighted_mean_path_ncust", missing),
            ip_mean_ps_length = get(r_IP_results, "mean_ps_length", missing),
            ip_weighted_mean_ps_length = get(r_IP_results, "weighted_mean_ps_length", missing),
            ip_utilization = get(r_IP_results, "utilization", missing),
            ip_driving_time_proportion = get(r_IP_results, "driving_time_proportion", missing),
            ip_charging_time_proportion = get(r_IP_results, "charging_time_proportion", missing),
        )
    ]
    CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
    println("Row $row_index processed: $(args_df[row_index, :])")
end