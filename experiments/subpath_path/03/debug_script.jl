include("../../../src/subpath_formulation.jl")
include("../../../src/path_formulation.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

using JuMP, Gurobi

const GRB_ENV = Gurobi.Env()

const TIME_LIMIT_SEC = 3600.0

println("Compilation complete.")

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

row_index = 126

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
        T = 60000,
        seed = seed,
        B = 30000,
        μ = μ,
        travel_cost_coeff = travel_cost_coeff,
        charge_cost_coeff = charge_cost_coeff,
        load_scale = load_scale,
        load_shape = load_shape,
        load_tolerance = load_tolerance,
        batch = batch,
        permissiveness = permissiveness,
        # data_dir = "../../../data/",
    )
    plot_instance(data)
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
        time_limit = TIME_LIMIT_SEC,
        christofides = christofides,
        ngroute = ngroute,
        ngroute_alt = ngroute_alt,
        ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
        verbose = true,
    )
            collect_path_solution_metrics!(r_LP_results, data, r_paths)
            collect_path_solution_metrics!(r_IP_results, data, r_paths)

    # if formulation == "subpath"
    #     (
    #         r_LP_results, r_IP_results, r_params, r_printlist, r_subpaths, r_charging_arcs
    #     ) = subpath_formulation_column_generation_integrated_from_paths(
    #         data, 
    #         Env = GRB_ENV, 
    #         method = method, 
    #         time_windows = use_time_windows,
    #         subpath_single_service = subpath_single_service,
    #         subpath_check_customers = subpath_check_customers,
    #         path_single_service = path_single_service,
    #         path_check_customers = path_check_customers,
    #         check_customers_accelerated = check_customers_accelerated,
    #         christofides = christofides,
    #         ngroute = ngroute,
    #         ngroute_alt = ngroute_alt,
    #         ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
    #         time_limit = TIME_LIMIT_SEC,
    #     )
    #     try
    #         collect_subpath_solution_metrics!(r_LP_results, data, r_subpaths, r_charging_arcs)
    #         collect_subpath_solution_metrics!(r_IP_results, data, r_subpaths, r_charging_arcs)
    #     catch
    #         nothing
    #     end
    # elseif formulation == "path"
    #     (
    #         r_LP_results, r_IP_results, r_params, r_printlist, r_paths
    #     ) = path_formulation_column_generation(
    #         data, 
    #         Env = GRB_ENV, 
    #         method = method, 
    #         time_windows = use_time_windows,
    #         subpath_single_service = subpath_single_service,
    #         subpath_check_customers = subpath_check_customers,
    #         path_single_service = path_single_service,
    #         path_check_customers = path_check_customers,
    #         time_limit = TIME_LIMIT_SEC,
    #         christofides = christofides,
    #         ngroute = ngroute,
    #         ngroute_alt = ngroute_alt,
    #         ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
    #     )
    #     try
    #         collect_path_solution_metrics!(r_LP_results, data, r_paths)
    #         collect_path_solution_metrics!(r_IP_results, data, r_paths)
    #     catch
    #         nothing
    #     end
    # end
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