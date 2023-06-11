include("../../../src/subpath_formulation.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

# simple test case to quickly compile 
begin
    data = generate_instance(
        n_depots = 4,
        n_customers = 10,
        n_charging = 9,
        n_vehicles = 6,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "grid",
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
    )
    G = construct_graph(data)
    (
        b_LP_results, b_IP_results, b_params, b_printlist, b_subpaths, b_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "benchmark",
    )
    (
        b_s_LP_results, b_s_IP_results, b_s_params, b_s_printlist, b_s_subpaths, b_s_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "benchmark", 
        path_single_service = true,
    )
    (
        b_sc_LP_results, b_sc_IP_results, b_sc_params, b_sc_printlist, b_sc_subpaths, b_sc_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "benchmark", 
        path_single_service = true, path_check_customers = true,
    )
    (
        b_sca_LP_results, b_sca_IP_results, b_sca_params, b_sca_printlist, b_sca_subpaths, b_sca_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "benchmark", 
        path_single_service = true, path_check_customers = true, check_customers_accelerated = true
    )
    (
        o_LP_results, o_IP_results, o_params, o_printlist, o_subpaths, o_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours",
    )
    (
        o_s_LP_results, o_s_IP_results, o_s_params, o_s_printlist, o_s_subpaths, o_s_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours", 
        subpath_single_service = true,
    )
    (
        o_sc_LP_results, o_sc_IP_results, o_sc_params, o_sc_printlist, o_sc_subpaths, o_sc_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours", 
        subpath_single_service = true, subpath_check_customers = true,
    )
    (
        o_sca_LP_results, o_sca_IP_results, o_sca_params, o_sca_printlist, o_sca_subpaths, o_sca_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours", 
        subpath_single_service = true, subpath_check_customers = true, 
        check_customers_accelerated = true,
    )
    (
        o_ss_LP_results, o_ss_IP_results, o_ss_params, o_ss_printlist, o_ss_subpaths, o_ss_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours", 
        subpath_single_service = true,
        path_single_service = true,
    )
    (
        o_scsc_LP_results, o_scsc_IP_results, o_scsc_params, o_scsc_printlist, o_scsc_subpaths, o_scsc_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours", 
        subpath_single_service = true, subpath_check_customers = true,
        path_single_service = true, path_check_customers = true,
    )
    (
        o_scsca_LP_results, o_scsca_IP_results, o_scsca_params, o_scsca_printlist, o_scsca_subpaths, o_scsca_charging_arcs
    ) = subpath_formulation_column_generation_integrated_from_paths(
        G, data, method = "ours", 
        subpath_single_service = true, subpath_check_customers = true, 
        path_single_service = true, path_check_customers = true,
        check_customers_accelerated = true,
    )
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
    collect_subpath_solution_metrics!(r_LP_results, data, r_subpaths, r_charging_arcs)
    collect_subpath_solution_metrics!(r_IP_results, data, r_subpaths, r_charging_arcs)
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
            lp_mean_subpath_length = r_LP_results["mean_subpath_length"],
            lp_weighted_mean_subpath_length = r_LP_results["weighted_mean_subpath_length"],
            lp_mean_path_length = r_LP_results["mean_path_length"],
            lp_weighted_mean_path_length = r_LP_results["weighted_mean_path_length"],
            lp_mean_ps_length = r_LP_results["mean_ps_length"],
            lp_weighted_mean_ps_length = r_LP_results["weighted_mean_ps_length"],
            ip_mean_subpath_length = r_IP_results["mean_subpath_length"],
            ip_weighted_mean_subpath_length = r_IP_results["weighted_mean_subpath_length"],
            ip_mean_path_length = r_IP_results["mean_path_length"],
            ip_weighted_mean_path_length = r_IP_results["weighted_mean_path_length"],
            ip_mean_ps_length = r_IP_results["mean_ps_length"],
            ip_weighted_mean_ps_length = r_IP_results["weighted_mean_ps_length"],
        )
    ]
    CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
    println("Row $row_index processed: $(args_df[row_index, :])")
end