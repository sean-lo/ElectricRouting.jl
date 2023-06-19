include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames

using Test

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 16,
    n_charging = 9,
    n_vehicles = 7,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "grid",
    shrinkage_depots = 1.0,
    shrinkage_charging = 0.7,
    T = 40000,
    seed = 5,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 4,
    permissiveness = 0.2,
)
G = construct_graph(data)
# plot_instance(data)

p_b_tw_LP_results, p_b_tw_IP_results, p_b_tw_params, p_b_tw_printlist, p_b_tw_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true)
p_b_tw_s_LP_results, p_b_tw_s_IP_results, p_b_tw_s_params, p_b_tw_s_printlist, p_b_tw_s_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, path_single_service = true, verbose = true)
p_b_tw_sc_LP_results, p_b_tw_sc_IP_results, p_b_tw_sc_params, p_b_tw_sc_printlist, p_b_tw_sc_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, path_single_service = true, path_check_customers = true, verbose = true)

p_b_tw_ch_LP_results, p_b_tw_ch_IP_results, p_b_tw_ch_params, p_b_tw_ch_printlist, p_b_tw_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, christofides = true, verbose = true)

p_b_tw_ie_LP_results, p_b_tw_ie_IP_results, p_b_tw_ie_params, p_b_tw_ie_printlist, p_b_tw_ie_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true,)
p_b_tw_iec_LP_results, p_b_tw_iec_IP_results, p_b_tw_iec_params, p_b_tw_iec_printlist, p_b_tw_iec_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true)
p_b_tw_iecw_LP_results, p_b_tw_iecw_IP_results, p_b_tw_iecw_params, p_b_tw_iecw_printlist, p_b_tw_iecw_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true)
p_b_tw_ie_hmaa_LP_results, p_b_tw_ie_hmaa_IP_results, p_b_tw_ie_hmaa_params, p_b_tw_ie_hmaa_printlist, p_b_tw_ie_hmaa_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, incremental_elementarity_rule = "hmaa")
p_b_tw_iec_hmaa_LP_results, p_b_tw_iec_hmaa_IP_results, p_b_tw_iec_hmaa_params, p_b_tw_iec_hmaa_printlist, p_b_tw_iec_hmaa_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa")
p_b_tw_iecw_hmaa_LP_results, p_b_tw_iecw_hmaa_IP_results, p_b_tw_iecw_hmaa_params, p_b_tw_iecw_hmaa_printlist, p_b_tw_iecw_hmaa_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa")

p_b_tw_ie_ch_LP_results, p_b_tw_ie_ch_IP_results, p_b_tw_ie_ch_params, p_b_tw_ie_ch_printlist, p_b_tw_ie_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, christofides = true)
p_b_tw_iec_ch_LP_results, p_b_tw_iec_ch_IP_results, p_b_tw_iec_ch_params, p_b_tw_iec_ch_printlist, p_b_tw_iec_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, christofides = true)
p_b_tw_iecw_ch_LP_results, p_b_tw_iecw_ch_IP_results, p_b_tw_iecw_ch_params, p_b_tw_iecw_ch_printlist, p_b_tw_iecw_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, christofides = true)
p_b_tw_ie_hmaa_ch_LP_results, p_b_tw_ie_hmaa_ch_IP_results, p_b_tw_ie_hmaa_ch_params, p_b_tw_ie_hmaa_ch_printlist, p_b_tw_ie_hmaa_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, incremental_elementarity_rule = "hmaa", christofides = true)
p_b_tw_iec_hmaa_ch_LP_results, p_b_tw_iec_hmaa_ch_IP_results, p_b_tw_iec_hmaa_ch_params, p_b_tw_iec_hmaa_ch_printlist, p_b_tw_iec_hmaa_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa", christofides = true)
p_b_tw_iecw_hmaa_ch_LP_results, p_b_tw_iecw_hmaa_ch_IP_results, p_b_tw_iecw_hmaa_ch_params, p_b_tw_iecw_hmaa_ch_printlist, p_b_tw_iecw_hmaa_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa", christofides = true)

p_b_LP_results, p_b_IP_results, p_b_params, p_b_printlist, p_b_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true)
p_b_s_LP_results, p_b_s_IP_results, p_b_s_params, p_b_s_printlist, p_b_s_some_paths = path_formulation_column_generation(G, data; method = "benchmark", path_single_service = true, verbose = true)
p_b_sc_LP_results, p_b_sc_IP_results, p_b_sc_params, p_b_sc_printlist, p_b_sc_some_paths = path_formulation_column_generation(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)

p_b_ch_LP_results, p_b_ch_IP_results, p_b_ch_params, p_b_ch_printlist, p_b_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, christofides = true)

p_b_ie_LP_results, p_b_ie_IP_results, p_b_ie_params, p_b_ie_printlist, p_b_ie_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true,)
p_b_iec_LP_results, p_b_iec_IP_results, p_b_iec_params, p_b_iec_printlist, p_b_iec_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true)
p_b_iecw_LP_results, p_b_iecw_IP_results, p_b_iecw_params, p_b_iecw_printlist, p_b_iecw_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true) # FIXME: debug this
p_b_ie_hmaa_LP_results, p_b_ie_hmaa_IP_results, p_b_ie_hmaa_params, p_b_ie_hmaa_printlist, p_b_ie_hmaa_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, incremental_elementarity_rule = "hmaa")
p_b_iec_hmaa_LP_results, p_b_iec_hmaa_IP_results, p_b_iec_hmaa_params, p_b_iec_hmaa_printlist, p_b_iec_hmaa_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa")
p_b_iecw_hmaa_LP_results, p_b_iecw_hmaa_IP_results, p_b_iecw_hmaa_params, p_b_iecw_hmaa_printlist, p_b_iecw_hmaa_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa")

p_b_ie_ch_LP_results, p_b_ie_ch_IP_results, p_b_ie_ch_params, p_b_ie_ch_printlist, p_b_ie_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, christofides = true)
p_b_iec_ch_LP_results, p_b_iec_ch_IP_results, p_b_iec_ch_params, p_b_iec_ch_printlist, p_b_iec_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, christofides = true)
p_b_iecw_ch_LP_results, p_b_iecw_ch_IP_results, p_b_iecw_ch_params, p_b_iecw_ch_printlist, p_b_iecw_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, christofides = true)
p_b_ie_hmaa_ch_LP_results, p_b_ie_hmaa_ch_IP_results, p_b_ie_hmaa_ch_params, p_b_ie_hmaa_ch_printlist, p_b_ie_hmaa_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, incremental_elementarity_rule = "hmaa", christofides = true)
p_b_iec_hmaa_ch_LP_results, p_b_iec_hmaa_ch_IP_results, p_b_iec_hmaa_ch_params, p_b_iec_hmaa_ch_printlist, p_b_iec_hmaa_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa", christofides = true)
p_b_iecw_hmaa_ch_LP_results, p_b_iecw_hmaa_ch_IP_results, p_b_iecw_hmaa_ch_params, p_b_iecw_hmaa_ch_printlist, p_b_iecw_hmaa_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa", christofides = true)

p_o_LP_results, p_o_IP_results, p_o_params, p_o_printlist, p_o_some_paths = path_formulation_column_generation(G, data; method = "ours", verbose = true)
p_o_s_LP_results, p_o_s_IP_results, p_o_s_params, p_o_s_printlist, p_o_s_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, verbose = true)
p_o_sc_LP_results, p_o_sc_IP_results, p_o_sc_params, p_o_sc_printlist, p_o_sc_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
p_o_ss_LP_results, p_o_ss_IP_results, p_o_ss_params, p_o_ss_printlist, p_o_ss_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
p_o_scsc_LP_results, p_o_scsc_IP_results, p_o_scsc_params, p_o_scsc_printlist, p_o_scsc_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)

p_o_ch_LP_results, p_o_ch_IP_results, p_o_ch_params, p_o_ch_printlist, p_o_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", verbose = true, christofides = true)
p_o_s_ch_LP_results, p_o_s_ch_IP_results, p_o_s_ch_params, p_o_s_ch_printlist, p_o_s_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, verbose = true, christofides = true)
p_o_sc_ch_LP_results, p_o_sc_ch_IP_results, p_o_sc_ch_params, p_o_sc_ch_printlist, p_o_sc_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true, christofides = true)
p_o_ss_ch_LP_results, p_o_ss_ch_IP_results, p_o_ss_ch_params, p_o_ss_ch_printlist, p_o_ss_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true, christofides = true)
p_o_scsc_ch_LP_results, p_o_scsc_ch_IP_results, p_o_scsc_ch_params, p_o_scsc_ch_printlist, p_o_scsc_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true, christofides = true)

sp_b_tw_LP_results, sp_b_tw_IP_results, sp_b_tw_params, sp_b_tw_printlist, sp_b_tw_some_subpaths, sp_b_tw_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, verbose = true)
sp_b_tw_s_LP_results, sp_b_tw_s_IP_results, sp_b_tw_s_params, sp_b_tw_s_printlist, sp_b_tw_s_some_subpaths, sp_b_tw_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, path_single_service = true, verbose = true)
sp_b_tw_sc_LP_results, sp_b_tw_sc_IP_results, sp_b_tw_sc_params, sp_b_tw_sc_printlist, sp_b_tw_sc_some_subpaths, sp_b_tw_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, path_single_service = true, path_check_customers = true, verbose = true)
sp_b_tw_sca_LP_results, sp_b_tw_sca_IP_results, sp_b_tw_sca_params, sp_b_tw_sca_printlist, sp_b_tw_sca_some_subpaths, sp_b_tw_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

sp_b_tw_ch_LP_results, sp_b_tw_ch_IP_results, sp_b_tw_ch_params, sp_b_tw_ch_printlist, sp_b_tw_ch_some_subpaths, sp_b_tw_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, verbose = true, christofides = true)

sp_b_tw_ie_LP_results, sp_b_tw_ie_IP_results, sp_b_tw_ie_params, sp_b_tw_ie_printlist, sp_b_tw_ie_some_subpaths, sp_b_tw_ie_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, verbose = true)
sp_b_tw_iec_LP_results, sp_b_tw_iec_IP_results, sp_b_tw_iec_params, sp_b_tw_iec_printlist, sp_b_tw_iec_some_subpaths, sp_b_tw_iec_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, verbose = true)
sp_b_tw_iecw_LP_results, sp_b_tw_iecw_IP_results, sp_b_tw_iecw_params, sp_b_tw_iecw_printlist, sp_b_tw_iecw_some_subpaths, sp_b_tw_iecw_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, verbose = true)
sp_b_tw_ie_hmaa_LP_results, sp_b_tw_ie_hmaa_IP_results, sp_b_tw_ie_hmaa_params, sp_b_tw_ie_hmaa_printlist, sp_b_tw_ie_hmaa_some_subpaths, sp_b_tw_ie_hmaa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, incremental_elementarity_rule = "hmaa", verbose = true)
sp_b_tw_iec_hmaa_LP_results, sp_b_tw_iec_hmaa_IP_results, sp_b_tw_iec_hmaa_params, sp_b_tw_iec_hmaa_printlist, sp_b_tw_iec_hmaa_some_subpaths, sp_b_tw_iec_hmaa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa", verbose = true)
sp_b_tw_iecw_hmaa_LP_results, sp_b_tw_iecw_hmaa_IP_results, sp_b_tw_iecw_hmaa_params, sp_b_tw_iecw_hmaa_printlist, sp_b_tw_iecw_hmaa_some_subpaths, sp_b_tw_iecw_hmaa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa", verbose = true)


sp_b_tw_ie_ch_LP_results, sp_b_tw_ie_ch_IP_results, sp_b_tw_ie_ch_params, sp_b_tw_ie_ch_printlist, sp_b_tw_ie_ch_some_subpaths, sp_b_tw_ie_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, verbose = true, christofides = true)
sp_b_tw_iec_ch_LP_results, sp_b_tw_iec_ch_IP_results, sp_b_tw_iec_ch_params, sp_b_tw_iec_ch_printlist, sp_b_tw_iec_ch_some_subpaths, sp_b_tw_iec_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, verbose = true, christofides = true)
sp_b_tw_iecw_ch_LP_results, sp_b_tw_iecw_ch_IP_results, sp_b_tw_iecw_ch_params, sp_b_tw_iecw_ch_printlist, sp_b_tw_iecw_ch_some_subpaths, sp_b_tw_iecw_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, verbose = true, christofides = true)
sp_b_tw_ie_hmaa_ch_LP_results, sp_b_tw_ie_hmaa_ch_IP_results, sp_b_tw_ie_hmaa_ch_params, sp_b_tw_ie_hmaa_ch_printlist, sp_b_tw_ie_hmaa_ch_some_subpaths, sp_b_tw_ie_hmaa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, incremental_elementarity_rule = "hmaa", verbose = true, christofides = true)
sp_b_tw_iec_hmaa_ch_LP_results, sp_b_tw_iec_hmaa_ch_IP_results, sp_b_tw_iec_hmaa_ch_params, sp_b_tw_iec_hmaa_ch_printlist, sp_b_tw_iec_hmaa_ch_some_subpaths, sp_b_tw_iec_hmaa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa", verbose = true, christofides = true)
sp_b_tw_iecw_hmaa_ch_LP_results, sp_b_tw_iecw_hmaa_ch_IP_results, sp_b_tw_iecw_hmaa_ch_params, sp_b_tw_iecw_hmaa_ch_printlist, sp_b_tw_iecw_hmaa_ch_some_subpaths, sp_b_tw_iecw_hmaa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa", verbose = true, christofides = true)

sp_b_LP_results, sp_b_IP_results, sp_b_params, sp_b_printlist, sp_b_some_subpaths, sp_b_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", verbose = true)
sp_b_s_LP_results, sp_b_s_IP_results, sp_b_s_params, sp_b_s_printlist, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, verbose = true)
sp_b_sc_LP_results, sp_b_sc_IP_results, sp_b_sc_params, sp_b_sc_printlist, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)
sp_b_sca_LP_results, sp_b_sca_IP_results, sp_b_sca_params, sp_b_sca_printlist, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

sp_b_ch_LP_results, sp_b_ch_IP_results, sp_b_ch_params, sp_b_ch_printlist, sp_b_ch_some_subpaths, sp_b_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", verbose = true, christofides = true)

sp_b_ie_LP_results, sp_b_ie_IP_results, sp_b_ie_params, sp_b_ie_printlist, sp_b_ie_some_subpaths, sp_b_ie_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, verbose = true)
sp_b_iec_LP_results, sp_b_iec_IP_results, sp_b_iec_params, sp_b_iec_printlist, sp_b_iec_some_subpaths, sp_b_iec_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, verbose = true)
sp_b_iecw_LP_results, sp_b_iecw_IP_results, sp_b_iecw_params, sp_b_iecw_printlist, sp_b_iecw_some_subpaths, sp_b_iecw_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, warm_start = true, verbose = true)
sp_b_ie_hmaa_LP_results, sp_b_ie_hmaa_IP_results, sp_b_ie_hmaa_params, sp_b_ie_hmaa_printlist, sp_b_ie_hmaa_some_subpaths, sp_b_ie_hmaa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, incremental_elementarity_rule = "hmaa", verbose = true)
sp_b_iec_hmaa_LP_results, sp_b_iec_hmaa_IP_results, sp_b_iec_hmaa_params, sp_b_iec_hmaa_printlist, sp_b_iec_hmaa_some_subpaths, sp_b_iec_hmaa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa", verbose = true)
sp_b_iecw_hmaa_LP_results, sp_b_iecw_hmaa_IP_results, sp_b_iecw_hmaa_params, sp_b_iecw_hmaa_printlist, sp_b_iecw_hmaa_some_subpaths, sp_b_iecw_hmaa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa", verbose = true)

sp_b_ie_ch_LP_results, sp_b_ie_ch_IP_results, sp_b_ie_ch_params, sp_b_ie_ch_printlist, sp_b_ie_ch_some_subpaths, sp_b_ie_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, verbose = true, christofides = true)
sp_b_iec_ch_LP_results, sp_b_iec_ch_IP_results, sp_b_iec_ch_params, sp_b_iec_ch_printlist, sp_b_iec_ch_some_subpaths, sp_b_iec_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, verbose = true, christofides = true)
sp_b_iecw_ch_LP_results, sp_b_iecw_ch_IP_results, sp_b_iecw_ch_params, sp_b_iecw_ch_printlist, sp_b_iecw_ch_some_subpaths, sp_b_iecw_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, warm_start = true, verbose = true, christofides = true)
sp_b_ie_hmaa_ch_LP_results, sp_b_ie_hmaa_ch_IP_results, sp_b_ie_hmaa_ch_params, sp_b_ie_hmaa_ch_printlist, sp_b_ie_hmaa_ch_some_subpaths, sp_b_ie_hmaa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, incremental_elementarity_rule = "hmaa", verbose = true, christofides = true)
sp_b_iec_hmaa_ch_LP_results, sp_b_iec_hmaa_ch_IP_results, sp_b_iec_hmaa_ch_params, sp_b_iec_hmaa_ch_printlist, sp_b_iec_hmaa_ch_some_subpaths, sp_b_iec_hmaa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, incremental_elementarity_rule = "hmaa", verbose = true, christofides = true)
sp_b_iecw_hmaa_ch_LP_results, sp_b_iecw_hmaa_ch_IP_results, sp_b_iecw_hmaa_ch_params, sp_b_iecw_hmaa_ch_printlist, sp_b_iecw_hmaa_ch_some_subpaths, sp_b_iecw_hmaa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", incremental_elementarity = true, path_check_customers = true, warm_start = true, incremental_elementarity_rule = "hmaa", verbose = true, christofides = true)


sp_o_LP_results, sp_o_IP_results, sp_o_params, sp_o_printlist, sp_o_some_subpaths, sp_o_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", verbose = true)
sp_o_s_LP_results, sp_o_s_IP_results, sp_o_s_params, sp_o_s_printlist, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, verbose = true)
sp_o_sc_LP_results, sp_o_sc_IP_results, sp_o_sc_params, sp_o_sc_printlist, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
sp_o_sca_LP_results, sp_o_sca_IP_results, sp_o_sca_params, sp_o_sca_printlist, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true)
sp_o_ss_LP_results, sp_o_ss_IP_results, sp_o_ss_params, sp_o_ss_printlist, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
sp_o_scsc_LP_results, sp_o_scsc_IP_results, sp_o_scsc_params, sp_o_scsc_printlist, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)
sp_o_scsca_LP_results, sp_o_scsca_IP_results, sp_o_scsca_params, sp_o_scsca_printlist, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

sp_o_ch_LP_results, sp_o_ch_IP_results, sp_o_ch_params, sp_o_ch_printlist, sp_o_ch_some_subpaths, sp_o_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", verbose = true, christofides = true)
sp_o_s_ch_LP_results, sp_o_s_ch_IP_results, sp_o_s_ch_params, sp_o_s_ch_printlist, sp_o_s_ch_some_subpaths, sp_o_s_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, verbose = true, christofides = true)
sp_o_sc_ch_LP_results, sp_o_sc_ch_IP_results, sp_o_sc_ch_params, sp_o_sc_ch_printlist, sp_o_sc_ch_some_subpaths, sp_o_sc_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true, christofides = true)
sp_o_sca_ch_LP_results, sp_o_sca_ch_IP_results, sp_o_sca_ch_params, sp_o_sca_ch_printlist, sp_o_sca_ch_some_subpaths, sp_o_sca_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true, christofides = true)
sp_o_ss_ch_LP_results, sp_o_ss_ch_IP_results, sp_o_ss_ch_params, sp_o_ss_ch_printlist, sp_o_ss_ch_some_subpaths, sp_o_ss_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true, christofides = true)
sp_o_scsc_ch_LP_results, sp_o_scsc_ch_IP_results, sp_o_scsc_ch_params, sp_o_scsc_ch_printlist, sp_o_scsc_ch_some_subpaths, sp_o_scsc_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true, christofides = true)
sp_o_scsca_ch_LP_results, sp_o_scsca_ch_IP_results, sp_o_scsca_ch_params, sp_o_scsca_ch_printlist, sp_o_scsca_ch_some_subpaths, sp_o_scsca_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true, christofides = true)

collect_path_solution_metrics!(p_b_tw_LP_results, data, p_b_tw_some_paths)
collect_path_solution_metrics!(p_b_tw_s_LP_results, data, p_b_tw_s_some_paths)
collect_path_solution_metrics!(p_b_tw_sc_LP_results, data, p_b_tw_sc_some_paths)
collect_path_solution_metrics!(p_b_tw_ch_LP_results, data, p_b_tw_ch_some_paths)

collect_path_solution_metrics!(p_b_tw_ie_LP_results, data, p_b_tw_ie_some_paths)
collect_path_solution_metrics!(p_b_tw_iec_LP_results, data, p_b_tw_iec_some_paths)
collect_path_solution_metrics!(p_b_tw_iecw_LP_results, data, p_b_tw_iecw_some_paths)
collect_path_solution_metrics!(p_b_tw_ie_hmaa_LP_results, data, p_b_tw_ie_hmaa_some_paths)
collect_path_solution_metrics!(p_b_tw_iec_hmaa_LP_results, data, p_b_tw_iec_hmaa_some_paths)
collect_path_solution_metrics!(p_b_tw_iecw_hmaa_LP_results, data, p_b_tw_iecw_hmaa_some_paths)
collect_path_solution_metrics!(p_b_tw_ie_ch_LP_results, data, p_b_tw_ie_ch_some_paths)
collect_path_solution_metrics!(p_b_tw_iec_ch_LP_results, data, p_b_tw_iec_ch_some_paths)
collect_path_solution_metrics!(p_b_tw_iecw_ch_LP_results, data, p_b_tw_iecw_ch_some_paths)
collect_path_solution_metrics!(p_b_tw_ie_hmaa_ch_LP_results, data, p_b_tw_ie_hmaa_ch_some_paths)
collect_path_solution_metrics!(p_b_tw_iec_hmaa_ch_LP_results, data, p_b_tw_iec_hmaa_ch_some_paths)
collect_path_solution_metrics!(p_b_tw_iecw_hmaa_ch_LP_results, data, p_b_tw_iecw_hmaa_ch_some_paths)

collect_path_solution_metrics!(p_b_LP_results, data, p_b_some_paths)
collect_path_solution_metrics!(p_b_s_LP_results, data, p_b_s_some_paths)
collect_path_solution_metrics!(p_b_sc_LP_results, data, p_b_sc_some_paths)
collect_path_solution_metrics!(p_b_ch_LP_results, data, p_b_ch_some_paths)

collect_path_solution_metrics!(p_b_ie_LP_results, data, p_b_ie_some_paths)
collect_path_solution_metrics!(p_b_iec_LP_results, data, p_b_iec_some_paths)
collect_path_solution_metrics!(p_b_iecw_LP_results, data, p_b_iecw_some_paths)
collect_path_solution_metrics!(p_b_ie_hmaa_LP_results, data, p_b_ie_hmaa_some_paths)
collect_path_solution_metrics!(p_b_iec_hmaa_LP_results, data, p_b_iec_hmaa_some_paths)
collect_path_solution_metrics!(p_b_iecw_hmaa_LP_results, data, p_b_iecw_hmaa_some_paths)
collect_path_solution_metrics!(p_b_ie_ch_LP_results, data, p_b_ie_ch_some_paths)
collect_path_solution_metrics!(p_b_iec_ch_LP_results, data, p_b_iec_ch_some_paths)
collect_path_solution_metrics!(p_b_iecw_ch_LP_results, data, p_b_iecw_ch_some_paths)
collect_path_solution_metrics!(p_b_ie_hmaa_ch_LP_results, data, p_b_ie_hmaa_ch_some_paths)
collect_path_solution_metrics!(p_b_iec_hmaa_ch_LP_results, data, p_b_iec_hmaa_ch_some_paths)
collect_path_solution_metrics!(p_b_iecw_hmaa_ch_LP_results, data, p_b_iecw_hmaa_ch_some_paths)

collect_path_solution_metrics!(p_o_LP_results, data, p_o_some_paths)
collect_path_solution_metrics!(p_o_s_LP_results, data, p_o_s_some_paths)
collect_path_solution_metrics!(p_o_sc_LP_results, data, p_o_sc_some_paths)
collect_path_solution_metrics!(p_o_ss_LP_results, data, p_o_ss_some_paths)
collect_path_solution_metrics!(p_o_scsc_LP_results, data, p_o_scsc_some_paths)

collect_path_solution_metrics!(p_o_ch_LP_results, data, p_o_ch_some_paths)
collect_path_solution_metrics!(p_o_s_ch_LP_results, data, p_o_s_ch_some_paths)
collect_path_solution_metrics!(p_o_sc_ch_LP_results, data, p_o_sc_ch_some_paths)
collect_path_solution_metrics!(p_o_ss_ch_LP_results, data, p_o_ss_ch_some_paths)
collect_path_solution_metrics!(p_o_scsc_ch_LP_results, data, p_o_scsc_ch_some_paths)

collect_subpath_solution_metrics!(sp_b_tw_LP_results, data, sp_b_tw_some_subpaths, sp_b_tw_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_s_LP_results, data, sp_b_tw_s_some_subpaths, sp_b_tw_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_sc_LP_results, data, sp_b_tw_sc_some_subpaths, sp_b_tw_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_sca_LP_results, data, sp_b_tw_sca_some_subpaths, sp_b_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_ch_LP_results, data, sp_b_tw_ch_some_subpaths, sp_b_tw_ch_some_charging_arcs)

collect_subpath_solution_metrics!(sp_b_tw_ie_LP_results, data, sp_b_tw_ie_some_subpaths, sp_b_tw_ie_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iec_LP_results, data, sp_b_tw_iec_some_subpaths, sp_b_tw_iec_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iecw_LP_results, data, sp_b_tw_iecw_some_subpaths, sp_b_tw_iecw_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_ie_hmaa_LP_results, data, sp_b_tw_ie_hmaa_some_subpaths, sp_b_tw_ie_hmaa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iec_hmaa_LP_results, data, sp_b_tw_iec_hmaa_some_subpaths, sp_b_tw_iec_hmaa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iecw_hmaa_LP_results, data, sp_b_tw_iecw_hmaa_some_subpaths, sp_b_tw_iecw_hmaa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_ie_ch_LP_results, data, sp_b_tw_ie_ch_some_subpaths, sp_b_tw_ie_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iec_ch_LP_results, data, sp_b_tw_iec_ch_some_subpaths, sp_b_tw_iec_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iecw_ch_LP_results, data, sp_b_tw_iecw_ch_some_subpaths, sp_b_tw_iecw_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_ie_hmaa_ch_LP_results, data, sp_b_tw_ie_hmaa_ch_some_subpaths, sp_b_tw_ie_hmaa_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iec_hmaa_ch_LP_results, data, sp_b_tw_iec_hmaa_ch_some_subpaths, sp_b_tw_iec_hmaa_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_iecw_hmaa_ch_LP_results, data, sp_b_tw_iecw_hmaa_ch_some_subpaths, sp_b_tw_iecw_hmaa_ch_some_charging_arcs)

collect_subpath_solution_metrics!(sp_b_LP_results, data, sp_b_some_subpaths, sp_b_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_s_LP_results, data, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sc_LP_results, data, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sca_LP_results, data, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ch_LP_results, data, sp_b_ch_some_subpaths, sp_b_ch_some_charging_arcs)

collect_subpath_solution_metrics!(sp_b_ie_LP_results, data, sp_b_ie_some_subpaths, sp_b_ie_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iec_LP_results, data, sp_b_iec_some_subpaths, sp_b_iec_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iecw_LP_results, data, sp_b_iecw_some_subpaths, sp_b_iecw_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ie_hmaa_LP_results, data, sp_b_ie_hmaa_some_subpaths, sp_b_ie_hmaa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iec_hmaa_LP_results, data, sp_b_iec_hmaa_some_subpaths, sp_b_iec_hmaa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iecw_hmaa_LP_results, data, sp_b_iecw_hmaa_some_subpaths, sp_b_iecw_hmaa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ie_ch_LP_results, data, sp_b_ie_ch_some_subpaths, sp_b_ie_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iec_ch_LP_results, data, sp_b_iec_ch_some_subpaths, sp_b_iec_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iecw_ch_LP_results, data, sp_b_iecw_ch_some_subpaths, sp_b_iecw_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ie_hmaa_ch_LP_results, data, sp_b_ie_hmaa_ch_some_subpaths, sp_b_ie_hmaa_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iec_hmaa_ch_LP_results, data, sp_b_iec_hmaa_ch_some_subpaths, sp_b_iec_hmaa_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_iecw_hmaa_ch_LP_results, data, sp_b_iecw_hmaa_ch_some_subpaths, sp_b_iecw_hmaa_ch_some_charging_arcs)

collect_subpath_solution_metrics!(sp_o_LP_results, data, sp_o_some_subpaths, sp_o_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_s_LP_results, data, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sc_LP_results, data, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sca_LP_results, data, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ss_LP_results, data, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsc_LP_results, data, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsca_LP_results, data, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs)

collect_subpath_solution_metrics!(sp_o_ch_LP_results, data, sp_o_ch_some_subpaths, sp_o_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_s_ch_LP_results, data, sp_o_s_ch_some_subpaths, sp_o_s_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sc_ch_LP_results, data, sp_o_sc_ch_some_subpaths, sp_o_sc_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sca_ch_LP_results, data, sp_o_sca_ch_some_subpaths, sp_o_sca_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ss_ch_LP_results, data, sp_o_ss_ch_some_subpaths, sp_o_ss_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsc_ch_LP_results, data, sp_o_scsc_ch_some_subpaths, sp_o_scsc_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsca_ch_LP_results, data, sp_o_scsca_ch_some_subpaths, sp_o_scsca_ch_some_charging_arcs)

p_b_tw_LP_results["objective"]
p_b_tw_s_LP_results["objective"]
p_b_tw_sc_LP_results["objective"]
p_b_tw_ch_LP_results["objective"]
p_b_tw_ie_LP_results["objective"]
p_b_tw_iec_LP_results["objective"]
p_b_tw_iecw_LP_results["objective"]
p_b_tw_ie_hmaa_LP_results["objective"]
p_b_tw_iec_hmaa_LP_results["objective"]
p_b_tw_iecw_hmaa_LP_results["objective"]
p_b_tw_ie_ch_LP_results["objective"]
p_b_tw_iec_ch_LP_results["objective"]
p_b_tw_iecw_ch_LP_results["objective"]
p_b_tw_ie_hmaa_ch_LP_results["objective"]
p_b_tw_iec_hmaa_ch_LP_results["objective"]
p_b_tw_iecw_hmaa_ch_LP_results["objective"]

p_b_LP_results["objective"]
p_b_s_LP_results["objective"]
p_b_sc_LP_results["objective"]
p_b_ch_LP_results["objective"]
p_b_ie_LP_results["objective"]
p_b_iec_LP_results["objective"]
p_b_iecw_LP_results["objective"]
p_b_ie_hmaa_LP_results["objective"]
p_b_iec_hmaa_LP_results["objective"]
p_b_iecw_hmaa_LP_results["objective"]
p_b_ie_ch_LP_results["objective"]
p_b_iec_ch_LP_results["objective"]
p_b_iecw_ch_LP_results["objective"]
p_b_ie_hmaa_ch_LP_results["objective"]
p_b_iec_hmaa_ch_LP_results["objective"]
p_b_iecw_hmaa_ch_LP_results["objective"]

p_o_LP_results["objective"]
p_o_s_LP_results["objective"]
p_o_sc_LP_results["objective"]
p_o_ss_LP_results["objective"]
p_o_scsc_LP_results["objective"]
p_o_ch_LP_results["objective"]
p_o_s_ch_LP_results["objective"]
p_o_sc_ch_LP_results["objective"]
p_o_ss_ch_LP_results["objective"]
p_o_scsc_ch_LP_results["objective"]

sp_b_tw_LP_results["objective"]
sp_b_tw_s_LP_results["objective"]
sp_b_tw_sc_LP_results["objective"]
sp_b_tw_sca_LP_results["objective"]
sp_b_tw_ch_LP_results["objective"]
sp_b_tw_ie_LP_results["objective"]
sp_b_tw_iec_LP_results["objective"]
sp_b_tw_iecw_LP_results["objective"]
sp_b_tw_ie_hmaa_LP_results["objective"]
sp_b_tw_iec_hmaa_LP_results["objective"]
sp_b_tw_iecw_hmaa_LP_results["objective"]
sp_b_tw_ie_ch_LP_results["objective"]
sp_b_tw_iec_ch_LP_results["objective"]
sp_b_tw_iecw_ch_LP_results["objective"]
sp_b_tw_ie_hmaa_ch_LP_results["objective"]
sp_b_tw_iec_hmaa_ch_LP_results["objective"]
sp_b_tw_iecw_hmaa_ch_LP_results["objective"]

sp_b_LP_results["objective"]
sp_b_s_LP_results["objective"]
sp_b_sc_LP_results["objective"]
sp_b_sca_LP_results["objective"]
sp_b_ch_LP_results["objective"]
sp_b_ie_LP_results["objective"]
sp_b_iec_LP_results["objective"]
sp_b_iecw_LP_results["objective"]
sp_b_ie_hmaa_LP_results["objective"]
sp_b_iec_hmaa_LP_results["objective"]
sp_b_iecw_hmaa_LP_results["objective"]
sp_b_ie_ch_LP_results["objective"]
sp_b_iec_ch_LP_results["objective"]
sp_b_iecw_ch_LP_results["objective"]
sp_b_ie_hmaa_ch_LP_results["objective"]
sp_b_iec_hmaa_ch_LP_results["objective"]
sp_b_iecw_hmaa_ch_LP_results["objective"]

sp_o_LP_results["objective"]
sp_o_s_LP_results["objective"]
sp_o_sc_LP_results["objective"]
sp_o_sca_LP_results["objective"]
sp_o_ss_LP_results["objective"]
sp_o_scsc_LP_results["objective"]
sp_o_scsca_LP_results["objective"]
sp_o_ch_LP_results["objective"]
sp_o_s_ch_LP_results["objective"]
sp_o_sc_ch_LP_results["objective"]
sp_o_sca_ch_LP_results["objective"]
sp_o_ss_ch_LP_results["objective"]
sp_o_scsc_ch_LP_results["objective"]
sp_o_scsca_ch_LP_results["objective"]

@test (
    p_b_tw_LP_results["objective"] 
    ≤ p_b_tw_ch_LP_results["objective"]
    ≤ p_b_tw_sc_LP_results["objective"]
)
@test p_b_tw_s_LP_results["objective"] ≥ p_b_tw_sc_LP_results["objective"]
@test (
    p_b_tw_ie_ch_LP_results["objective"]
    ≥ p_b_tw_ie_LP_results["objective"]
    ≥ p_b_tw_iec_LP_results["objective"]
    ≈ p_b_tw_iecw_LP_results["objective"]
    ≈ p_b_tw_iec_ch_LP_results["objective"]
    ≈ p_b_tw_iecw_ch_LP_results["objective"]
)
@test(
    p_b_tw_sc_LP_results["objective"]
    ≈ p_b_tw_iec_LP_results["objective"]
    ≈ p_b_tw_iecw_LP_results["objective"]
    ≈ p_b_tw_iec_ch_LP_results["objective"]
    ≈ p_b_tw_iecw_ch_LP_results["objective"]
)

@test (
    p_b_LP_results["objective"] 
    ≤ p_b_ch_LP_results["objective"]
    ≤ p_b_sc_LP_results["objective"]
)
@test p_b_s_LP_results["objective"] ≥ p_b_sc_LP_results["objective"]
@test (
    p_b_ie_ch_LP_results["objective"]
    ≥ p_b_ie_LP_results["objective"]
    ≥ p_b_iec_LP_results["objective"]
    ≈ p_b_iecw_LP_results["objective"]
    ≈ p_b_iec_ch_LP_results["objective"]
    ≈ p_b_iecw_ch_LP_results["objective"]
)
@test(
    p_b_sc_LP_results["objective"]
    ≈ p_b_iec_LP_results["objective"]
    ≈ p_b_iecw_LP_results["objective"]
    ≈ p_b_iec_ch_LP_results["objective"]
    ≈ p_b_iecw_ch_LP_results["objective"]
)

@test(
    p_o_LP_results["objective"] 
    ≤ p_o_sc_LP_results["objective"]
    ≤ p_o_scsc_LP_results["objective"]
    ≈ p_b_sc_LP_results["objective"]
)
@test p_o_s_LP_results["objective"] ≥ p_o_sc_LP_results["objective"]
@test p_o_ss_LP_results["objective"] ≥ p_o_scsc_LP_results["objective"]
@test(
    p_o_ch_LP_results["objective"] 
    ≤ p_o_sc_ch_LP_results["objective"]
    ≤ p_o_scsc_ch_LP_results["objective"]
)
@test p_o_s_ch_LP_results["objective"] ≥ p_o_sc_ch_LP_results["objective"]
@test p_o_ss_ch_LP_results["objective"] ≥ p_o_scsc_ch_LP_results["objective"]

@test p_b_LP_results["objective"] ≈ p_o_LP_results["objective"]

@test (
    sp_b_tw_LP_results["objective"] 
    ≤ sp_b_tw_ch_LP_results["objective"]
    ≤ sp_b_tw_sc_LP_results["objective"]
)
@test sp_b_tw_s_LP_results["objective"] ≥ sp_b_tw_sc_LP_results["objective"]
@test (
    sp_b_tw_ie_ch_LP_results["objective"]
    ≥ sp_b_tw_ie_LP_results["objective"]
    ≥ sp_b_tw_iec_LP_results["objective"]
    ≈ sp_b_tw_iecw_LP_results["objective"]
    ≈ sp_b_tw_iec_ch_LP_results["objective"]
    ≈ sp_b_tw_iecw_ch_LP_results["objective"]
)
@test(
    sp_b_tw_sc_LP_results["objective"]
    ≈ sp_b_tw_iec_LP_results["objective"]
    ≈ sp_b_tw_iecw_LP_results["objective"]
    ≈ sp_b_tw_iec_ch_LP_results["objective"]
    ≈ sp_b_tw_iecw_ch_LP_results["objective"]
)

@test (
    sp_b_LP_results["objective"] 
    ≤ sp_b_ch_LP_results["objective"]
    ≤ sp_b_sc_LP_results["objective"]
)
@test sp_b_s_LP_results["objective"] ≥ sp_b_sc_LP_results["objective"]
@test (
    sp_b_ie_ch_LP_results["objective"]
    ≥ sp_b_ie_LP_results["objective"]
    ≥ sp_b_iec_LP_results["objective"]
    ≈ sp_b_iecw_LP_results["objective"]
    ≈ sp_b_iec_ch_LP_results["objective"]
    ≈ sp_b_iecw_ch_LP_results["objective"]
)
@test(
    sp_b_sc_LP_results["objective"]
    ≈ sp_b_iec_LP_results["objective"]
    ≈ sp_b_iecw_LP_results["objective"]
    ≈ sp_b_iec_ch_LP_results["objective"]
    ≈ sp_b_iecw_ch_LP_results["objective"]
)


@test(
    sp_o_LP_results["objective"] 
    ≤ sp_o_sc_LP_results["objective"]
    ≤ sp_o_scsc_LP_results["objective"]
    ≈ sp_b_sc_LP_results["objective"]
)
@test (
    sp_o_s_LP_results["objective"] 
    ≥ sp_o_sc_LP_results["objective"]
    ≈ sp_o_sca_LP_results["objective"]
)
@test (
    sp_o_ss_LP_results["objective"] 
    ≥ sp_o_scsc_LP_results["objective"]
    ≈ sp_o_scsca_LP_results["objective"]
)
@test(
    sp_o_ch_LP_results["objective"] 
    ≤ sp_o_sc_ch_LP_results["objective"]
    ≤ sp_o_scsc_ch_LP_results["objective"]
    ≈ sp_b_sc_LP_results["objective"]
)
@test (
    sp_o_s_ch_LP_results["objective"] 
    ≥ sp_o_sc_ch_LP_results["objective"]
)
@test (
    sp_o_ss_ch_LP_results["objective"] 
    ≥ sp_o_scsc_ch_LP_results["objective"]
    ≈ sp_o_scsca_ch_LP_results["objective"]
)
@test sp_o_LP_results["objective"] ≤ sp_o_ch_LP_results["objective"]
@test sp_o_s_LP_results["objective"] ≤ sp_o_s_ch_LP_results["objective"]
@test sp_o_sc_LP_results["objective"] ≤ sp_o_sc_LP_results["objective"]
@test sp_o_sca_LP_results["objective"] ≤ sp_o_sca_ch_LP_results["objective"]
@test sp_o_ss_LP_results["objective"] ≤ sp_o_ss_ch_LP_results["objective"]
@test sp_o_scsc_LP_results["objective"] ≤ sp_o_scsc_ch_LP_results["objective"]
@test sp_o_scsca_LP_results["objective"] ≤ sp_o_scsca_ch_LP_results["objective"]

@test sp_b_LP_results["objective"] ≈ sp_o_LP_results["objective"]

@test p_b_LP_results["objective"] ≈ sp_b_LP_results["objective"]
@test p_b_sc_LP_results["objective"] ≈ sp_b_sc_LP_results["objective"]
@test p_o_LP_results["objective"] ≈ sp_o_LP_results["objective"]
@test p_o_s_LP_results["objective"] ≈ sp_o_s_LP_results["objective"]
@test p_o_sc_LP_results["objective"] ≈ sp_o_sc_LP_results["objective"]
@test p_o_scsc_LP_results["objective"] ≈ sp_o_scsc_LP_results["objective"]

@test p_b_ch_LP_results["objective"] ≥ sp_b_ch_LP_results["objective"]
@test p_o_ch_LP_results["objective"] ≥ sp_o_ch_LP_results["objective"]
@test p_o_s_ch_LP_results["objective"] ≥ sp_o_s_ch_LP_results["objective"]
@test p_o_sc_ch_LP_results["objective"] ≥ sp_o_sc_ch_LP_results["objective"]
@test p_o_scsc_ch_LP_results["objective"] ≥ sp_o_scsc_ch_LP_results["objective"]

p_b_tw_params["time_taken"]
p_b_tw_s_params["time_taken"]
p_b_tw_sc_params["time_taken"]
p_b_tw_ch_params["time_taken"]
p_b_tw_ie_params["time_taken"]
p_b_tw_iec_params["time_taken"]
p_b_tw_iecw_params["time_taken"]
p_b_tw_ie_hmaa_params["time_taken"]
p_b_tw_iec_hmaa_params["time_taken"]
p_b_tw_iecw_hmaa_params["time_taken"]
p_b_tw_ie_ch_params["time_taken"]
p_b_tw_iec_ch_params["time_taken"]
p_b_tw_iecw_ch_params["time_taken"]
p_b_tw_ie_hmaa_ch_params["time_taken"]
p_b_tw_iec_hmaa_ch_params["time_taken"]
p_b_tw_iecw_hmaa_ch_params["time_taken"]

p_b_params["time_taken"]
p_b_s_params["time_taken"]
p_b_sc_params["time_taken"]
p_b_ch_params["time_taken"]
p_b_ie_params["time_taken"]
p_b_iec_params["time_taken"]
p_b_iecw_params["time_taken"]
p_b_ie_hmaa_params["time_taken"]
p_b_iec_hmaa_params["time_taken"]
p_b_iecw_hmaa_params["time_taken"]
p_b_ie_ch_params["time_taken"]
p_b_iec_ch_params["time_taken"]
p_b_iecw_ch_params["time_taken"]
p_b_ie_hmaa_ch_params["time_taken"]
p_b_iec_hmaa_ch_params["time_taken"]
p_b_iecw_hmaa_ch_params["time_taken"]

p_o_params["time_taken"]
p_o_s_params["time_taken"]
p_o_sc_params["time_taken"]
p_o_ss_params["time_taken"]
p_o_scsc_params["time_taken"]
p_o_ch_params["time_taken"]
p_o_s_ch_params["time_taken"]
p_o_sc_ch_params["time_taken"]
p_o_ss_ch_params["time_taken"]
p_o_scsc_ch_params["time_taken"]

sp_b_tw_params["time_taken"]
sp_b_tw_s_params["time_taken"]
sp_b_tw_sc_params["time_taken"]
sp_b_tw_sca_params["time_taken"]
sp_b_tw_ch_params["time_taken"]
sp_b_tw_ie_params["time_taken"]
sp_b_tw_iec_params["time_taken"]
sp_b_tw_iecw_params["time_taken"]
sp_b_tw_ie_hmaa_params["time_taken"]
sp_b_tw_iec_hmaa_params["time_taken"]
sp_b_tw_iecw_hmaa_params["time_taken"]
sp_b_tw_ie_ch_params["time_taken"]
sp_b_tw_iec_ch_params["time_taken"]
sp_b_tw_iecw_ch_params["time_taken"]
sp_b_tw_ie_hmaa_ch_params["time_taken"]
sp_b_tw_iec_hmaa_ch_params["time_taken"]
sp_b_tw_iecw_hmaa_ch_params["time_taken"]

sp_b_params["time_taken"]
sp_b_s_params["time_taken"]
sp_b_sc_params["time_taken"]
sp_b_sca_params["time_taken"]
sp_b_ch_params["time_taken"]
sp_b_ie_params["time_taken"]
sp_b_iec_params["time_taken"]
sp_b_iecw_params["time_taken"]
sp_b_ie_hmaa_params["time_taken"]
sp_b_iec_hmaa_params["time_taken"]
sp_b_iecw_hmaa_params["time_taken"]
sp_b_ie_ch_params["time_taken"]
sp_b_iec_ch_params["time_taken"]
sp_b_iecw_ch_params["time_taken"]
sp_b_ie_hmaa_ch_params["time_taken"]
sp_b_iec_hmaa_ch_params["time_taken"]
sp_b_iecw_hmaa_ch_params["time_taken"]

sp_o_params["time_taken"]
sp_o_s_params["time_taken"]
sp_o_sc_params["time_taken"]
sp_o_sca_params["time_taken"]
sp_o_ss_params["time_taken"]
sp_o_scsc_params["time_taken"]
sp_o_scsca_params["time_taken"]
sp_o_ch_params["time_taken"]
sp_o_s_ch_params["time_taken"]
sp_o_sc_ch_params["time_taken"]
sp_o_sca_ch_params["time_taken"]
sp_o_ss_ch_params["time_taken"]
sp_o_scsc_ch_params["time_taken"]
sp_o_scsca_ch_params["time_taken"]

### Printouts

@printf("                                                       \t\tno 2-cycles\n")
@printf("path, benchmark:                                       %8.3f\t%8.3f\n", p_b_params["time_taken"], p_b_ch_params["time_taken"])
@printf("path, benchmark, elementary:                           %8.3f\t    ----\n", p_b_sc_params["time_taken"])
@printf("path, benchmark, boland:                               %8.3f\t%8.3f\n", p_b_iec_params["time_taken"], p_b_iec_ch_params["time_taken"])
@printf("path, ours:                                            %8.3f\t%8.3f\n", p_o_params["time_taken"], p_o_ch_params["time_taken"])
@printf("path, ours, elementary subpaths:                       %8.3f\t%8.3f\n", p_o_sc_params["time_taken"], p_o_sc_ch_params["time_taken"])
@printf("path, ours, elementary subpaths & paths:               %8.3f\t%8.3f\n", p_o_scsc_params["time_taken"], p_o_scsc_ch_params["time_taken"])

@printf("                                                       \t\tno 2-cycles\n")
@printf("subpath, benchmark:                                    %8.3f\t%8.3f\n", sp_b_params["time_taken"], sp_b_ch_params["time_taken"])
@printf("subpath, benchmark, elementary:                        %8.3f\t    ----\n", sp_b_sc_params["time_taken"])
@printf("subpath, benchmark, elementary (accel):                %8.3f\t    ----\n", sp_b_sca_params["time_taken"])
@printf("subpath, benchmark, boland:                            %8.3f\t%8.3f\n", sp_b_iec_params["time_taken"], sp_b_iec_ch_params["time_taken"])
@printf("subpath, ours:                                         %8.3f\t%8.3f\n", sp_o_params["time_taken"], sp_o_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths:                    %8.3f\t%8.3f\n", sp_o_sc_params["time_taken"], sp_o_sc_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths (accel):            %8.3f\t%8.3f\n", sp_o_sca_params["time_taken"], sp_o_sca_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths & paths:            %8.3f\t%8.3f\n", sp_o_scsc_params["time_taken"], sp_o_scsc_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths & paths (accel):    %8.3f\t%8.3f\n", sp_o_scsca_params["time_taken"], sp_o_scsca_ch_params["time_taken"])


@printf("                                                       \t\t\tno 2-cycles\n")
@printf("path, benchmark:                                       %8.1f\t\t%8.1f\n", p_b_LP_results["objective"], p_b_ch_LP_results["objective"])
@printf("path, benchmark, elementary:                           %8.1f\t\t    ----\n", p_b_sc_LP_results["objective"])
@printf("path, benchmark, boland:                               %8.1f\t\t%8.1f\n", p_b_iec_LP_results["objective"], p_b_iec_ch_LP_results["objective"])
@printf("path, ours:                                            %8.1f\t\t%8.1f\n", p_o_LP_results["objective"], p_o_ch_LP_results["objective"])
@printf("path, ours, elementary subpaths:                       %8.1f\t\t%8.1f\n", p_o_sc_LP_results["objective"], p_o_sc_ch_LP_results["objective"])
@printf("path, ours, elementary subpaths & paths:               %8.1f\t\t%8.1f\n", p_o_scsc_LP_results["objective"], p_o_scsc_ch_LP_results["objective"])

@printf("                                                       \t\t\tno 2-cycles\n")
@printf("subpath, benchmark:                                    %8.1f\t\t%8.1f\n", sp_b_LP_results["objective"], sp_b_ch_LP_results["objective"])
@printf("subpath, benchmark, elementary:                        %8.1f\t\t    ----\n", sp_b_sc_LP_results["objective"])
@printf("subpath, benchmark, elementary (accel):                %8.1f\t\t    ----\n", sp_b_sca_LP_results["objective"])
@printf("subpath, benchmark, boland:                            %8.1f\t\t%8.1f\n", sp_b_iec_LP_results["objective"], sp_b_iec_ch_LP_results["objective"])
@printf("subpath, ours:                                         %8.1f\t\t%8.1f\n", sp_o_LP_results["objective"], sp_o_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths:                    %8.1f\t\t%8.1f\n", sp_o_sc_LP_results["objective"], sp_o_sc_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths (accel):            %8.1f\t\t%8.1f\n", sp_o_sca_LP_results["objective"], sp_o_sca_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths:            %8.1f\t\t%8.1f\n", sp_o_scsc_LP_results["objective"], sp_o_scsc_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths (accel):    %8.1f\t\t%8.1f\n", sp_o_scsca_LP_results["objective"], sp_o_scsca_ch_LP_results["objective"])

### Scratch work

include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames

using Test

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 16,
    n_charging = 9,
    n_vehicles = 7,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "grid",
    shrinkage_depots = 1.0,
    shrinkage_charging = 0.7,
    T = 40000,
    seed = 5,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 4,
    permissiveness = 0.2,
)
G = construct_graph(data)
# plot_instance(data)

Env = nothing
method = "benchmark"
time_windows = false
subpath_single_service = false
subpath_check_customers = false
path_single_service = false
path_check_customers = false
incremental_elementarity = true
warm_start = false
christofides = true
verbose = true

### Beginning debug of path_formulation_column_generation()
compute_minimum_time_to_nearest_depot!(data, G)
compute_minimum_charge_to_nearest_depot_charging_station!(data, G)

some_paths = generate_artificial_paths(data)
path_costs = compute_path_costs(
    data, 
    some_paths,
)
path_service = compute_path_service(
    data, 
    some_paths,
)
mp_results = Dict()
params = Dict()
params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
params["objective"] = Float64[]
params["κ"] = Dict{Int, Float64}[]
params["μ"] = Dict{Int, Float64}[]
params["ν"] = Vector{Float64}[]
params["lp_relaxation_solution_time_taken"] = Float64[]
params["sp_base_time_taken"] = Float64[]
params["sp_full_time_taken"] = Float64[]
params["sp_total_time_taken"] = Float64[]
params["lp_relaxation_constraint_time_taken"] = Float64[]
params["number_of_new_paths"] = Int[]

printlist = String[]
counter = 0
converged = false

if isnothing(Env)
    mp_model = @suppress Model(Gurobi.Optimizer)
else
    mp_model = @suppress Model(() -> Gurobi.Optimizer(Env))
end
JuMP.set_attribute(mp_model, "MIPGapAbs", 1e-3)
JuMP.set_string_names_on_creation(mp_model, false)
z = Dict{
    Tuple{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Int
    }, 
    VariableRef
}(
    (key, p) => @variable(mp_model, lower_bound = 0)
    for key in keys(some_paths)
        for p in 1:length(some_paths[key])
)
@constraint(
    mp_model,
    κ[i in data["N_depots"]],
    sum(
        sum(
            z[((i,0,data["B"]),state2),p]
            for p in 1:length(some_paths[((i,0,data["B"]),state2)])
        )        
        for (state1, state2) in keys(some_paths)
            if state1[1] == i && state1[2] == 0 && state1[3] == data["B"]
    )
    == data["v_start"][findfirst(x -> (x == i), data["N_depots"])]
)
@constraint(
    mp_model,
    μ[n2 in data["N_depots"]],
    sum(
        sum(
            z[(state1, state2),p]
            for p in 1:length(some_paths[(state1, state2)])
        )
        for (state1, state2) in keys(some_paths)
            if state2[1] == n2
    ) ≥ data["v_end"][n2]
)
@constraint(
    mp_model,
    ν[j in data["N_customers"]],
    sum(
        sum(
            path_service[((state1, state2),j)][p] * z[(state1, state2),p]
            for p in 1:length(some_paths[(state1, state2)])
        )
        for (state1, state2) in keys(some_paths)
    ) == 1
)
@expression(
    mp_model,
    path_costs_expr,
    sum(
        sum(
            path_costs[state_pair][p] * z[state_pair,p]
            for p in 1:length(some_paths[state_pair])
        )
        for state_pair in keys(some_paths)
    )
)
@objective(mp_model, Min, path_costs_expr)

# while (
#     !converged
#     && time_limit ≥ (time() - start_time)
# )
begin
    counter += 1
    mp_solution_start_time = time()
    @suppress optimize!(mp_model)
    mp_solution_end_time = time()
    mp_results = Dict(
        "model" => mp_model,
        "objective" => objective_value(mp_model),
        "z" => Dict(
            (key, p) => value.(z[(key, p)])
            for (key, p) in keys(z)
        ),
        "κ" => Dict(zip(data["N_depots"], dual.(mp_model[:κ]).data)),
        "μ" => Dict(zip(data["N_depots"], dual.(mp_model[:μ]).data)),
        "ν" => dual.(mp_model[:ν]).data,
        "solution_time_taken" => round(mp_solution_end_time - mp_solution_start_time, digits = 3),
    )
    push!(params["objective"], mp_results["objective"])
    push!(params["κ"], mp_results["κ"])
    push!(params["μ"], mp_results["μ"])
    push!(params["ν"], mp_results["ν"])
    push!(params["lp_relaxation_solution_time_taken"], mp_results["solution_time_taken"])




    if method == "ours"
        (negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
            G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
            ;
            subpath_single_service = subpath_single_service,
            subpath_check_customers = subpath_check_customers,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
        )
        (generated_paths) = get_paths_from_negative_path_labels(
            data, negative_full_labels,
        )
        push!(
            params["sp_base_time_taken"],
            round(base_labels_time, digits=3)
        )
        push!(
            params["sp_full_time_taken"],
            round(full_labels_time, digits=3)
        )
        push!(
            params["sp_total_time_taken"],
            round(base_labels_time + full_labels_time, digits=3)
        )
    elseif method == "benchmark"
        if incremental_elementarity
            (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark_incremental_elementarity(
                G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                ;
                time_windows = time_windows,
                path_check_customers = path_check_customers,
                warm_start = warm_start,
                christofides = christofides,
                verbose = verbose,
                rule = "hmaa",
            )
        else
            (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
                G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
                ;
                time_windows = time_windows,
                path_single_service = path_single_service,
                path_check_customers = path_check_customers,
                christofides = christofides,
            )
        end
        generated_paths = get_paths_from_negative_pure_path_labels(
            data, negative_pure_path_labels,
        )
        push!(
            params["sp_base_time_taken"],
            0.0
        )
        push!(
            params["sp_full_time_taken"],
            round(pure_path_labels_time, digits=3)
        )
        push!(
            params["sp_total_time_taken"],
            round(pure_path_labels_time, digits=3)
        )
    end
    
    if length(generated_paths) == 0
        push!(params["number_of_new_paths"], 0)
        converged = true
    else
        push!(
            params["number_of_new_paths"],
            sum(length(v) for v in values(generated_paths))
        )
    end

    mp_constraint_start_time = time()
    for state_pair in keys(generated_paths)
        if !(state_pair in keys(some_paths))
            some_paths[state_pair] = []
            path_costs[state_pair] = []
            for i in 1:data["n_customers"]
                path_service[(state_pair, i)] = []
            end
            count = 0
        else
            count = length(some_paths[state_pair])
        end
        for p_new in generated_paths[state_pair]
            if state_pair in keys(some_paths)
                add = !any(isequal(p_new, s) for s in some_paths[state_pair])
            else
                add = true
            end
            if add
                # 1: include in some_paths
                push!(some_paths[state_pair], p_new)
                # 2: add path cost
                push!(
                    path_costs[state_pair], 
                    compute_path_cost(data, p_new)
                )
                # 3: add path service
                for i in 1:data["n_customers"]
                    push!(path_service[(state_pair, i)], p_new.served[i])
                end
                # 4: create variable
                count += 1
                z[(state_pair, count)] = @variable(mp_model, lower_bound = 0)
                (state1, state2) = state_pair
                # 5: modify constraints starting from depot, ending at depot
                set_normalized_coefficient(κ[state1[1]], z[state_pair,count], 1)
                set_normalized_coefficient(μ[state2[1]], z[state_pair,count], 1)
                # 6: modify customer service constraints
                for l in data["N_customers"]
                    set_normalized_coefficient(ν[l], z[state_pair, count], p_new.served[l])
                end
                # 7: modify objective
                set_objective_coefficient(mp_model, z[state_pair, count], path_costs[state_pair][count])
            end
        end
    end
    mp_constraint_end_time = time()

    push!(
        params["number_of_paths"], 
        sum(length(v) for v in values(some_paths))
    )
    push!(
        params["lp_relaxation_constraint_time_taken"],
        round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
    )
end
params["number_of_paths"]
params["sp_full_time_taken"]


### Starting debug of subproblem_iteration_benchmark_incremental_elementarity

κ, μ, ν = mp_results["κ"], mp_results["μ"], mp_results["ν"]
time_windows = false
path_check_customers = true
verbose = true 
warm_start = false 
christofides = true

start_time = time()
S = Int[]
initial_pure_path_labels = nothing

pure_path_labels = find_nondominated_paths_S(
    G, data, κ, μ, ν,
    ;
    S = S,
    time_windows = time_windows, 
    check_customers = path_check_customers,
    christofides = christofides,
    initial_pure_path_labels = initial_pure_path_labels
)
pure_path_labels_count = sum(
    length(pure_path_labels[starting_node][end_node]) 
    for starting_node in data["N_depots"]
        for end_node in keys(pure_path_labels[starting_node]) 
)
# filter for path labels: (i) ending at depot, 
# (ii) w/ negative reduced cost,
# (iii) elementary
d_pure_path_labels = get_depot_pure_path_labels(data, pure_path_labels);
n_d_pure_path_labels = get_negative_pure_path_labels_from_pure_path_labels(data, d_pure_path_labels);
n_d_pure_path_labels_count = sum(
    length(n_d_pure_path_labels[starting_node][end_node]) 
    for starting_node in data["N_depots"]
        for end_node in keys(n_d_pure_path_labels[starting_node]) 
)
if n_d_pure_path_labels_count == 0
    # this activates if the overall CG converges
    return (n_d_pure_path_labels, 0, time() - start_time)
end
(e_n_d_pure_path_labels, ne_n_d_pure_path_labels) = get_elementary_nonelementary_pure_path_labels(
    data, n_d_pure_path_labels; 
    S = data["N_customers"]
);
e_n_d_pure_path_labels_count = sum(
    length(e_n_d_pure_path_labels[starting_node][end_node]) 
    for starting_node in data["N_depots"]
        for end_node in keys(e_n_d_pure_path_labels[starting_node]) 
)
if e_n_d_pure_path_labels_count > 0
    # this activates if there is some negative path found (elementary w.r.t. all customers)
    return (e_n_d_pure_path_labels, e_n_d_pure_path_labels_count, time() - start_time)
end

# else, expand S
# println(length(S))
if length(S) == data["n_customers"]
    error()
end

max_freq = 0
max_custs = Int[]
for starting_node in data["N_depots"]
    for end_node in keys(n_d_pure_path_labels[starting_node])
        for (key, path_label) in pairs(n_d_pure_path_labels[starting_node][end_node])
            freq = maximum(path_label.served)
            custs = findall(==(freq), path_label.served)
            if freq < max_freq
                continue
            elseif freq > max_freq
                max_freq = freq
                max_custs = custs 
            else
                max_custs = union(max_custs, custs)
            end
            println("$max_freq, $max_custs")
        end
    end
end


S = sort(union(S, max_custs))