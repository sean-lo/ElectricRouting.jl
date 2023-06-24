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
    seed = 6,
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

p_b_LP_results, p_b_IP_results, p_b_params, p_b_printlist, p_b_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true)
p_b_s_LP_results, p_b_s_IP_results, p_b_s_params, p_b_s_printlist, p_b_s_some_paths = path_formulation_column_generation(G, data; method = "benchmark", path_single_service = true, verbose = true)
p_b_sc_LP_results, p_b_sc_IP_results, p_b_sc_params, p_b_sc_printlist, p_b_sc_some_paths = path_formulation_column_generation(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)

p_b_ch_LP_results, p_b_ch_IP_results, p_b_ch_params, p_b_ch_printlist, p_b_ch_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, christofides = true)

p_o_LP_results, p_o_IP_results, p_o_params, p_o_printlist, p_o_some_paths = path_formulation_column_generation(G, data; method = "ours", verbose = true)
p_o_s_LP_results, p_o_s_IP_results, p_o_s_params, p_o_s_printlist, p_o_s_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, verbose = true)
p_o_sc_LP_results, p_o_sc_IP_results, p_o_sc_params, p_o_sc_printlist, p_o_sc_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
p_o_ss_LP_results, p_o_ss_IP_results, p_o_ss_params, p_o_ss_printlist, p_o_ss_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
p_o_scsc_LP_results, p_o_scsc_IP_results, p_o_scsc_params, p_o_scsc_printlist, p_o_scsc_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)
p_o_ngs_LP_results, p_o_ngs_IP_results, p_o_ngs_params, p_o_ngs_printlist, p_o_ngs_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
p_o_ngl_LP_results, p_o_ngl_IP_results, p_o_ngl_params, p_o_ngl_printlist, p_o_ngl_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)
p_o_ngsa_LP_results, p_o_ngsa_IP_results, p_o_ngsa_params, p_o_ngsa_printlist, p_o_ngsa_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
p_o_ngla_LP_results, p_o_ngla_IP_results, p_o_ngla_params, p_o_ngla_printlist, p_o_ngla_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)


p_o_ch_LP_results, p_o_ch_IP_results, p_o_ch_params, p_o_ch_printlist, p_o_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", verbose = true, christofides = true)
p_o_s_ch_LP_results, p_o_s_ch_IP_results, p_o_s_ch_params, p_o_s_ch_printlist, p_o_s_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, verbose = true, christofides = true)
p_o_sc_ch_LP_results, p_o_sc_ch_IP_results, p_o_sc_ch_params, p_o_sc_ch_printlist, p_o_sc_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true, christofides = true)
p_o_ss_ch_LP_results, p_o_ss_ch_IP_results, p_o_ss_ch_params, p_o_ss_ch_printlist, p_o_ss_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true, christofides = true)
p_o_scsc_ch_LP_results, p_o_scsc_ch_IP_results, p_o_scsc_ch_params, p_o_scsc_ch_printlist, p_o_scsc_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true, christofides = true)
p_o_ngs_ch_LP_results, p_o_ngs_ch_IP_results, p_o_ngs_ch_params, p_o_ngs_ch_printlist, p_o_ngs_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true)
p_o_ngl_ch_LP_results, p_o_ngl_ch_IP_results, p_o_ngl_ch_params, p_o_ngl_ch_printlist, p_o_ngl_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true)
p_o_ngsa_ch_LP_results, p_o_ngsa_ch_IP_results, p_o_ngsa_ch_params, p_o_ngsa_ch_printlist, p_o_ngsa_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true)
p_o_ngla_ch_LP_results, p_o_ngla_ch_IP_results, p_o_ngla_ch_params, p_o_ngla_ch_printlist, p_o_ngla_ch_some_paths = path_formulation_column_generation(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true)

sp_b_LP_results, sp_b_IP_results, sp_b_params, sp_b_printlist, sp_b_some_subpaths, sp_b_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", verbose = true)
sp_b_s_LP_results, sp_b_s_IP_results, sp_b_s_params, sp_b_s_printlist, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, verbose = true)
sp_b_sc_LP_results, sp_b_sc_IP_results, sp_b_sc_params, sp_b_sc_printlist, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)
sp_b_sca_LP_results, sp_b_sca_IP_results, sp_b_sca_params, sp_b_sca_printlist, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

sp_b_ch_LP_results, sp_b_ch_IP_results, sp_b_ch_params, sp_b_ch_printlist, sp_b_ch_some_subpaths, sp_b_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", verbose = true, christofides = true)

sp_o_LP_results, sp_o_IP_results, sp_o_params, sp_o_printlist, sp_o_some_subpaths, sp_o_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", verbose = true)
sp_o_s_LP_results, sp_o_s_IP_results, sp_o_s_params, sp_o_s_printlist, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, verbose = true)
sp_o_sc_LP_results, sp_o_sc_IP_results, sp_o_sc_params, sp_o_sc_printlist, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
sp_o_sca_LP_results, sp_o_sca_IP_results, sp_o_sca_params, sp_o_sca_printlist, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true)
sp_o_ss_LP_results, sp_o_ss_IP_results, sp_o_ss_params, sp_o_ss_printlist, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
sp_o_scsc_LP_results, sp_o_scsc_IP_results, sp_o_scsc_params, sp_o_scsc_printlist, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)
sp_o_scsca_LP_results, sp_o_scsca_IP_results, sp_o_scsca_params, sp_o_scsca_printlist, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)
sp_o_ngs_LP_results, sp_o_ngs_IP_results, sp_o_ngs_params, sp_o_ngs_printlist, sp_o_ngs_some_subpaths, sp_o_ngs_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
sp_o_ngl_LP_results, sp_o_ngl_IP_results, sp_o_ngl_params, sp_o_ngl_printlist, sp_o_ngl_some_subpaths, sp_o_ngl_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)
sp_o_ngsa_LP_results, sp_o_ngsa_IP_results, sp_o_ngsa_params, sp_o_ngsa_printlist, sp_o_ngsa_some_subpaths, sp_o_ngsa_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
sp_o_ngla_LP_results, sp_o_ngla_IP_results, sp_o_ngla_params, sp_o_ngla_printlist, sp_o_ngla_some_subpaths, sp_o_ngla_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)

sp_o_ch_LP_results, sp_o_ch_IP_results, sp_o_ch_params, sp_o_ch_printlist, sp_o_ch_some_subpaths, sp_o_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", verbose = true, christofides = true)
sp_o_s_ch_LP_results, sp_o_s_ch_IP_results, sp_o_s_ch_params, sp_o_s_ch_printlist, sp_o_s_ch_some_subpaths, sp_o_s_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, verbose = true, christofides = true)
sp_o_sc_ch_LP_results, sp_o_sc_ch_IP_results, sp_o_sc_ch_params, sp_o_sc_ch_printlist, sp_o_sc_ch_some_subpaths, sp_o_sc_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true, christofides = true)
sp_o_sca_ch_LP_results, sp_o_sca_ch_IP_results, sp_o_sca_ch_params, sp_o_sca_ch_printlist, sp_o_sca_ch_some_subpaths, sp_o_sca_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true, christofides = true)
sp_o_ss_ch_LP_results, sp_o_ss_ch_IP_results, sp_o_ss_ch_params, sp_o_ss_ch_printlist, sp_o_ss_ch_some_subpaths, sp_o_ss_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true, christofides = true)
sp_o_scsc_ch_LP_results, sp_o_scsc_ch_IP_results, sp_o_scsc_ch_params, sp_o_scsc_ch_printlist, sp_o_scsc_ch_some_subpaths, sp_o_scsc_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true, christofides = true)
sp_o_scsca_ch_LP_results, sp_o_scsca_ch_IP_results, sp_o_scsca_ch_params, sp_o_scsca_ch_printlist, sp_o_scsca_ch_some_subpaths, sp_o_scsca_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true, christofides = true)
sp_o_ngs_ch_LP_results, sp_o_ngs_ch_IP_results, sp_o_ngs_ch_params, sp_o_ngs_ch_printlist, sp_o_ngs_ch_some_subpaths, sp_o_ngs_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", christofides = true, verbose = true)
sp_o_ngl_ch_LP_results, sp_o_ngl_ch_IP_results, sp_o_ngl_ch_params, sp_o_ngl_ch_printlist, sp_o_ngl_ch_some_subpaths, sp_o_ngl_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", christofides = true, verbose = true)
sp_o_ngsa_ch_LP_results, sp_o_ngsa_ch_IP_results, sp_o_ngsa_ch_params, sp_o_ngsa_ch_printlist, sp_o_ngsa_ch_some_subpaths, sp_o_ngsa_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "small", christofides = true, verbose = true)
sp_o_ngla_ch_LP_results, sp_o_ngla_ch_IP_results, sp_o_ngla_ch_params, sp_o_ngla_ch_printlist, sp_o_ngla_ch_some_subpaths, sp_o_ngla_ch_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", ngroute = true, ngroute_alt = true, ngroute_neighborhood_charging_depots_size = "large", christofides = true, verbose = true)


collect_path_solution_metrics!(p_b_LP_results, data, p_b_some_paths)
collect_path_solution_metrics!(p_b_s_LP_results, data, p_b_s_some_paths)
collect_path_solution_metrics!(p_b_sc_LP_results, data, p_b_sc_some_paths)
collect_path_solution_metrics!(p_b_ch_LP_results, data, p_b_ch_some_paths)

collect_path_solution_metrics!(p_o_LP_results, data, p_o_some_paths)
collect_path_solution_metrics!(p_o_s_LP_results, data, p_o_s_some_paths)
collect_path_solution_metrics!(p_o_sc_LP_results, data, p_o_sc_some_paths)
collect_path_solution_metrics!(p_o_ss_LP_results, data, p_o_ss_some_paths)
collect_path_solution_metrics!(p_o_scsc_LP_results, data, p_o_scsc_some_paths)
collect_path_solution_metrics!(p_o_ngs_LP_results, data, p_o_ngs_some_paths)
collect_path_solution_metrics!(p_o_ngl_LP_results, data, p_o_ngl_some_paths)
collect_path_solution_metrics!(p_o_ngsa_LP_results, data, p_o_ngsa_some_paths)
collect_path_solution_metrics!(p_o_ngla_LP_results, data, p_o_ngla_some_paths)

collect_path_solution_metrics!(p_o_ch_LP_results, data, p_o_ch_some_paths)
collect_path_solution_metrics!(p_o_s_ch_LP_results, data, p_o_s_ch_some_paths)
collect_path_solution_metrics!(p_o_sc_ch_LP_results, data, p_o_sc_ch_some_paths)
collect_path_solution_metrics!(p_o_ss_ch_LP_results, data, p_o_ss_ch_some_paths)
collect_path_solution_metrics!(p_o_scsc_ch_LP_results, data, p_o_scsc_ch_some_paths)
collect_path_solution_metrics!(p_o_ngs_ch_LP_results, data, p_o_ngs_ch_some_paths)
collect_path_solution_metrics!(p_o_ngl_ch_LP_results, data, p_o_ngl_ch_some_paths)
collect_path_solution_metrics!(p_o_ngsa_ch_LP_results, data, p_o_ngsa_ch_some_paths)
collect_path_solution_metrics!(p_o_ngla_ch_LP_results, data, p_o_ngla_ch_some_paths)

collect_subpath_solution_metrics!(sp_b_LP_results, data, sp_b_some_subpaths, sp_b_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_s_LP_results, data, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sc_LP_results, data, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sca_LP_results, data, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ch_LP_results, data, sp_b_ch_some_subpaths, sp_b_ch_some_charging_arcs)

collect_subpath_solution_metrics!(sp_o_LP_results, data, sp_o_some_subpaths, sp_o_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_s_LP_results, data, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sc_LP_results, data, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sca_LP_results, data, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ss_LP_results, data, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsc_LP_results, data, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsca_LP_results, data, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngs_LP_results, data, sp_o_ngs_some_subpaths, sp_o_ngs_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngl_LP_results, data, sp_o_ngl_some_subpaths, sp_o_ngl_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngsa_LP_results, data, sp_o_ngsa_some_subpaths, sp_o_ngsa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngla_LP_results, data, sp_o_ngla_some_subpaths, sp_o_ngla_some_charging_arcs)

collect_subpath_solution_metrics!(sp_o_ch_LP_results, data, sp_o_ch_some_subpaths, sp_o_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_s_ch_LP_results, data, sp_o_s_ch_some_subpaths, sp_o_s_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sc_ch_LP_results, data, sp_o_sc_ch_some_subpaths, sp_o_sc_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sca_ch_LP_results, data, sp_o_sca_ch_some_subpaths, sp_o_sca_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ss_ch_LP_results, data, sp_o_ss_ch_some_subpaths, sp_o_ss_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsc_ch_LP_results, data, sp_o_scsc_ch_some_subpaths, sp_o_scsc_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsca_ch_LP_results, data, sp_o_scsca_ch_some_subpaths, sp_o_scsca_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngs_ch_LP_results, data, sp_o_ngs_ch_some_subpaths, sp_o_ngs_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngl_ch_LP_results, data, sp_o_ngl_ch_some_subpaths, sp_o_ngl_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngsa_ch_LP_results, data, sp_o_ngsa_ch_some_subpaths, sp_o_ngsa_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ngla_ch_LP_results, data, sp_o_ngla_ch_some_subpaths, sp_o_ngla_ch_some_charging_arcs)


p_b_LP_results["objective"]
p_b_s_LP_results["objective"]
p_b_sc_LP_results["objective"]
p_b_ch_LP_results["objective"]

p_o_LP_results["objective"]
p_o_s_LP_results["objective"]
p_o_sc_LP_results["objective"]
p_o_ss_LP_results["objective"]
p_o_scsc_LP_results["objective"]
p_o_ngs_LP_results["objective"]
p_o_ngl_LP_results["objective"]
p_o_ngsa_LP_results["objective"]
p_o_ngla_LP_results["objective"]
p_o_ch_LP_results["objective"]
p_o_s_ch_LP_results["objective"]
p_o_sc_ch_LP_results["objective"]
p_o_ss_ch_LP_results["objective"]
p_o_scsc_ch_LP_results["objective"]
p_o_ngs_ch_LP_results["objective"]
p_o_ngl_ch_LP_results["objective"]
p_o_ngsa_ch_LP_results["objective"]
p_o_ngla_ch_LP_results["objective"]


sp_b_LP_results["objective"]
sp_b_s_LP_results["objective"]
sp_b_sc_LP_results["objective"]
sp_b_sca_LP_results["objective"]
sp_b_ch_LP_results["objective"]

sp_o_LP_results["objective"]
sp_o_s_LP_results["objective"]
sp_o_sc_LP_results["objective"]
sp_o_sca_LP_results["objective"]
sp_o_ss_LP_results["objective"]
sp_o_scsc_LP_results["objective"]
sp_o_scsca_LP_results["objective"]
sp_o_ngs_LP_results["objective"]
sp_o_ngl_LP_results["objective"]
sp_o_ngsa_LP_results["objective"]
sp_o_ngla_LP_results["objective"]
sp_o_ch_LP_results["objective"]
sp_o_s_ch_LP_results["objective"]
sp_o_sc_ch_LP_results["objective"]
sp_o_sca_ch_LP_results["objective"]
sp_o_ss_ch_LP_results["objective"]
sp_o_scsc_ch_LP_results["objective"]
sp_o_scsca_ch_LP_results["objective"]
sp_o_ngs_ch_LP_results["objective"]
sp_o_ngl_ch_LP_results["objective"]
sp_o_ngsa_ch_LP_results["objective"]
sp_o_ngla_ch_LP_results["objective"]


@test (
    p_b_LP_results["objective"] 
    ≤ p_b_ch_LP_results["objective"]
    ≤ p_b_sc_LP_results["objective"]
)
@test p_b_s_LP_results["objective"] ≥ p_b_sc_LP_results["objective"]

@test(
    p_o_LP_results["objective"] 
    ≤ p_o_sc_LP_results["objective"]
    ≤ p_o_scsc_LP_results["objective"]
    ≈ p_b_sc_LP_results["objective"]
)
@test(
    p_o_LP_results["objective"] 
    ≤ p_o_ngs_LP_results["objective"] 
    ≤ p_o_ngl_LP_results["objective"] 
    ≤ p_o_scsc_LP_results["objective"]
)
@test p_o_ngs_LP_results["objective"] ≈ p_o_ngsa_LP_results["objective"]
@test p_o_ngl_LP_results["objective"] ≈ p_o_ngla_LP_results["objective"]
@test p_o_s_LP_results["objective"] ≥ p_o_sc_LP_results["objective"]
@test p_o_ss_LP_results["objective"] ≥ p_o_scsc_LP_results["objective"]
@test(
    p_o_ch_LP_results["objective"] 
    ≤ p_o_sc_ch_LP_results["objective"]
    ≤ p_o_scsc_ch_LP_results["objective"]
)
@test(
    p_o_ch_LP_results["objective"] 
    ≤ p_o_ngs_ch_LP_results["objective"]
    ≤ p_o_ngl_ch_LP_results["objective"]
    ≤ p_o_scsc_ch_LP_results["objective"]
)
@test p_o_ngs_ch_LP_results["objective"] ≈ p_o_ngsa_ch_LP_results["objective"]
@test p_o_ngl_ch_LP_results["objective"] ≈ p_o_ngla_ch_LP_results["objective"]
@test p_o_s_ch_LP_results["objective"] ≥ p_o_sc_ch_LP_results["objective"]
@test p_o_ss_ch_LP_results["objective"] ≥ p_o_scsc_ch_LP_results["objective"]

@test p_b_LP_results["objective"] ≈ p_o_LP_results["objective"]


@test (
    sp_b_LP_results["objective"] 
    ≤ sp_b_ch_LP_results["objective"]
    ≤ sp_b_sc_LP_results["objective"]
)
@test sp_b_s_LP_results["objective"] ≥ sp_b_sc_LP_results["objective"]


@test(
    sp_o_LP_results["objective"] 
    ≤ sp_o_sc_LP_results["objective"]
    ≤ sp_o_scsc_LP_results["objective"]
    ≈ sp_b_sc_LP_results["objective"]
)
@test(
    sp_o_LP_results["objective"] 
    ≤ sp_o_ngs_LP_results["objective"]
    ≤ sp_o_ngl_LP_results["objective"]
    ≤ sp_o_scsc_LP_results["objective"]
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
    ≈ sp_b_sc_ch_LP_results["objective"]
)
@test(
    sp_o_ch_LP_results["objective"] 
    ≤ sp_o_ngs_ch_LP_results["objective"]
    ≤ sp_o_ngl_ch_LP_results["objective"]
    ≤ sp_o_scsc_ch_LP_results["objective"]
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
@test sp_o_ngs_LP_results["objective"] ≈ sp_o_ngsa_LP_results["objective"]
@test sp_o_ngl_LP_results["objective"] ≈ sp_o_ngla_LP_results["objective"]
@test sp_o_ngs_ch_LP_results["objective"] ≈ sp_o_ngsa_ch_LP_results["objective"]
@test sp_o_ngl_ch_LP_results["objective"] ≈ sp_o_ngla_ch_LP_results["objective"]
@test sp_o_LP_results["objective"] ≤ sp_o_ch_LP_results["objective"]
@test sp_o_s_LP_results["objective"] ≤ sp_o_s_ch_LP_results["objective"]
@test sp_o_sc_LP_results["objective"] ≤ sp_o_sc_ch_LP_results["objective"]
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
@test p_o_ngs_LP_results["objective"] ≈ sp_o_ngs_LP_results["objective"]
@test p_o_ngl_LP_results["objective"] ≈ sp_o_ngl_LP_results["objective"]


@test p_b_ch_LP_results["objective"] ≥ sp_b_ch_LP_results["objective"]
@test p_o_ch_LP_results["objective"] ≥ sp_o_ch_LP_results["objective"]
@test p_o_s_ch_LP_results["objective"] ≥ sp_o_s_ch_LP_results["objective"]
@test p_o_sc_ch_LP_results["objective"] ≥ sp_o_sc_ch_LP_results["objective"]
@test p_o_scsc_ch_LP_results["objective"] ≥ sp_o_scsc_ch_LP_results["objective"]
@test p_o_ngs_ch_LP_results["objective"] ≥ sp_o_ngs_ch_LP_results["objective"]
@test p_o_ngl_ch_LP_results["objective"] ≥ sp_o_ngl_ch_LP_results["objective"]


p_b_params["time_taken"]
p_b_s_params["time_taken"]
p_b_sc_params["time_taken"]
p_b_ch_params["time_taken"]

p_o_params["time_taken"]
p_o_s_params["time_taken"]
p_o_sc_params["time_taken"]
p_o_ss_params["time_taken"]
p_o_scsc_params["time_taken"]
p_o_ngs_params["time_taken"]
p_o_ngl_params["time_taken"]
p_o_ngsa_params["time_taken"]
p_o_ngla_params["time_taken"]

p_o_ch_params["time_taken"]
p_o_s_ch_params["time_taken"]
p_o_sc_ch_params["time_taken"]
p_o_ss_ch_params["time_taken"]
p_o_scsc_ch_params["time_taken"]
p_o_ngs_ch_params["time_taken"]
p_o_ngl_ch_params["time_taken"]
p_o_ngsa_ch_params["time_taken"]
p_o_ngla_ch_params["time_taken"]


sp_b_params["time_taken"]
sp_b_s_params["time_taken"]
sp_b_sc_params["time_taken"]
sp_b_sca_params["time_taken"]
sp_b_ch_params["time_taken"]

sp_o_params["time_taken"]
sp_o_s_params["time_taken"]
sp_o_sc_params["time_taken"]
sp_o_sca_params["time_taken"]
sp_o_ss_params["time_taken"]
sp_o_scsc_params["time_taken"]
sp_o_scsca_params["time_taken"]
sp_o_ngs_params["time_taken"]
sp_o_ngl_params["time_taken"]
sp_o_ngsa_params["time_taken"]
sp_o_ngla_params["time_taken"]

sp_o_ch_params["time_taken"]
sp_o_s_ch_params["time_taken"]
sp_o_sc_ch_params["time_taken"]
sp_o_sca_ch_params["time_taken"]
sp_o_ss_ch_params["time_taken"]
sp_o_scsc_ch_params["time_taken"]
sp_o_scsca_ch_params["time_taken"]
sp_o_ngs_ch_params["time_taken"]
sp_o_ngl_ch_params["time_taken"]
sp_o_ngsa_ch_params["time_taken"]
sp_o_ngla_ch_params["time_taken"]

### Printouts

@printf("                                                               \t\tno 2-cycles\n")
@printf("path, benchmark:                                               %8.3f\t%8.3f\n", p_b_params["time_taken"], p_b_ch_params["time_taken"])
@printf("path, benchmark, elementary:                                   %8.3f\t    ----\n", p_b_sc_params["time_taken"])
@printf("path, ours:                                                    %8.3f\t%8.3f\n", p_o_params["time_taken"], p_o_ch_params["time_taken"])
@printf("path, ours, elementary subpaths:                               %8.3f\t%8.3f\n", p_o_sc_params["time_taken"], p_o_sc_ch_params["time_taken"])
@printf("path, ours, elementary subpaths & paths:                       %8.3f\t%8.3f\n", p_o_scsc_params["time_taken"], p_o_scsc_ch_params["time_taken"])
@printf("path, ours, ng-route relaxation (small N at depots/CS):        %8.3f\t%8.3f\n", p_o_ngs_params["time_taken"], p_o_ngs_ch_params["time_taken"])
@printf("path, ours, ng-route relaxation (large N at depots/CS):        %8.3f\t%8.3f\n", p_o_ngl_params["time_taken"], p_o_ngl_ch_params["time_taken"])
@printf("path, ours, ng-route relaxation (small N at depots/CS, alt):   %8.3f\t%8.3f\n", p_o_ngsa_params["time_taken"], p_o_ngsa_ch_params["time_taken"])
@printf("path, ours, ng-route relaxation (large N at depots/CS, alt):   %8.3f\t%8.3f\n", p_o_ngla_params["time_taken"], p_o_ngla_ch_params["time_taken"])

@printf("                                                               \t\tno 2-cycles\n")
@printf("subpath, benchmark:                                            %8.3f\t%8.3f\n", sp_b_params["time_taken"], sp_b_ch_params["time_taken"])
@printf("subpath, benchmark, elementary:                                %8.3f\t    ----\n", sp_b_sc_params["time_taken"])
@printf("subpath, benchmark, elementary (accel):                        %8.3f\t    ----\n", sp_b_sca_params["time_taken"])
@printf("subpath, ours:                                                 %8.3f\t%8.3f\n", sp_o_params["time_taken"], sp_o_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths:                            %8.3f\t%8.3f\n", sp_o_sc_params["time_taken"], sp_o_sc_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths (accel):                    %8.3f\t%8.3f\n", sp_o_sca_params["time_taken"], sp_o_sca_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths & paths:                    %8.3f\t%8.3f\n", sp_o_scsc_params["time_taken"], sp_o_scsc_ch_params["time_taken"])
@printf("subpath, ours, elementary subpaths & paths (accel):            %8.3f\t%8.3f\n", sp_o_scsca_params["time_taken"], sp_o_scsca_ch_params["time_taken"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS):     %8.3f\t%8.3f\n", sp_o_ngs_params["time_taken"], sp_o_ngs_ch_params["time_taken"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS):     %8.3f\t%8.3f\n", sp_o_ngl_params["time_taken"], sp_o_ngl_ch_params["time_taken"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS, alt):%8.3f\t%8.3f\n", sp_o_ngsa_params["time_taken"], sp_o_ngsa_ch_params["time_taken"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS, alt):%8.3f\t%8.3f\n", sp_o_ngla_params["time_taken"], sp_o_ngla_ch_params["time_taken"])


@printf("                                                               \t\tno 2-cycles\n")
@printf("path, benchmark:                                               %8.1f\t%8.1f\n", p_b_LP_results["objective"], p_b_ch_LP_results["objective"])
@printf("path, benchmark, elementary:                                   %8.1f\t    ----\n", p_b_sc_LP_results["objective"])
@printf("path, ours:                                                    %8.1f\t%8.1f\n", p_o_LP_results["objective"], p_o_ch_LP_results["objective"])
@printf("path, ours, elementary subpaths:                               %8.1f\t%8.1f\n", p_o_sc_LP_results["objective"], p_o_sc_ch_LP_results["objective"])
@printf("path, ours, elementary subpaths & paths:                       %8.1f\t%8.1f\n", p_o_scsc_LP_results["objective"], p_o_scsc_ch_LP_results["objective"])
@printf("path, ours, ng-route relaxation (small N at depots/CS):        %8.1f\t%8.1f\n", p_o_ngs_LP_results["objective"], p_o_ngs_ch_LP_results["objective"])
@printf("path, ours, ng-route relaxation (large N at depots/CS):        %8.1f\t%8.1f\n", p_o_ngl_LP_results["objective"], p_o_ngl_ch_LP_results["objective"])
@printf("path, ours, ng-route relaxation (small N at depots/CS, alt):   %8.1f\t%8.1f\n", p_o_ngsa_LP_results["objective"], p_o_ngsa_ch_LP_results["objective"])
@printf("path, ours, ng-route relaxation (large N at depots/CS, alt):   %8.1f\t%8.1f\n", p_o_ngla_LP_results["objective"], p_o_ngla_ch_LP_results["objective"])

@printf("                                                               \t\tno 2-cycles\n")
@printf("subpath, benchmark:                                            %8.1f\t%8.1f\n", sp_b_LP_results["objective"], sp_b_ch_LP_results["objective"])
@printf("subpath, benchmark, elementary:                                %8.1f\t    ----\n", sp_b_sc_LP_results["objective"])
@printf("subpath, benchmark, elementary (accel):                        %8.1f\t    ----\n", sp_b_sca_LP_results["objective"])
@printf("subpath, ours:                                                 %8.1f\t%8.1f\n", sp_o_LP_results["objective"], sp_o_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths:                            %8.1f\t%8.1f\n", sp_o_sc_LP_results["objective"], sp_o_sc_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths (accel):                    %8.1f\t%8.1f\n", sp_o_sca_LP_results["objective"], sp_o_sca_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths:                    %8.1f\t%8.1f\n", sp_o_scsc_LP_results["objective"], sp_o_scsc_ch_LP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths (accel):            %8.1f\t%8.1f\n", sp_o_scsca_LP_results["objective"], sp_o_scsca_ch_LP_results["objective"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS):     %8.1f\t%8.1f\n", sp_o_ngs_LP_results["objective"], sp_o_ngs_ch_LP_results["objective"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS):     %8.1f\t%8.1f\n", sp_o_ngl_LP_results["objective"], sp_o_ngl_ch_LP_results["objective"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS, alt):%8.1f\t%8.1f\n", sp_o_ngsa_LP_results["objective"], sp_o_ngsa_ch_LP_results["objective"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS, alt):%8.1f\t%8.1f\n", sp_o_ngla_LP_results["objective"], sp_o_ngla_ch_LP_results["objective"])


@printf("                                                               \t\tno 2-cycles\n")
@printf("path, benchmark:                                               %8.1f\t%8.1f\n", p_b_IP_results["objective"], p_b_ch_IP_results["objective"])
@printf("path, benchmark, elementary:                                   %8.1f\t    ----\n", p_b_sc_IP_results["objective"])
@printf("path, ours:                                                    %8.1f\t%8.1f\n", p_o_IP_results["objective"], p_o_ch_IP_results["objective"])
@printf("path, ours, elementary subpaths:                               %8.1f\t%8.1f\n", p_o_sc_IP_results["objective"], p_o_sc_ch_IP_results["objective"])
@printf("path, ours, elementary subpaths & paths:                       %8.1f\t%8.1f\n", p_o_scsc_IP_results["objective"], p_o_scsc_ch_IP_results["objective"])
@printf("path, ours, ng-route relaxation (small N at depots/CS):        %8.1f\t%8.1f\n", p_o_ngs_IP_results["objective"], p_o_ngs_ch_IP_results["objective"])
@printf("path, ours, ng-route relaxation (large N at depots/CS):        %8.1f\t%8.1f\n", p_o_ngl_IP_results["objective"], p_o_ngl_ch_IP_results["objective"])
@printf("path, ours, ng-route relaxation (small N at depots/CS, alt):   %8.1f\t%8.1f\n", p_o_ngsa_IP_results["objective"], p_o_ngsa_ch_IP_results["objective"])
@printf("path, ours, ng-route relaxation (large N at depots/CS, alt):   %8.1f\t%8.1f\n", p_o_ngla_IP_results["objective"], p_o_ngla_ch_IP_results["objective"])

@printf("                                                               \t\tno 2-cycles\n")
@printf("subpath, benchmark:                                            %8.1f\t%8.1f\n", sp_b_IP_results["objective"], sp_b_ch_IP_results["objective"])
@printf("subpath, benchmark, elementary:                                %8.1f\t    ----\n", sp_b_sc_IP_results["objective"])
@printf("subpath, benchmark, elementary (accel):                        %8.1f\t    ----\n", sp_b_sca_IP_results["objective"])
@printf("subpath, ours:                                                 %8.1f\t%8.1f\n", sp_o_IP_results["objective"], sp_o_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths:                            %8.1f\t%8.1f\n", sp_o_sc_IP_results["objective"], sp_o_sc_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths (accel):                    %8.1f\t%8.1f\n", sp_o_sca_IP_results["objective"], sp_o_sca_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths:                    %8.1f\t%8.1f\n", sp_o_scsc_IP_results["objective"], sp_o_scsc_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths (accel):            %8.1f\t%8.1f\n", sp_o_scsca_IP_results["objective"], sp_o_scsca_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS):     %8.1f\t%8.1f\n", sp_o_ngs_IP_results["objective"], sp_o_ngs_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS):     %8.1f\t%8.1f\n", sp_o_ngl_IP_results["objective"], sp_o_ngl_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS, alt):%8.1f\t%8.1f\n", sp_o_ngsa_IP_results["objective"], sp_o_ngsa_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS, alt):%8.1f\t%8.1f\n", sp_o_ngla_IP_results["objective"], sp_o_ngla_ch_IP_results["objective"])




## trial

print.(p_o_ngs_printlist);
print.(p_o_ngsa_printlist);
print.(p_o_ngs_ch_printlist);
print.(p_o_ngsa_ch_printlist);
print.(p_o_ngl_printlist);
print.(p_o_ngla_printlist);
print.(p_o_ngl_ch_printlist);
print.(p_o_ngla_ch_printlist);

print.(sp_o_ngs_printlist);
print.(sp_o_ngsa_printlist);
print.(sp_o_ngs_ch_printlist);
print.(sp_o_ngsa_ch_printlist);
print.(sp_o_ngl_printlist);
print.(sp_o_ngla_printlist);
print.(sp_o_ngl_ch_printlist);
print.(sp_o_ngla_ch_printlist);



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
method = "ours"
time_windows = false
subpath_single_service = false
subpath_check_customers = false
path_single_service = false
path_check_customers = false
christofides = true
verbose = true

### Beginning debug of path_formulation_column_generation()
compute_minimum_time_to_nearest_depot!(data, G)
compute_minimum_charge_to_nearest_depot_charging_station!(data, G)
compute_ngroute_neighborhoods!(data, Int(floor(sqrt(data["n_customers"]))))

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
        (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
            G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
            ;
            time_windows = time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
        )
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

### Starting debug of generate_base_labels_ngroute
include("utils.jl")
include("subpath_stitching.jl")
compute_ngroute_neighborhoods!(data, 4)
# compute_ngroute_neighborhoods!(data, Int(floor(sqrt(data["n_customers"]))))

base_labels = @time generate_base_labels_ngroute(G, data, κ, μ, ν, christofides = true)
full_labels = @time find_nondominated_paths_notimewindows_ngroute(data, base_labels, κ, μ, christofides = true)


base_labels = @time generate_base_labels_singleservice(G, data, κ, μ, ν, check_customers = true, christofides = true)
full_labels = @time find_nondominated_paths_notimewindows(data, base_labels, κ, μ, single_service = true, check_customers = true, christofides = true)

base_labels = @time generate_base_labels_nonsingleservice(G, data, κ, μ, ν, christofides = true)
full_labels = @time find_nondominated_paths_notimewindows(data, base_labels, κ, μ, single_service = false, check_customers = false, christofides = true)

[p.served for set in keys(full_labels[20][20]) for p in values(full_labels[20][20][set])]

κ = mp_results["κ"]
μ = mp_results["μ"]
ν = mp_results["ν"]

base_labels[19][22]
base_labels[22][21]
sort(data["neighborhoods"])

data["neighborhoods"]
collect((1,2,))

function ngroute_extend_partial_path_check(
    data,
    set::Tuple{Vararg{Int}},
    s::BaseSubpathLabel,
)
    new_set = collect(set)
    for next_node in s.nodes[2:end]
        if next_node in new_set
            return (nothing, false)
        end
        new_set = [
            node for node in new_set
                if node in data["neighborhoods"][next_node]
        ]
        push!(new_set, next_node)
        println("$next_node, $new_set")
    end
    return (Tuple(sort(unique(new_set))), true)
end

set = (8,16,22)
s = base_labels[22][21][(7,21)][end]
ngroute_extend_partial_path_check(data, set, s)

# function ngroute_customer_check(
#     data, 
#     nodes::Vector{Int}, 
#     next_node::Int,
# )
#     # Returns true if a subpath defined by nodes can be extended to customer next_node 
#     # based on the ng-route rules.
#     if !(next_node in nodes)
#         return true
#     end
#     ind = findlast(x -> x == next_node, nodes)
#     for j in nodes[ind+1:end]
#         if !(next_node in data["neighborhoods"][j])
#             return true
#         end
#     end
#     return false
# end

function ngroute_create_set(
    data, 
    set::Tuple{Vararg{Int}}, 
    next_node::Int,
)
    new_set = Int[
        node for node in set
            if node in data["neighborhoods"][next_node]
    ]
    push!(new_set, next_node) 
    return Tuple(sort(unique(new_set)))
end

function add_subpath_longlabel_to_collection!(
    collection::SortedDict{
        Int,
        BaseSubpathLabel,
    },
    k1::Int,
    v1::BaseSubpathLabel,
    ;
    verbose::Bool = false,
)
    added = true
    for (k2, v2) in pairs(collection)
        # println(k2)
        # check if v2 dominates v1
        if v2.cost ≤ v1.cost
            if all(k2 .≤ k1)
                added = false
                if verbose
                    println("$(k1), $(v1.cost) dominated by $(k2), $(v2.cost)")
                end
                break
            end
        end
        # check if v1 dominates v2
        if v1.cost ≤ v2.cost
            if all(k1 .≤ k2)
                if verbose
                    println("$(k1), $(v1.cost) dominates $(k2), $(v2.cost)")
                end
                pop!(collection, k2)
            end
        end
    end
    if added
        if verbose
            println("$(k1), $(v1.cost) added!")
        end
        insert!(collection, k1, v1)
    end
    return added
end

modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
for j in data["N_customers"]
    for i in data["N_nodes"]
        modified_costs[i,j] -= ν[j]
    end
end


base_labels = Dict(
    starting_node => Dict(
        current_node => Dict{
            Tuple{Vararg{Int}}, 
            SortedDict{Int, BaseSubpathLabel},
        }()
        for current_node in data["N_nodes"]
    )
    for starting_node in union(data["N_depots"], data["N_charging"])
)
for node in union(data["N_depots"], data["N_charging"])
    base_labels[node][node][(node,)] = SortedDict(
        0 => BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, data["n_customers"]),
        )
    )
end

unexplored_states = SortedSet(
    [
        (0.0, node, node)
        for node in union(data["N_depots"], data["N_charging"])
    ]
);


# begin
@time @suppress while length(unexplored_states) > 0
    state = pop!(unexplored_states)
    starting_node = state[end-1]
    current_node = state[end]
    for set in keys(base_labels[starting_node][current_node])
        if !(state[1] in keys(base_labels[starting_node][current_node][set]))
            println("$starting_node, $current_node, $set, $(state[1]) not found!")
            continue
        end
        current_subpath = base_labels[starting_node][current_node][set][state[1]]
        println("$starting_node, $current_node, $set, $(state[1]): $(current_subpath.nodes)")
        for next_node in setdiff(outneighbors(G, current_node), current_node)
            println("exploring next node $next_node")
            # if next_node is not a customer, proceed
            if next_node in data["N_customers"]
                # if next_node is a customer not yet visited, proceed
                # only if one can extend current_subpath along next_node according to ng-route rules
                if !ngroute_customer_check(data, current_subpath.nodes, next_node)
                    println("ngroute check failed")
                    continue
                end
                # if next_node is a customer, only proceed if 
                # it is not part of a 2-cycle (Christofides 1981)
                if christofides
                    if length(current_subpath.nodes) ≥ 2
                        if next_node == current_subpath.nodes[end-1]
                            continue
                        end
                    end
                end
            end
            # time and charge feasibility
            # if current_subpath.time_taken + data["t"][current_node, next_node] + data["min_t"][next_node] > data["T"]
            #     continue
            # end 
            if current_subpath.charge_taken + data["q"][current_node, next_node] + data["min_q"][next_node] > data["B"]
                println("not charge feasible")
                continue
            end
            new_subpath = copy(current_subpath)
            new_subpath.time_taken += data["t"][current_node, next_node]
            new_subpath.charge_taken += data["q"][current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in data["N_customers"]
                new_subpath.served[next_node] += 1
            end
            new_set = ngroute_create_set(data, set, next_node)
            println("set: $set, next_node: $next_node, new_set: $new_set")
            if !(new_set in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][new_set] = SortedDict{
                    Tuple{Vararg{Int}},
                    BaseSubpathLabel,
                }()
            end
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][new_set],
                new_subpath.time_taken, new_subpath,
                ;
                verbose = true,
            )
            println("added = $added")
            if added && next_node in data["N_customers"]
                new_state = (new_subpath.time_taken, starting_node, next_node)
                push!(unexplored_states, new_state)
                println("added next state: $new_state")
            end
        end
    end
end

unexplored_states
base_labels

for starting_node in vcat(data["N_depots"], data["N_charging"])
    for end_node in data["N_customers"]
        delete!(base_labels[starting_node], end_node)
    end
end

for starting_node in data["N_depots"]
    for end_node in vcat(data["N_depots"], data["N_charging"])
        for set in keys(base_labels[starting_node][end_node])
            for v in values(base_labels[starting_node][end_node][set])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
end
for end_node in data["N_depots"]
    for starting_node in vcat(data["N_depots"], data["N_charging"])
        for set in keys(base_labels[starting_node][end_node])
            for v in values(base_labels[starting_node][end_node][set])
                v.cost = v.cost - μ[end_node]
            end
        end
    end
end

# remove self-loops with nonnegative cost
for node in union(data["N_depots"], data["N_charging"])
    for set in keys(base_labels[node][node])
        for (k, v) in pairs(base_labels[node][node][set])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node][set], k)
            end
        end
    end
end

*("3", "3", "4")

for starting_node in union(data["N_depots"], data["N_charging"])
    vals = [
        length(keys(base_labels[starting_node][end_node]))
        for end_node in data["N_nodes"]
    ]
    # vals = [
    #     sum(
    #         [length(base_labels[starting_node][end_node][set])
    #         for set in keys(base_labels[starting_node][end_node])],
    #         init = 0,
    #     )
    #     for end_node in data["N_nodes"]
    # ]
    println(*(["$val\t" for val in vals]...))
end


base_labels[28][24]
sum(
    length(base_labels[starting_node][end_node][set])
    for starting_node in union(data["N_depots"], data["N_charging"])
        for end_node in union(data["N_depots"], data["N_charging"])
            for set in keys(base_labels[starting_node][end_node])
)
sum(
    length(base_labels[starting_node][end_node])
    for starting_node in union(data["N_depots"], data["N_charging"])
        for end_node in union(data["N_depots"], data["N_charging"])
)

plot_instance(data)

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