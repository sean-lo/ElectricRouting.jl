include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames

using Test

data, G = generate_instance(
    ;
    n_depots = 4,
    n_customers = 10,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 6,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
plot_instance(data)


a_LP_results, a_LP_params = arc_formulation(data, graph, with_charging = true, integral = false)
a_IP_results, a_IP_params = arc_formulation(data, graph, with_charging = true, integral = true, time_limit = 600)

p_b_LP_results, p_b_IP_results, p_b_params, p_b_printlist, p_b_some_paths = path_formulation_column_generation(data; method = "benchmark", verbose = true);
p_b_s_LP_results, p_b_s_IP_results, p_b_s_params, p_b_s_printlist, p_b_s_some_paths = path_formulation_column_generation(data; method = "benchmark", path_single_service = true, verbose = true)
p_b_sc_LP_results, p_b_sc_IP_results, p_b_sc_params, p_b_sc_printlist, p_b_sc_some_paths = path_formulation_column_generation(data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)
p_b_ngs_LP_results, p_b_ngs_IP_results, p_b_ngs_params, p_b_ngs_printlist, p_b_ngs_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
p_b_ngl_LP_results, p_b_ngl_IP_results, p_b_ngl_params, p_b_ngl_printlist, p_b_ngl_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)
p_b_ngsa_LP_results, p_b_ngsa_IP_results, p_b_ngsa_params, p_b_ngsa_printlist, p_b_ngsa_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
p_b_ngla_LP_results, p_b_ngla_IP_results, p_b_ngla_params, p_b_ngla_printlist, p_b_ngla_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)

p_b_ch_LP_results, p_b_ch_IP_results, p_b_ch_params, p_b_ch_printlist, p_b_ch_some_paths = path_formulation_column_generation(data; method = "benchmark", verbose = true, christofides = true)
p_b_ngs_ch_LP_results, p_b_ngs_ch_IP_results, p_b_ngs_ch_params, p_b_ngs_ch_printlist, p_b_ngs_ch_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true)
p_b_ngl_ch_LP_results, p_b_ngl_ch_IP_results, p_b_ngl_ch_params, p_b_ngl_ch_printlist, p_b_ngl_ch_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true)
p_b_ngsa_ch_LP_results, p_b_ngsa_ch_IP_results, p_b_ngsa_ch_params, p_b_ngsa_ch_printlist, p_b_ngsa_ch_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true)
p_b_ngla_ch_LP_results, p_b_ngla_ch_IP_results, p_b_ngla_ch_params, p_b_ngla_ch_printlist, p_b_ngla_ch_some_paths = path_formulation_column_generation(data; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true)

p_o_LP_results, p_o_IP_results, p_o_params, p_o_printlist, p_o_some_paths = path_formulation_column_generation(data; method = "ours", verbose = true)
p_o_s_LP_results, p_o_s_IP_results, p_o_s_params, p_o_s_printlist, p_o_s_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, verbose = true)
p_o_sc_LP_results, p_o_sc_IP_results, p_o_sc_params, p_o_sc_printlist, p_o_sc_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
p_o_ss_LP_results, p_o_ss_IP_results, p_o_ss_params, p_o_ss_printlist, p_o_ss_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
p_o_scsc_LP_results, p_o_scsc_IP_results, p_o_scsc_params, p_o_scsc_printlist, p_o_scsc_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)
p_o_ngs_LP_results, p_o_ngs_IP_results, p_o_ngs_params, p_o_ngs_printlist, p_o_ngs_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
p_o_ngl_LP_results, p_o_ngl_IP_results, p_o_ngl_params, p_o_ngl_printlist, p_o_ngl_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)
p_o_ngsa_LP_results, p_o_ngsa_IP_results, p_o_ngsa_params, p_o_ngsa_printlist, p_o_ngsa_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true)
p_o_ngla_LP_results, p_o_ngla_IP_results, p_o_ngla_params, p_o_ngla_printlist, p_o_ngla_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true)


p_o_ch_LP_results, p_o_ch_IP_results, p_o_ch_params, p_o_ch_printlist, p_o_ch_some_paths = path_formulation_column_generation(data; method = "ours", verbose = true, christofides = true)
p_o_s_ch_LP_results, p_o_s_ch_IP_results, p_o_s_ch_params, p_o_s_ch_printlist, p_o_s_ch_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, verbose = true, christofides = true)
p_o_sc_ch_LP_results, p_o_sc_ch_IP_results, p_o_sc_ch_params, p_o_sc_ch_printlist, p_o_sc_ch_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true, christofides = true)
p_o_ss_ch_LP_results, p_o_ss_ch_IP_results, p_o_ss_ch_params, p_o_ss_ch_printlist, p_o_ss_ch_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true, christofides = true)
p_o_scsc_ch_LP_results, p_o_scsc_ch_IP_results, p_o_scsc_ch_params, p_o_scsc_ch_printlist, p_o_scsc_ch_some_paths = path_formulation_column_generation(data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true, christofides = true)
p_o_ngs_ch_LP_results, p_o_ngs_ch_IP_results, p_o_ngs_ch_params, p_o_ngs_ch_printlist, p_o_ngs_ch_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true)
p_o_ngl_ch_LP_results, p_o_ngl_ch_IP_results, p_o_ngl_ch_params, p_o_ngl_ch_printlist, p_o_ngl_ch_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true)
p_o_ngsa_ch_LP_results, p_o_ngsa_ch_IP_results, p_o_ngsa_ch_params, p_o_ngsa_ch_printlist, p_o_ngsa_ch_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true)
p_o_ngla_ch_LP_results, p_o_ngla_ch_IP_results, p_o_ngla_ch_params, p_o_ngla_ch_printlist, p_o_ngla_ch_some_paths = path_formulation_column_generation(data; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true)

sp_b_LP_results, sp_b_IP_results, sp_b_params, sp_b_printlist, sp_b_some_subpaths, sp_b_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", verbose = true);
sp_b_s_LP_results, sp_b_s_IP_results, sp_b_s_params, sp_b_s_printlist, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", path_single_service = true, verbose = true);
sp_b_sc_LP_results, sp_b_sc_IP_results, sp_b_sc_params, sp_b_sc_printlist, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true);
sp_b_sca_LP_results, sp_b_sca_IP_results, sp_b_sca_params, sp_b_sca_printlist, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true);
sp_b_ngs_LP_results, sp_b_ngs_IP_results, sp_b_ngs_params, sp_b_ngs_printlist, sp_b_ngs_some_subpaths, sp_b_ngs_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true);
sp_b_ngl_LP_results, sp_b_ngl_IP_results, sp_b_ngl_params, sp_b_ngl_printlist, sp_b_ngl_some_subpaths, sp_b_ngl_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true);
sp_b_ngsa_LP_results, sp_b_ngsa_IP_results, sp_b_ngsa_params, sp_b_ngsa_printlist, sp_b_ngsa_some_subpaths, sp_b_ngsa_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true);
sp_b_ngla_LP_results, sp_b_ngla_IP_results, sp_b_ngla_params, sp_b_ngla_printlist, sp_b_ngla_some_subpaths, sp_b_ngla_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true);

sp_b_ch_LP_results, sp_b_ch_IP_results, sp_b_ch_params, sp_b_ch_printlist, sp_b_ch_some_subpaths, sp_b_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", verbose = true, christofides = true);
sp_b_ngs_ch_LP_results, sp_b_ngs_ch_IP_results, sp_b_ngs_ch_params, sp_b_ngs_ch_printlist, sp_b_ngs_ch_some_subpaths, sp_b_ngs_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
sp_b_ngl_ch_LP_results, sp_b_ngl_ch_IP_results, sp_b_ngl_ch_params, sp_b_ngl_ch_printlist, sp_b_ngl_ch_some_subpaths, sp_b_ngl_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true);
sp_b_ngsa_ch_LP_results, sp_b_ngsa_ch_IP_results, sp_b_ngsa_ch_params, sp_b_ngsa_ch_printlist, sp_b_ngsa_ch_some_subpaths, sp_b_ngsa_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true, christofides = true);
sp_b_ngla_ch_LP_results, sp_b_ngla_ch_IP_results, sp_b_ngla_ch_params, sp_b_ngla_ch_printlist, sp_b_ngla_ch_some_subpaths, sp_b_ngla_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "benchmark", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true, christofides = true);


sp_o_LP_results, sp_o_IP_results, sp_o_params, sp_o_printlist, sp_o_some_subpaths, sp_o_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", verbose = true);
sp_o_s_LP_results, sp_o_s_IP_results, sp_o_s_params, sp_o_s_printlist, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, verbose = true);
sp_o_sc_LP_results, sp_o_sc_IP_results, sp_o_sc_params, sp_o_sc_printlist, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true);
sp_o_sca_LP_results, sp_o_sca_IP_results, sp_o_sca_params, sp_o_sca_printlist, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true);
sp_o_ss_LP_results, sp_o_ss_IP_results, sp_o_ss_params, sp_o_ss_printlist, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true);
sp_o_scsc_LP_results, sp_o_scsc_IP_results, sp_o_scsc_params, sp_o_scsc_printlist, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true);
sp_o_scsca_LP_results, sp_o_scsca_IP_results, sp_o_scsca_params, sp_o_scsca_printlist, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true);
sp_o_ngs_LP_results, sp_o_ngs_IP_results, sp_o_ngs_params, sp_o_ngs_printlist, sp_o_ngs_some_subpaths, sp_o_ngs_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true);
sp_o_ngl_LP_results, sp_o_ngl_IP_results, sp_o_ngl_params, sp_o_ngl_printlist, sp_o_ngl_some_subpaths, sp_o_ngl_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true);
sp_o_ngsa_LP_results, sp_o_ngsa_IP_results, sp_o_ngsa_params, sp_o_ngsa_printlist, sp_o_ngsa_some_subpaths, sp_o_ngsa_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", verbose = true);
sp_o_ngla_LP_results, sp_o_ngla_IP_results, sp_o_ngla_params, sp_o_ngla_printlist, sp_o_ngla_some_subpaths, sp_o_ngla_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", verbose = true);

sp_o_ch_LP_results, sp_o_ch_IP_results, sp_o_ch_params, sp_o_ch_printlist, sp_o_ch_some_subpaths, sp_o_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", verbose = true, christofides = true);
sp_o_s_ch_LP_results, sp_o_s_ch_IP_results, sp_o_s_ch_params, sp_o_s_ch_printlist, sp_o_s_ch_some_subpaths, sp_o_s_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, verbose = true, christofides = true);
sp_o_sc_ch_LP_results, sp_o_sc_ch_IP_results, sp_o_sc_ch_params, sp_o_sc_ch_printlist, sp_o_sc_ch_some_subpaths, sp_o_sc_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true, christofides = true);
sp_o_sca_ch_LP_results, sp_o_sca_ch_IP_results, sp_o_sca_ch_params, sp_o_sca_ch_printlist, sp_o_sca_ch_some_subpaths, sp_o_sca_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true, christofides = true);
sp_o_ss_ch_LP_results, sp_o_ss_ch_IP_results, sp_o_ss_ch_params, sp_o_ss_ch_printlist, sp_o_ss_ch_some_subpaths, sp_o_ss_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true, christofides = true);
sp_o_scsc_ch_LP_results, sp_o_scsc_ch_IP_results, sp_o_scsc_ch_params, sp_o_scsc_ch_printlist, sp_o_scsc_ch_some_subpaths, sp_o_scsc_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true, christofides = true);
sp_o_scsca_ch_LP_results, sp_o_scsca_ch_IP_results, sp_o_scsca_ch_params, sp_o_scsca_ch_printlist, sp_o_scsca_ch_some_subpaths, sp_o_scsca_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true, christofides = true);
sp_o_ngs_ch_LP_results, sp_o_ngs_ch_IP_results, sp_o_ngs_ch_params, sp_o_ngs_ch_printlist, sp_o_ngs_ch_some_subpaths, sp_o_ngs_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", christofides = true, verbose = true);
sp_o_ngl_ch_LP_results, sp_o_ngl_ch_IP_results, sp_o_ngl_ch_params, sp_o_ngl_ch_printlist, sp_o_ngl_ch_some_subpaths, sp_o_ngl_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", christofides = true, verbose = true);
sp_o_ngsa_ch_LP_results, sp_o_ngsa_ch_IP_results, sp_o_ngsa_ch_params, sp_o_ngsa_ch_printlist, sp_o_ngsa_ch_some_subpaths, sp_o_ngsa_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "small", christofides = true, verbose = true);
sp_o_ngla_ch_LP_results, sp_o_ngla_ch_IP_results, sp_o_ngla_ch_params, sp_o_ngla_ch_printlist, sp_o_ngla_ch_some_subpaths, sp_o_ngla_ch_some_charging_arcs = @time subpath_formulation_column_generation_integrated_from_paths(data, G; method = "ours", ngroute = true, ngroute_neighborhood_charging_depots_size = "large", christofides = true, verbose = true);


collect_path_solution_metrics!(p_b_LP_results, data, graph, p_b_some_paths)
collect_path_solution_metrics!(p_b_s_LP_results, data, graph, p_b_s_some_paths)
collect_path_solution_metrics!(p_b_sc_LP_results, data, graph, p_b_sc_some_paths)
collect_path_solution_metrics!(p_b_ngs_LP_results, data, graph, p_b_ngs_some_paths)
collect_path_solution_metrics!(p_b_ngsa_LP_results, data, graph, p_b_ngsa_some_paths)
collect_path_solution_metrics!(p_b_ngl_LP_results, data, graph, p_b_ngl_some_paths)
collect_path_solution_metrics!(p_b_ngla_LP_results, data, graph, p_b_ngla_some_paths)
collect_path_solution_metrics!(p_b_ch_LP_results, data, graph, p_b_ch_some_paths)
collect_path_solution_metrics!(p_b_ngs_ch_LP_results, data, graph, p_b_ngs_ch_some_paths)
collect_path_solution_metrics!(p_b_ngsa_ch_LP_results, data, graph, p_b_ngsa_ch_some_paths)
collect_path_solution_metrics!(p_b_ngl_ch_LP_results, data, graph, p_b_ngl_ch_some_paths)
collect_path_solution_metrics!(p_b_ngla_ch_LP_results, data, graph, p_b_ngla_ch_some_paths)

collect_path_solution_metrics!(p_o_LP_results, data, graph, p_o_some_paths)
collect_path_solution_metrics!(p_o_s_LP_results, data, graph, p_o_s_some_paths)
collect_path_solution_metrics!(p_o_sc_LP_results, data, graph, p_o_sc_some_paths)
collect_path_solution_metrics!(p_o_ss_LP_results, data, graph, p_o_ss_some_paths)
collect_path_solution_metrics!(p_o_scsc_LP_results, data, graph, p_o_scsc_some_paths)
collect_path_solution_metrics!(p_o_ngs_LP_results, data, graph, p_o_ngs_some_paths)
collect_path_solution_metrics!(p_o_ngl_LP_results, data, graph, p_o_ngl_some_paths)
collect_path_solution_metrics!(p_o_ngsa_LP_results, data, graph, p_o_ngsa_some_paths)
collect_path_solution_metrics!(p_o_ngla_LP_results, data, graph, p_o_ngla_some_paths)

collect_path_solution_metrics!(p_o_ch_LP_results, data, graph, p_o_ch_some_paths)
collect_path_solution_metrics!(p_o_s_ch_LP_results, data, graph, p_o_s_ch_some_paths)
collect_path_solution_metrics!(p_o_sc_ch_LP_results, data, graph, p_o_sc_ch_some_paths)
collect_path_solution_metrics!(p_o_ss_ch_LP_results, data, graph, p_o_ss_ch_some_paths)
collect_path_solution_metrics!(p_o_scsc_ch_LP_results, data, graph, p_o_scsc_ch_some_paths)
collect_path_solution_metrics!(p_o_ngs_ch_LP_results, data, graph, p_o_ngs_ch_some_paths)
collect_path_solution_metrics!(p_o_ngl_ch_LP_results, data, graph, p_o_ngl_ch_some_paths)
collect_path_solution_metrics!(p_o_ngsa_ch_LP_results, data, graph, p_o_ngsa_ch_some_paths)
collect_path_solution_metrics!(p_o_ngla_ch_LP_results, data, graph, p_o_ngla_ch_some_paths)

collect_subpath_solution_metrics!(sp_b_LP_results, data, sp_b_some_subpaths, sp_b_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_s_LP_results, data, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sc_LP_results, data, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sca_LP_results, data, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs)

collect_subpath_solution_metrics!(sp_b_ngs_LP_results, data, sp_b_ngs_some_subpaths, sp_b_ngs_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngl_LP_results, data, sp_b_ngl_some_subpaths, sp_b_ngl_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngsa_LP_results, data, sp_b_ngsa_some_subpaths, sp_b_ngsa_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngla_LP_results, data, sp_b_ngla_some_subpaths, sp_b_ngla_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ch_LP_results, data, sp_b_ch_some_subpaths, sp_b_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngs_ch_LP_results, data, sp_b_ngs_ch_some_subpaths, sp_b_ngs_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngl_ch_LP_results, data, sp_b_ngl_ch_some_subpaths, sp_b_ngl_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngsa_ch_LP_results, data, sp_b_ngsa_ch_some_subpaths, sp_b_ngsa_ch_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_ngla_ch_LP_results, data, sp_b_ngla_ch_some_subpaths, sp_b_ngla_ch_some_charging_arcs)

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

a_LP_results["objective"]

p_b_LP_results["objective"]
p_b_s_LP_results["objective"]
p_b_sc_LP_results["objective"]
p_b_ngs_LP_results["objective"]
p_b_ngl_LP_results["objective"]
p_b_ngsa_LP_results["objective"]
p_b_ngla_LP_results["objective"]

p_b_ch_LP_results["objective"]
p_b_ngs_ch_LP_results["objective"]
p_b_ngl_ch_LP_results["objective"]
p_b_ngsa_ch_LP_results["objective"]
p_b_ngla_ch_LP_results["objective"]

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
sp_b_ngs_LP_results["objective"]
sp_b_ngl_LP_results["objective"]
sp_b_ngsa_LP_results["objective"]
sp_b_ngla_LP_results["objective"]

sp_b_ch_LP_results["objective"]
sp_b_ngs_ch_LP_results["objective"]
sp_b_ngl_ch_LP_results["objective"]
sp_b_ngsa_ch_LP_results["objective"]
sp_b_ngla_ch_LP_results["objective"]

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
    p_b_LP_results["objective"] 
    ≤ p_b_ngs_LP_results["objective"] 
    ≤ p_b_ngl_LP_results["objective"] 
    ≤ p_b_sc_LP_results["objective"]
)
@test p_b_ngs_LP_results["objective"] ≈ p_b_ngsa_LP_results["objective"]
@test p_b_ngl_LP_results["objective"] ≈ p_b_ngla_LP_results["objective"]
@test(
    p_b_ch_LP_results["objective"] 
    ≤ p_b_ngs_ch_LP_results["objective"]
    ≤ p_b_ngl_ch_LP_results["objective"]
    ≤ p_b_scsc_ch_LP_results["objective"]
)
@test p_b_ngs_ch_LP_results["objective"] ≈ p_b_ngsa_ch_LP_results["objective"]
@test p_b_ngl_ch_LP_results["objective"] ≈ p_b_ngla_ch_LP_results["objective"]

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
    sp_b_LP_results["objective"] 
    ≤ sp_b_ngs_LP_results["objective"] 
    ≤ sp_b_ngl_LP_results["objective"] 
    ≤ sp_b_sc_LP_results["objective"]
)
@test sp_b_ngs_LP_results["objective"] ≈ sp_b_ngsa_LP_results["objective"]
@test sp_b_ngl_LP_results["objective"] ≈ sp_b_ngla_LP_results["objective"]
@test(
    sp_b_ch_LP_results["objective"] 
    ≤ sp_b_ngs_ch_LP_results["objective"]
    ≤ sp_b_ngl_ch_LP_results["objective"]
    ≤ sp_b_scsc_ch_LP_results["objective"]
)
@test sp_b_ngs_ch_LP_results["objective"] ≈ sp_b_ngsa_ch_LP_results["objective"]
@test sp_b_ngl_ch_LP_results["objective"] ≈ sp_b_ngla_ch_LP_results["objective"]

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
p_b_ngs_params["time_taken"]
p_b_ngl_params["time_taken"]
p_b_ngsa_params["time_taken"]
p_b_ngla_params["time_taken"]

p_b_ch_params["time_taken"]
p_b_ngs_ch_params["time_taken"]
p_b_ngl_ch_params["time_taken"]
p_b_ngsa_ch_params["time_taken"]
p_b_ngla_ch_params["time_taken"]

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
sp_b_ngs_params["time_taken"]
sp_b_ngl_params["time_taken"]
sp_b_ngsa_params["time_taken"]
sp_b_ngla_params["time_taken"]

sp_b_ch_params["time_taken"]
sp_b_ngs_ch_params["time_taken"]
sp_b_ngl_ch_params["time_taken"]
sp_b_ngsa_ch_params["time_taken"]
sp_b_ngla_ch_params["time_taken"]

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
@printf("path, benc, ng-route relaxation (small N at depots/CS):        %8.3f\t%8.3f\n", p_b_ngs_params["time_taken"], p_b_ngs_ch_params["time_taken"])
@printf("path, benc, ng-route relaxation (large N at depots/CS):        %8.3f\t%8.3f\n", p_b_ngl_params["time_taken"], p_b_ngl_ch_params["time_taken"])
@printf("path, benc, ng-route relaxation (small N at depots/CS, alt):   %8.3f\t%8.3f\n", p_b_ngsa_params["time_taken"], p_b_ngsa_ch_params["time_taken"])
@printf("path, benc, ng-route relaxation (large N at depots/CS, alt):   %8.3f\t%8.3f\n", p_b_ngla_params["time_taken"], p_b_ngla_ch_params["time_taken"])
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
@printf("subpath, benc, ng-route relaxation (small N at depots/CS):     %8.3f\t%8.3f\n", sp_b_ngs_params["time_taken"], sp_b_ngs_ch_params["time_taken"])
@printf("subpath, benc, ng-route relaxation (large N at depots/CS):     %8.3f\t%8.3f\n", sp_b_ngl_params["time_taken"], sp_b_ngl_ch_params["time_taken"])
@printf("subpath, benc, ng-route relaxation (small N at depots/CS, alt):%8.3f\t%8.3f\n", sp_b_ngsa_params["time_taken"], sp_b_ngsa_ch_params["time_taken"])
@printf("subpath, benc, ng-route relaxation (large N at depots/CS, alt):%8.3f\t%8.3f\n", sp_b_ngla_params["time_taken"], sp_b_ngla_ch_params["time_taken"])
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
@printf("arc formulation                                                %8.1f\t    ----\n", a_LP_results["objective"])
@printf("path, benchmark:                                               %8.1f\t%8.1f\n", p_b_LP_results["objective"], p_b_ch_LP_results["objective"])
@printf("path, benchmark, elementary:                                   %8.1f\t    ----\n", p_b_sc_LP_results["objective"])
@printf("path, benc, ng-route relaxation (small N at depots/CS):        %8.1f\t%8.1f\n", p_b_ngs_LP_results["objective"], p_b_ngs_ch_LP_results["objective"])
@printf("path, benc, ng-route relaxation (large N at depots/CS):        %8.1f\t%8.1f\n", p_b_ngl_LP_results["objective"], p_b_ngl_ch_LP_results["objective"])
@printf("path, benc, ng-route relaxation (small N at depots/CS, alt):   %8.1f\t%8.1f\n", p_b_ngsa_LP_results["objective"], p_b_ngsa_ch_LP_results["objective"])
@printf("path, benc, ng-route relaxation (large N at depots/CS, alt):   %8.1f\t%8.1f\n", p_b_ngla_LP_results["objective"], p_b_ngla_ch_LP_results["objective"])
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
@printf("subpath, benc, ng-route relaxation (small N at depots/CS):     %8.1f\t%8.1f\n", sp_b_ngs_LP_results["objective"], sp_b_ngs_ch_LP_results["objective"])
@printf("subpath, benc, ng-route relaxation (large N at depots/CS):     %8.1f\t%8.1f\n", sp_b_ngl_LP_results["objective"], sp_b_ngl_ch_LP_results["objective"])
@printf("subpath, benc, ng-route relaxation (small N at depots/CS, alt):%8.1f\t%8.1f\n", sp_b_ngsa_LP_results["objective"], sp_b_ngsa_ch_LP_results["objective"])
@printf("subpath, benc, ng-route relaxation (large N at depots/CS, alt):%8.1f\t%8.1f\n", sp_b_ngla_LP_results["objective"], sp_b_ngla_ch_LP_results["objective"])
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
@printf("arc formulation                                                %8.1f\t    ----\n", a_IP_results["objective"])
@printf("path, benchmark:                                               %8.1f\t%8.1f\n", p_b_IP_results["objective"], p_b_ch_IP_results["objective"])
@printf("path, benchmark, elementary:                                   %8.1f\t    ----\n", p_b_sc_IP_results["objective"])
@printf("path, benc, ng-route relaxation (small N at depots/CS):        %8.1f\t%8.1f\n", p_b_ngs_IP_results["objective"], p_b_ngs_ch_IP_results["objective"])
@printf("path, benc, ng-route relaxation (large N at depots/CS):        %8.1f\t%8.1f\n", p_b_ngl_IP_results["objective"], p_b_ngl_ch_IP_results["objective"])
@printf("path, benc, ng-route relaxation (small N at depots/CS, alt):   %8.1f\t%8.1f\n", p_b_ngsa_IP_results["objective"], p_b_ngsa_ch_IP_results["objective"])
@printf("path, benc, ng-route relaxation (large N at depots/CS, alt):   %8.1f\t%8.1f\n", p_b_ngla_IP_results["objective"], p_b_ngla_ch_IP_results["objective"])
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
@printf("subpath, benc, ng-route relaxation (small N at depots/CS):     %8.1f\t%8.1f\n", sp_b_ngs_IP_results["objective"], sp_b_ngs_ch_IP_results["objective"])
@printf("subpath, benc, ng-route relaxation (large N at depots/CS):     %8.1f\t%8.1f\n", sp_b_ngl_IP_results["objective"], sp_b_ngl_ch_IP_results["objective"])
@printf("subpath, benc, ng-route relaxation (small N at depots/CS, alt):%8.1f\t%8.1f\n", sp_b_ngsa_IP_results["objective"], sp_b_ngsa_ch_IP_results["objective"])
@printf("subpath, benc, ng-route relaxation (large N at depots/CS, alt):%8.1f\t%8.1f\n", sp_b_ngla_IP_results["objective"], sp_b_ngla_ch_IP_results["objective"])
@printf("subpath, ours:                                                 %8.1f\t%8.1f\n", sp_o_IP_results["objective"], sp_o_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths:                            %8.1f\t%8.1f\n", sp_o_sc_IP_results["objective"], sp_o_sc_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths (accel):                    %8.1f\t%8.1f\n", sp_o_sca_IP_results["objective"], sp_o_sca_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths:                    %8.1f\t%8.1f\n", sp_o_scsc_IP_results["objective"], sp_o_scsc_ch_IP_results["objective"])
@printf("subpath, ours, elementary subpaths & paths (accel):            %8.1f\t%8.1f\n", sp_o_scsca_IP_results["objective"], sp_o_scsca_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS):     %8.1f\t%8.1f\n", sp_o_ngs_IP_results["objective"], sp_o_ngs_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS):     %8.1f\t%8.1f\n", sp_o_ngl_IP_results["objective"], sp_o_ngl_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (small N at depots/CS, alt):%8.1f\t%8.1f\n", sp_o_ngsa_IP_results["objective"], sp_o_ngsa_ch_IP_results["objective"])
@printf("subpath, ours, ng-route relaxation (large N at depots/CS, alt):%8.1f\t%8.1f\n", sp_o_ngla_IP_results["objective"], sp_o_ngla_ch_IP_results["objective"])



## Plotting

include("subpath_formulation.jl")
plot_subpath_solution(
    sp_o_ngs_ch_LP_results, 
    data, 
    sp_o_ngs_ch_some_subpaths, 
    sp_o_ngs_ch_some_charging_arcs,
)
sp_o_ngs_ch_LP_results["paths"]
## trial

print.(p_b_ngs_printlist);
print.(p_b_ngsa_printlist);
print.(p_b_ngs_ch_printlist);
print.(p_b_ngsa_ch_printlist);
print.(p_b_ngl_printlist);
print.(p_b_ngla_printlist);
print.(p_b_ngl_ch_printlist);
print.(p_b_ngla_ch_printlist);

print.(sp_b_ngs_printlist);
print.(sp_b_ngsa_printlist);
print.(sp_b_ngs_ch_printlist);
print.(sp_b_ngsa_ch_printlist);
print.(sp_b_ngl_printlist);
print.(sp_b_ngla_printlist);
print.(sp_b_ngl_ch_printlist);
print.(sp_b_ngla_ch_printlist);

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

### Testing

include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames

using Test

all_results = []
for (n_customers, n_vehicles, T, n_charging) in Iterators.product(
    20:4:40,
    [6,8],
    40000:10000:80000,
    [13],
)
    data, G = generate_instance(
        ;
        n_depots = 4,
        n_customers = n_customers,
        n_charging = n_charging,
        n_vehicles = n_vehicles,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "circular_packing",
        shrinkage_depots = 1.0,
        shrinkage_charging = 1.0,
        T = T,
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
    # plot_instance(data)
    LP_results, IP_results, cgparams, printlist, some_subpaths, some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(
        data, G;
        method = "ours",
        time_windows = false,
        ngroute = true,
        christofides = true,
    )
    collect_subpath_solution_metrics!(LP_results, data, some_subpaths, some_charging_arcs)
    collect_subpath_solution_metrics!(IP_results, data, some_subpaths, some_charging_arcs)
    push!(all_results, (
        n_customers = n_customers,
        n_vehicles = n_vehicles,
        T = T,
        CG_time_taken = cgparams["time_taken"],
        CG_lp_relaxation_time_taken_mean = cgparams["lp_relaxation_time_taken_mean"],
        CG_sp_base_time_taken_mean = cgparams["sp_base_time_taken_mean"],
        CG_sp_full_time_taken_mean = cgparams["sp_full_time_taken_mean"],
        CG_sp_time_taken_mean = cgparams["sp_time_taken_mean"],
        LP_objective = LP_results["objective"],
        LP_mean_subpath_length = LP_results["mean_subpath_length"],
        LP_mean_subpath_ncust = LP_results["mean_subpath_ncust"],
        LP_mean_path_length = LP_results["mean_path_length"],
        LP_mean_path_ncust = LP_results["mean_path_ncust"],
        LP_mean_ps_length = LP_results["mean_ps_length"],
        LP_weighted_mean_subpath_length = LP_results["weighted_mean_subpath_length"],
        LP_weighted_mean_subpath_ncust = LP_results["weighted_mean_subpath_ncust"],
        LP_weighted_mean_path_length = LP_results["weighted_mean_path_length"],
        LP_weighted_mean_path_ncust = LP_results["weighted_mean_path_ncust"],
        LP_weighted_mean_ps_length = LP_results["weighted_mean_ps_length"],
        LP_utilization = LP_results["utilization"],
        LP_driving_time_proportion = LP_results["driving_time_proportion"],
        LP_charging_time_proportion = LP_results["charging_time_proportion"],
        IP_objective = IP_results["objective"],
        IP_mean_subpath_length = IP_results["mean_subpath_length"],
        IP_mean_subpath_ncust = IP_results["mean_subpath_ncust"],
        IP_mean_path_length = IP_results["mean_path_length"],
        IP_mean_path_ncust = IP_results["mean_path_ncust"],
        IP_mean_ps_length = IP_results["mean_ps_length"],
        IP_weighted_mean_subpath_length = IP_results["weighted_mean_subpath_length"],
        IP_weighted_mean_subpath_ncust = IP_results["weighted_mean_subpath_ncust"],
        IP_weighted_mean_path_length = IP_results["weighted_mean_path_length"],
        IP_weighted_mean_path_ncust = IP_results["weighted_mean_path_ncust"],
        IP_weighted_mean_ps_length = IP_results["weighted_mean_ps_length"],
        IP_utilization = IP_results["utilization"],
        IP_driving_time_proportion = IP_results["driving_time_proportion"],
        IP_charging_time_proportion = IP_results["charging_time_proportion"],
    ))
end
df = DataFrame(all_results) 
df |>
    x -> sort(x, [:n_vehicles, :n_customers, :T]) |>
    x -> select(x, [
        :n_customers, 
        :n_vehicles,
        :T,
        :CG_time_taken,
        :CG_lp_relaxation_time_taken_mean,
        :CG_sp_base_time_taken_mean,
        :CG_sp_full_time_taken_mean,
        :CG_sp_time_taken_mean,
        :LP_objective, :IP_objective,
        :LP_weighted_mean_subpath_length,
        :LP_weighted_mean_path_length,
        :LP_weighted_mean_ps_length,
    ]) |>
    x -> show(x, allrows = true)

df |>
    x -> filter(r -> r.n_vehicles == 6, df) |>
    x -> unstack(x, :n_customers, :T, :LP_weighted_mean_subpath_length)
df |>
    x -> filter(r -> r.n_vehicles == 6, df) |>
    x -> unstack(x, :n_customers, :T, :LP_weighted_mean_ps_length)
df |>
    x -> filter(r -> r.n_vehicles == 6, df) |>
    x -> unstack(x, :n_customers, :T, :LP_weighted_mean_path_length)



time_feas_results = []
for (n_customers, n_vehicles, T, n_charging) in Iterators.product(
    [16],
    [8],
    [40000],
    [7],
)
    data, G = generate_instance(
        ;
        n_depots = 4,
        n_customers = n_customers,
        n_charging = n_charging,
        n_vehicles = n_vehicles,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "circular_packing",
        shrinkage_depots = 1.0,
        shrinkage_charging = 1.0,
        T = T,
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
    # plot_instance(data)
    LP_results, IP_results, cgparams, printlist, some_paths = path_formulation_column_generation(data;
        method = "benchmark",
        path_single_service = true,
        path_check_customers = true,
        time_windows = false,
        christofides = false,
    )
    collect_path_solution_metrics!(LP_results, data, graph, some_paths)
    collect_path_solution_metrics!(IP_results, data, graph, some_paths)
    push!(time_feas_results, (
        n_customers = n_customers,
        n_vehicles = n_vehicles,
        T = T,
        CG_time_taken = cgparams["time_taken"],
        CG_lp_relaxation_time_taken_mean = cgparams["lp_relaxation_time_taken_mean"],
        CG_sp_base_time_taken_mean = cgparams["sp_base_time_taken_mean"],
        CG_sp_full_time_taken_mean = cgparams["sp_full_time_taken_mean"],
        CG_sp_time_taken_mean = cgparams["sp_time_taken_mean"],
        LP_objective = LP_results["objective"],
        LP_mean_subpath_length = LP_results["mean_subpath_length"],
        LP_mean_subpath_ncust = LP_results["mean_subpath_ncust"],
        LP_mean_path_length = LP_results["mean_path_length"],
        LP_mean_path_ncust = LP_results["mean_path_ncust"],
        LP_mean_ps_length = LP_results["mean_ps_length"],
        LP_weighted_mean_subpath_length = LP_results["weighted_mean_subpath_length"],
        LP_weighted_mean_subpath_ncust = LP_results["weighted_mean_subpath_ncust"],
        LP_weighted_mean_path_length = LP_results["weighted_mean_path_length"],
        LP_weighted_mean_path_ncust = LP_results["weighted_mean_path_ncust"],
        LP_weighted_mean_ps_length = LP_results["weighted_mean_ps_length"],
        LP_utilization = LP_results["utilization"],
        LP_driving_time_proportion = LP_results["driving_time_proportion"],
        LP_charging_time_proportion = LP_results["charging_time_proportion"],
        IP_objective = IP_results["objective"],
        IP_mean_subpath_length = IP_results["mean_subpath_length"],
        IP_mean_subpath_ncust = IP_results["mean_subpath_ncust"],
        IP_mean_path_length = IP_results["mean_path_length"],
        IP_mean_path_ncust = IP_results["mean_path_ncust"],
        IP_mean_ps_length = IP_results["mean_ps_length"],
        IP_weighted_mean_subpath_length = IP_results["weighted_mean_subpath_length"],
        IP_weighted_mean_subpath_ncust = IP_results["weighted_mean_subpath_ncust"],
        IP_weighted_mean_path_length = IP_results["weighted_mean_path_length"],
        IP_weighted_mean_path_ncust = IP_results["weighted_mean_path_ncust"],
        IP_weighted_mean_ps_length = IP_results["weighted_mean_ps_length"],
        IP_utilization = IP_results["utilization"],
        IP_driving_time_proportion = IP_results["driving_time_proportion"],
        IP_charging_time_proportion = IP_results["charging_time_proportion"],
    ))
end
df = DataFrame(all_results) 
df |>
    x -> sort(x, [:n_vehicles, :n_customers, :T]) |>
    x -> select(x, [
        :n_customers, 
        :n_vehicles,
        :T,
        :CG_time_taken,
        :CG_lp_relaxation_time_taken_mean,
        :CG_sp_base_time_taken_mean,
        :CG_sp_full_time_taken_mean,
        :CG_sp_time_taken_mean,
        :LP_objective, :IP_objective,
        :LP_weighted_mean_subpath_length,
        :LP_weighted_mean_path_length,
        :LP_weighted_mean_ps_length,
    ]) |>
    x -> show(x, allrows = true)

df |>
    x -> filter(r -> r.n_vehicles == 6, df) |>
    x -> unstack(x, :n_customers, :T, :LP_weighted_mean_subpath_length)
df |>
    x -> filter(r -> r.n_vehicles == 6, df) |>
    x -> unstack(x, :n_customers, :T, :LP_weighted_mean_ps_length)
df |>
    x -> filter(r -> r.n_vehicles == 6, df) |>
    x -> unstack(x, :n_customers, :T, :LP_weighted_mean_path_length)

include("path_formulation.jl")
include("subpath_formulation.jl")
include("utils.jl")
TIME_LIMIT_SEC = 3.0
# simple test case to quickly compile 

    sample_data, sample_G = generate_instance(
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
        # data_dir = "../../../data/",
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
        # ngroute_neighborhood_charging_depots_size
        ("path", "benchmark", false, false, false, false, false,  true, false, "none"),
        ("path", "benchmark", false, false,  true,  true, false, false, false, "none"),
        ("path", "benchmark", false, false, false, false, false,  true,  true, "small"),
        ("path", "benchmark", false, false, false, false, false,  true,  true, "large"),
        ("path", "ours", false, false, false, false, false,  true, false, "none"),
        ("path", "ours",  true,  true, false, false, false,  true, false, "none"),
        ("path", "ours",  true,  true,  true,  true, false,  true, false, "none"),
        ("path", "ours", false, false, false, false, false,  true,  true,  "small"),
        ("path", "ours", false, false, false, false, false,  true,  true,  "large"),
        ("subpath", "benchmark", false, false, false, false, false,  true, false, "none"),
        ("subpath", "benchmark", false, false,  true,  true, false, false, false, "none"),
        ("subpath", "benchmark", false, false,  true,  true,  true, false, false, "none"),
        ("subpath", "benchmark", false, false, false, false, false,  true,  true, "small"),
        ("subpath", "benchmark", false, false, false, false, false,  true,  true, "large"),
        ("subpath", "ours", false, false, false, false, false,  true, false, "none"),
        ("subpath", "ours",  true,  true, false, false, false,  true, false, "none"),
        ("subpath", "ours",  true,  true, false, false,  true,  true, false, "none"),
        ("subpath", "ours",  true,  true,  true,  true, false,  true, false, "none"),
        ("subpath", "ours",  true,  true,  true,  true,  true,  true, false, "none"),
        ("subpath", "ours", false, false, false, false, false,  true,  true, "small"),
        ("subpath", "ours", false, false, false, false, false,  true,  true, "large"),
    ]
    for method_param in method_params
        if method_param[1] == "path"
            (
                LP_results, IP_results, cgparams, printlist, paths
            ) = path_formulation_column_generation(
                sample_data,
                ;
                method = method_param[2],
                subpath_single_service = method_param[3],
                subpath_check_customers = method_param[4],
                path_single_service = method_param[5],
                path_check_customers = method_param[6],
                christofides = method_param[8],
                ngroute = method_param[9],
                ngroute_neighborhood_charging_depots_size = method_param[11],
                time_limit = TIME_LIMIT_SEC,
            )
        elseif method_param[1] == "subpath"
            (
                LP_results, IP_results, cgparams, printlist, subpaths, charging_arcs
            ) = subpath_formulation_column_generation_integrated_from_paths(
                sample_data, sample_G,
                ;
                method = method_param[2],
                subpath_single_service = method_param[3],
                subpath_check_customers = method_param[4],
                path_single_service = method_param[5],
                path_check_customers = method_param[6],
                christofides = method_param[8],
                ngroute = method_param[9],
                ngroute_neighborhood_charging_depots_size = method_param[11],
                time_limit = TIME_LIMIT_SEC,
            )
        end
    end



### Scratch work

include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("decomposition_heuristic.jl")
include("utils.jl")

using CSV, DataFrames

using Test

data, G = generate_instance(
    ;
    n_depots = 4,
    n_customers = 10,
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
    batch = 1,
    permissiveness = 0.2,
)
plot_instance(data)

for (time_windows, path_single_service, path_check_customers, christofides, ngroute) in [
    (false, false, false, false, false),
    (false,  true,  true, false, false),
    (false, false, false,  true, false),
    (false,  true,  true,  true, false),
    (false, false, false, false,  true),
    (false, false, false,  true,  true),
    ( true, false, false, false, false),
    ( true,  true,  true, false, false),
    ( true, false, false,  true, false),
    ( true,  true,  true,  true, false),
    ( true, false, false, false,  true),
    ( true, false, false,  true,  true),
]
    println("$time_windows, $path_single_service, $path_check_customers, $christofides, $ngroute")
    @suppress path_formulation_column_generation_nocharge(
        data, G,
        ;
        time_windows = time_windows,
        path_single_service = path_single_service,
        path_check_customers = path_check_customers,
        christofides = christofides,
        ngroute = ngroute,
    )
end

CGLP_results, CGIP_results, CG_params, CG_printlist, CG_some_paths = path_formulation_column_generation(
    data;
    time_windows = false,
    path_single_service = false,
    path_check_customers = false,
    christofides = true,
    ngroute = true,
)
CGLP_results["objective"]
DH_results_paths = path_formulation_decomposition_heuristic(
    data, G;
    time_windows = false,
    path_single_service = false,
    path_check_customers = false,
    christofides = true,
    ngroute = true,
    use_integer_paths = true,
)
compute_objective_from_path_solution(DH_results_paths, data, graph)


time_windows = false
path_single_service = false
path_check_customers = false
christofides = false
ngroute = false
ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))) 
verbose = true
time_limit = Inf
max_iters = Inf

time_heuristic_slack = 1.0
infeasible = false
CGLP_results, CGIP_results, CG_params, printlist, some_paths = path_formulation_column_generation_nocharge(
    data, G,
    ;
    time_heuristic_slack = 0.3,
    time_windows = time_windows,
    path_single_service = path_single_service,
    path_check_customers = path_check_customers,
    christofides = christofides,
    ngroute = ngroute,
    ngroute_neighborhood_size = ngroute_neighborhood_size,
    verbose = verbose,
    time_limit = time_limit,
    max_iters = max_iters,
);
results_paths = collect_path_solution_support(CGLP_results, some_paths, data, graph)
results_paths_withcharge = Tuple{Float64, Path}[]
feasible = true
for (val, p) in results_paths
    if p.artificial 
        error("Infeasible solution.")
    end
    nodelist = vcat(p.subpaths[1].arcs[1][1], [a[2] for a in p.subpaths[1].arcs])
    if length(nodelist) == 2 && p.subpaths[1].current_time == 0
        push!(results_paths_withcharge, (val, p))
    else
        (p_feasible, pure_path_label) = get_postcharge_shortest_pure_path_label(
            data, G, nodelist,
            ;
            time_windows = time_windows,
        )
        println(p_feasible)
        feasible = feasible && p_feasible
        if !p_feasible
            break
        end
        p_ = convert_pure_path_label_to_path(pure_path_label, graph)
        push!(results_paths_withcharge, (val, p_))
    end
end
(val, p) = results_paths
nodelist = vcat(p.subpaths[1].arcs[1][1], [a[2] for a in p.subpaths[1].arcs])
infeasible = false
p.subpaths[1].current_time
results_paths_withcharge

time_heuristic_slack = time_heuristic_slack - 0.05




CGLP_results, CGIP_results, CG_params, printlist, some_paths = path_formulation_column_generation_nocharge(
    data, G,
    ;
    time_windows = false,
    path_single_service = false,
    path_check_customers = false,
    christofides = false,
    ngroute = false,
    # max_iters = 1.0,
)
results_paths = collect_path_solution_support(CGLP_results, some_paths, data, graph)
nodelists = [
    vcat(p.subpaths[1].arcs[1][1], [a[2] for a in p.subpaths[1].arcs])
    for (val, p) in results_paths
]
results_paths
time_windows = false



start_time = time()
modified_costs = compute_arc_modified_costs(graph, data, zeros(Float64, graph.n_customers))
t = graph.t
B = graph.B
q = graph.q
if time_windows
    α = graph.α
    β = graph.β
else
    α = zeros(Int, graph.n_nodes)
    β = repeat([graph.T], graph.n_nodes)
end

nodelist = nodelists[4]

pure_path_labels = Dict(
    (nodelist[1:j]...,) => Dict(
        current_node => SortedDict{
            Tuple{Vararg{Int}}, 
            PurePathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for current_node in graph.N_nodes
    )
    for j in eachindex(nodelist)
)
key = (0, -graph.B, graph.B)
starting_node = nodelist[1]
nodeseq = (starting_node,)
pure_path_labels[nodeseq][starting_node][key] = PurePathLabel(
    0.0, 
    [starting_node],
    Int[],
    Int[],
    0,
    0,
    graph.B,
    graph.B,
    false,
    zeros(Int, graph.n_customers),
    false,
)
unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
push!(unexplored_states, (key..., starting_node, nodeseq...))

nodelist
pure_path_labels[(12,)]
pure_path_labels[(12, 1,)]
pure_path_labels[(12, 1, 7,)]
pure_path_labels[(12, 1, 7, 8,)]
pure_path_labels[(12, 1, 7, 8, 2,)]
pure_path_labels[(12, 1, 7, 8, 2, 13)]

first(pure_path_labels[(12, 1, 7, 8, 2, 13)][13])[2]
unexplored_states
while length(unexplored_states) > 0
# begin
    # if time_limit < time() - start_time
    #     throw(TimeLimitException())
    # end
    state = pop!(unexplored_states)
    current_node = state[4]
    current_nodeseq = state[5:end]
    current_key = state[1:3]
    if !(current_key in keys(pure_path_labels[current_nodeseq][current_node]))
        continue
    end
    current_path = pure_path_labels[current_nodeseq][current_node][current_key]
    println("current_path: $(current_path.nodes)")
    eventual_next_node = nodelist[length(current_nodeseq) + 1]
    for next_node in setdiff(vcat(eventual_next_node, graph.N_charging), current_node)
        # feasibility checks
        # (1) battery
        excess = max(
            0, 
            q[current_node,next_node] - current_path.charge_mincharge 
        )
        # (2) time windows
        if current_path.time_mincharge + excess + t[current_node,next_node] > β[next_node]
            # println("$(current_path.time_mincharge), $excess, $(t[current_node,next_node]), $(β[next_node])")
            # println("not time windows feasible")
            continue
        end
        if current_path.time_mincharge + excess + t[current_node,next_node] + graph.min_t[next_node] > graph.T
            continue
        end
        # (3) charge interval 
        if (
            (current_node in graph.N_charging && excess > max(B - current_path.charge_mincharge, 0))
            || 
            (!(current_node in graph.N_charging) && excess > max(current_path.charge_maxcharge - current_path.charge_mincharge, 0))
        )
            # if current_node in graph.N_charging
            #     println("$excess, $(B), $(current_path.charge_mincharge)")
            # else
            #     println("$excess, $(current_path.charge_maxcharge), $(current_path.charge_mincharge)")
            # end
            # println("not charge feasible")
            continue
        end
        
        new_path = copy(current_path)
        push!(new_path.nodes, next_node)
        if next_node in graph.N_customers
            new_path.served[next_node] += 1
        end

        push!(new_path.excesses, excess)
        new_path.time_mincharge = max(
            α[next_node],
            current_path.time_mincharge + t[current_node,next_node] + excess
        )
        if current_node in graph.N_charging
            slack = max(
                # floating point accuracy
                0, 
                min(
                    new_path.time_mincharge - (current_path.time_mincharge + t[current_node,next_node] + excess),
                    B - (current_path.charge_mincharge + excess),
                )
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[next_node],
                max(
                    α[next_node],
                    current_path.time_mincharge + (B - current_path.charge_mincharge) + t[current_node,next_node],
                )
            )
        else
            slack = max(
                # floating point accuracy
                0, 
                min(
                    new_path.time_mincharge - (current_path.time_mincharge + t[current_node,next_node] + excess),
                    current_path.charge_maxcharge - (current_path.charge_mincharge + excess),
                )
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[next_node],
                max(
                    α[next_node],
                    current_path.time_maxcharge + t[current_node,next_node],
                )
            )
        end
        
        new_path.charge_mincharge = (
            current_path.charge_mincharge 
            + excess 
            + slack
            - q[current_node,next_node]
        )
        new_path.charge_maxcharge = (
            new_path.charge_mincharge 
            + new_path.time_maxcharge 
            - new_path.time_mincharge
        )

        new_path.cost += modified_costs[current_node,next_node]
        new_path.cost += data.charge_cost_coeff * (slack + excess)

        # add new_path to collection
        if next_node in union(graph.N_customers, graph.N_depots)
            new_nodeseq = (current_nodeseq..., next_node)
        else
            new_nodeseq = current_nodeseq
        end
        new_key = (
            new_path.time_mincharge, 
            - new_path.charge_maxcharge, 
            new_path.time_mincharge - new_path.charge_mincharge,
        )
        println("$next_node, $new_nodeseq, $new_key")
        added = add_pure_path_label_to_collection!(
            pure_path_labels[new_nodeseq][next_node], 
            new_key, new_path, 
            ;
            verbose = true,
        )
        if added && !(next_node in graph.N_depots)
            new_state = (new_key..., next_node, new_nodeseq...,)
            push!(unexplored_states, new_state)
        end
    end
    break
end

println(nodelist)
pure_path_labels = Dict(
    (nodelist[1:j]...,) => Dict(
        current_node => SortedDict{
            Tuple{Vararg{Int}}, 
            PurePathLabel,
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for current_node in graph.N_nodes
    )
    for j in eachindex(nodelist)
)
key = (0, -graph.B, graph.B, -zeros(Int, graph.n_customers)...)
unexplored_states = SortedSet{Tuple{Vararg{Int}}}()
starting_node = nodelist[1]
pure_path_labels[starting_node][starting_node][key] = PurePathLabel(
    0.0, 
    [starting_node],
    Int[],
    Int[],
    0,
    0,
    graph.B,
    graph.B,
    false,
    zeros(Int, graph.n_customers),
    false,
)
push!(unexplored_states, (key..., starting_node, 1))


while length(unexplored_states) > 0
# begin
    # if time_limit < time() - start_time
    #     throw(TimeLimitException())
    # end
    state = pop!(unexplored_states)
    current_node = state[end-1]
    current_key = state[1:end-2]
    if !(current_key in keys(pure_path_labels[starting_node][current_node]))
        continue
    end
    current_path = pure_path_labels[starting_node][current_node][current_key]
    eventual_next_node = nodelist[state[end] + 1]
    for next_node in vcat(eventual_next_node, graph.N_charging)
        # feasibility checks
        # (1) battery
        excess = max(
            0, 
            q[current_node,next_node] - current_path.charge_mincharge 
        )
        # (2) time windows
        if current_path.time_mincharge + excess + t[current_node,next_node] > β[next_node]
            # println("$(current_path.time_mincharge), $excess, $(t[current_node,next_node]), $(β[next_node])")
            # println("not time windows feasible")
            continue
        end
        if current_path.time_mincharge + excess + t[current_node,next_node] + graph.min_t[next_node] > graph.T
            continue
        end
        # (3) charge interval 
        if (
            (current_node in graph.N_charging && excess > max(B - current_path.charge_mincharge, 0))
            || 
            (!(current_node in graph.N_charging) && excess > max(current_path.charge_maxcharge - current_path.charge_mincharge, 0))
        )
            # if current_node in graph.N_charging
            #     println("$excess, $(B), $(current_path.charge_mincharge)")
            # else
            #     println("$excess, $(current_path.charge_maxcharge), $(current_path.charge_mincharge)")
            # end
            # println("not charge feasible")
            continue
        end
        
        new_path = copy(current_path)
        push!(new_path.nodes, next_node)
        if next_node in graph.N_customers
            new_path.served[next_node] += 1
        end

        push!(new_path.excesses, excess)
        new_path.time_mincharge = max(
            α[next_node],
            current_path.time_mincharge + t[current_node,next_node] + excess
        )
        if current_node in graph.N_charging
            slack = max(
                # floating point accuracy
                0, 
                min(
                    new_path.time_mincharge - (current_path.time_mincharge + t[current_node,next_node] + excess),
                    B - (current_path.charge_mincharge + excess),
                )
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[next_node],
                max(
                    α[next_node],
                    current_path.time_mincharge + (B - current_path.charge_mincharge) + t[current_node,next_node],
                )
            )
        else
            slack = max(
                # floating point accuracy
                0, 
                min(
                    new_path.time_mincharge - (current_path.time_mincharge + t[current_node,next_node] + excess),
                    current_path.charge_maxcharge - (current_path.charge_mincharge + excess),
                )
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[next_node],
                max(
                    α[next_node],
                    current_path.time_maxcharge + t[current_node,next_node],
                )
            )
        end
        
        new_path.charge_mincharge = (
            current_path.charge_mincharge 
            + excess 
            + slack
            - q[current_node,next_node]
        )
        new_path.charge_maxcharge = (
            new_path.charge_mincharge 
            + new_path.time_maxcharge 
            - new_path.time_mincharge
        )

        new_path.cost += modified_costs[current_node,next_node]
        new_path.cost += data.charge_cost_coeff * (slack + excess)

        # add new_path to collection
        new_key = (
            new_path.time_mincharge, 
            - new_path.charge_maxcharge, 
            new_path.time_mincharge - new_path.charge_mincharge, 
            - new_path.served...,
        )
        added = add_pure_path_label_to_collection!(
            pure_path_labels[starting_node][next_node], 
            new_key, new_path, 
            ;
            verbose = false,
        )
        if added && !(next_node in graph.N_depots)
            if next_node in graph.N_customers
                new_state = (new_key..., next_node, state[end] + 1)
            else
                new_state = (new_key..., next_node, state[end])
            end
            push!(unexplored_states, new_state)
        end
    end
end




include("utils.jl")
include("path_formulation.jl")
include("subpath_formulation.jl")
include("decomposition_heuristic.jl")

using Profile


data_params = collect(Iterators.product(
    [4], # n_depots
    20:4:40, # n_customers
    [7], # n_charging
    [8], # n_vehicles
    40000:10000:80000, # T
    [15000], # B
))
method_params = [
    # method
    # path_single_service
    # path_check_customers
    # christofides
    # ngroute_check_create_fset
    ("ours", false, false, false, false),
    ("ours", false, false,  true, false),
    ("ours", false, false, false,  true),
    ("ours", false, false,  true,  true),
    ("ours", false, false, false,  true),
    ("ours", false, false,  true,  true),
    ("ours",  true,  true, false, false),
]

data, G = generate_instance(
    ;
    n_depots = 4,
    n_customers = 20,
    n_charging = 7,
    n_vehicles = 8,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
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
)
for method_param in method_params
    use_time_windows = false
    method, path_single_service, path_check_customers, christofides, ngroute = method_param
    a = @allocated @suppress begin 
        CGLP_results, CGIP_results, CG_params, CG_printlist, CG_some_paths = path_formulation_column_generation(
            data,
            ;
            method = method,
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
        )
        DHL_results_paths = path_formulation_decomposition_heuristic(
            data, G,
            ;
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            use_integer_paths = false,
        )
        DHI_results_paths = path_formulation_decomposition_heuristic(
            data, G,
            ;
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            use_integer_paths = true,
        )
        DHL_objective = compute_objective_from_path_solution(DHL_results_paths, data, graph)
        DHI_objective = compute_objective_from_path_solution(DHL_results_paths, data, graph)
    end
    println("$method_param: $a")
end

data
(method, path_single_service, path_check_customers, christofides, ngroute) = ("ours", false, false, false, false)
use_time_windows = false

@profile path_formulation_column_generation(
    data,
    ;
    method = method,
    time_windows = use_time_windows,
    path_single_service = path_single_service,
    path_check_customers = path_check_customers,
    christofides = christofides,
    ngroute = ngroute,
)

Profile.print(format = :flat)

@code_warntype  path_formulation_column_generation(
    data,
    ;
    method = method,
    time_windows = use_time_windows,
    path_single_service = path_single_service,
    path_check_customers = path_check_customers,
    christofides = christofides,
    ngroute = ngroute,
)


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
neighborhoods = compute_ngroute_neighborhoods(graph, Int(floor(sqrt(graph.n_customers))))

some_paths = generate_artificial_paths(data, graph)
path_costs = compute_path_costs(
    data, graph,
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)
CGLP_results = Dict()
CG_params = Dict()
CG_params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
CG_params["objective"] = Float64[]
CG_params["κ"] = Dict{Int, Float64}[]
CG_params["μ"] = Dict{Int, Float64}[]
CG_params["ν"] = Vector{Float64}[]
CG_params["lp_relaxation_solution_time_taken"] = Float64[]
CG_params["sp_base_time_taken"] = Float64[]
CG_params["sp_full_time_taken"] = Float64[]
CG_params["sp_total_time_taken"] = Float64[]
CG_params["lp_relaxation_constraint_time_taken"] = Float64[]
CG_params["number_of_new_paths"] = Int[]

printlist = String[]
counter = 0
converged = false

if isnothing(Env)
    model = @suppress Model(Gurobi.Optimizer)
else
    model = @suppress Model(() -> Gurobi.Optimizer(Env))
end
JuMP.set_attribute(model, "MIPGapAbs", 1e-3)
JuMP.set_string_names_on_creation(model, false)
z = Dict{
    Tuple{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Int,
    }, 
    VariableRef
}(
    (key, p) => @variable(model, lower_bound = 0)
    for key in keys(some_paths)
        for p in 1:length(some_paths[key])
)
@constraint(
    model,
    κ[i in graph.N_depots],
    sum(
        sum(
            z[((i,0,graph.B),state2),p]
            for p in 1:length(some_paths[((i,0,graph.B),state2)])
        )        
        for (state1, state2) in keys(some_paths)
            if state1[1] == i && state1[2] == 0 && state1[3] == graph.B
    )
    == data.v_start[findfirst(x -> (x == i), graph.N_depots)]
)
@constraint(
    model,
    μ[n2 in graph.N_depots],
    sum(
        sum(
            z[(state1, state2),p]
            for p in 1:length(some_paths[(state1, state2)])
        )
        for (state1, state2) in keys(some_paths)
            if state2[1] == n2
    ) ≥ data.v_end[n2]
)
@constraint(
    model,
    ν[j in graph.N_customers],
    sum(
        sum(
            path_service[((state1, state2),j)][p] * z[(state1, state2),p]
            for p in 1:length(some_paths[(state1, state2)])
        )
        for (state1, state2) in keys(some_paths)
    ) == 1
)
@expression(
    model,
    path_costs_expr,
    sum(
        sum(
            path_costs[state_pair][p] * z[state_pair,p]
            for p in 1:length(some_paths[state_pair])
        )
        for state_pair in keys(some_paths)
    )
)
@objective(model, Min, path_costs_expr)

# while (
#     !converged
#     && time_limit ≥ (time() - start_time)
# )
begin
    counter += 1
    mp_solution_start_time = time()
    @suppress optimize!(model)
    mp_solution_end_time = time()
    CGLP_results = Dict(
        "model" => model,
        "objective" => objective_value(model),
        "z" => Dict(
            (key, p) => value.(z[(key, p)])
            for (key, p) in keys(z)
        ),
        "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
        "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
        "ν" => dual.(model[:ν]).data,
        "solution_time_taken" => round(mp_solution_end_time - mp_solution_start_time, digits = 3),
    )
    push!(CG_params["objective"], CGLP_results["objective"])
    push!(CG_params["κ"], CGLP_results["κ"])
    push!(CG_params["μ"], CGLP_results["μ"])
    push!(CG_params["ν"], CGLP_results["ν"])
    push!(CG_params["lp_relaxation_solution_time_taken"], CGLP_results["solution_time_taken"])

    if method == "ours"
        (negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
            data, G,
            CGLP_results["κ"], 
            CGLP_results["μ"], 
            CGLP_results["ν"], 
            Dict{NTuple{3, Int}, Float64}(),
            ;
            elementary = elementary,
        )
        (generated_paths) = get_paths_from_negative_path_labels(
            graph, negative_full_labels,
        )
        push!(
            CG_params["sp_base_time_taken"],
            round(base_labels_time, digits=3)
        )
        push!(
            CG_params["sp_full_time_taken"],
            round(full_labels_time, digits=3)
        )
        push!(
            CG_params["sp_total_time_taken"],
            round(base_labels_time + full_labels_time, digits=3)
        )
    elseif method == "benchmark"
        (negative_pure_path_labels, _, pure_path_labels_time) = subproblem_iteration_benchmark(
            data, graph, 
            CGLP_results["κ"], 
            CGLP_results["μ"], 
            CGLP_results["ν"], 
            Dict{NTuple{3, Int}, Float64}(),
            ;
            time_windows = time_windows,
        )
        generated_paths = get_paths_from_negative_pure_path_labels(
            graph, negative_pure_path_labels,
        )
        push!(
            CG_params["sp_base_time_taken"],
            0.0
        )
        push!(
            CG_params["sp_full_time_taken"],
            round(pure_path_labels_time, digits=3)
        )
        push!(
            CG_params["sp_total_time_taken"],
            round(pure_path_labels_time, digits=3)
        )
    end

    if length(generated_paths) == 0
        push!(CG_params["number_of_new_paths"], 0)
        converged = true
    else
        push!(
            CG_params["number_of_new_paths"],
            sum(length(v) for v in values(generated_paths))
        )
    end

    mp_constraint_start_time = time()
    for state_pair in keys(generated_paths)
        if !(state_pair in keys(some_paths))
            some_paths[state_pair] = []
            path_costs[state_pair] = []
            for i in 1:graph.n_customers
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
                    compute_path_cost(data, graph, p_new)
                )
                # 3: add path service
                for i in 1:graph.n_customers
                    push!(path_service[(state_pair, i)], p_new.served[i])
                end
                # 4: create variable
                count += 1
                z[(state_pair, count)] = @variable(model, lower_bound = 0)
                (state1, state2) = state_pair
                # 5: modify constraints starting from depot, ending at depot
                set_normalized_coefficient(κ[state1[1]], z[state_pair,count], 1)
                set_normalized_coefficient(μ[state2[1]], z[state_pair,count], 1)
                # 6: modify customer service constraints
                for l in graph.N_customers
                    set_normalized_coefficient(ν[l], z[state_pair, count], p_new.served[l])
                end
                # 7: modify objective
                set_objective_coefficient(model, z[state_pair, count], path_costs[state_pair][count])
            end
        end
    end
    mp_constraint_end_time = time()

    push!(
        CG_params["number_of_paths"], 
        sum(length(v) for v in values(some_paths))
    )
    push!(
        CG_params["lp_relaxation_constraint_time_taken"],
        round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
    )
end
CG_params["number_of_paths"]
CG_params["sp_full_time_taken"]

### Starting debug of generate_base_labels_ngroute
include("utils.jl")
include("subpath_stitching.jl")
neighborhoods = compute_ngroute_neighborhoods(graph, 4)
# neighborhoods = compute_ngroute_neighborhoods(graph, Int(floor(sqrt(graph.n_customers))))

base_labels = @time generate_base_labels_ngroute(data, graph, neighborhoods, κ, μ, ν, christofides = true)
full_labels = @time find_nondominated_paths_notimewindows_ngroute(data, graph, neighorhoods, base_labels, κ, μ, christofides = true)


base_labels = @time generate_base_labels_singleservice(graph, κ, μ, ν, check_customers = true, christofides = true)
full_labels = @time find_nondominated_paths_notimewindows(data, graph, base_labels, κ, μ, single_service = true, check_customers = true, christofides = true)

base_labels = @time generate_base_labels_nonsingleservice(graph, κ, μ, ν, christofides = true)
full_labels = @time find_nondominated_paths_notimewindows(data, graph, base_labels, κ, μ, single_service = false, check_customers = false, christofides = true)

[p.served for set in keys(full_labels[20][20]) for p in values(full_labels[20][20][set])]

κ = CGLP_results["κ"]
μ = CGLP_results["μ"]
ν = CGLP_results["ν"]

base_labels[19][22]
base_labels[22][21]
sort(neighborhoods)

neighborhoods
collect((1,2,))

function ngroute_extend_partial_path_check(
    neighborhoods::BitMatrix,
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
                if neighborhoods[node, next_node]
        ]
        push!(new_set, next_node)
        println("$next_node, $new_set")
    end
    return (Tuple(sort(unique(new_set))), true)
end

set = (8,16,22)
s = base_labels[22][21][(7,21)][end]
ngroute_extend_partial_path_check(neighborhoods, set, s)

function add_subpath_longlabel_to_collection!(
    collection::SortedDict{
        Int,
        BaseSubpathLabel,
        Base.Order.ForwardOrdering,
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

modified_costs = compute_arc_modified_costs(graph, data, ν)

base_labels = Dict(
    starting_node => Dict(
        current_node => Dict{
            Tuple{Vararg{Int}}, 
            SortedDict{Int, BaseSubpathLabel},
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for current_node in graph.N_nodes
    )
    for starting_node in graph.N_depots_charging
)
for node in graph.N_depots_charging
    base_labels[node][node][(node,)] = SortedDict{
        Int, 
        BaseSubpathLabel,
    }(
        Base.Order.ForwardOrdering(),
        0 => BaseSubpathLabel(
            0, 0, 0.0, [node,], zeros(Int, graph.n_customers),
        )
    )
end

unexplored_states = SortedSet(
    [
        (0.0, node, node)
        for node in graph.N_depots_charging
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
            if next_node in graph.N_customers
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
            # if current_subpath.time_taken + graph.t[current_node, next_node] + graph.min_t[next_node] > graph.T
            #     continue
            # end 
            if current_subpath.charge_taken + graph.q[current_node, next_node] + graph.min_q[next_node] > graph.B
                println("not charge feasible")
                continue
            end
            new_subpath = copy(current_subpath)
            new_subpath.time_taken += graph.t[current_node, next_node]
            new_subpath.charge_taken += graph.q[current_node, next_node]
            new_subpath.cost += modified_costs[current_node, next_node]
            push!(new_subpath.nodes, next_node)
            if next_node in graph.N_customers
                new_subpath.served[next_node] += 1
            end
            new_set = ngroute_create_set(neighborhoods, set, next_node)
            println("set: $set, next_node: $next_node, new_set: $new_set")
            if !(new_set in keys(base_labels[starting_node][next_node]))
                base_labels[starting_node][next_node][new_set] = SortedDict{
                    Tuple{Vararg{Int}},
                    BaseSubpathLabel,
                    Base.Order.ForwardOrdering,
                }(Base.Order.ForwardOrdering())
            end
            added = add_subpath_longlabel_to_collection!(
                base_labels[starting_node][next_node][new_set],
                new_subpath.time_taken, new_subpath,
                ;
                verbose = true,
            )
            println("added = $added")
            if added && next_node in graph.N_customers
                new_state = (new_subpath.time_taken, starting_node, next_node)
                push!(unexplored_states, new_state)
                println("added next state: $new_state")
            end
        end
    end
end

unexplored_states
base_labels

for starting_node in graph.N_depots_charging
    for end_node in graph.N_customers
        delete!(base_labels[starting_node], end_node)
    end
end

for starting_node in graph.N_depots
    for end_node in graph.N_depots_charging
        for set in keys(base_labels[starting_node][end_node])
            for v in values(base_labels[starting_node][end_node][set])
                v.cost = v.cost - κ[starting_node]
            end
        end
    end
end
for end_node in graph.N_depots
    for starting_node in graph.N_depots_charging
        for set in keys(base_labels[starting_node][end_node])
            for v in values(base_labels[starting_node][end_node][set])
                v.cost = v.cost - μ[end_node]
            end
        end
    end
end

# remove self-loops with nonnegative cost
for node in graph.N_depots_charging
    for set in keys(base_labels[node][node])
        for (k, v) in pairs(base_labels[node][node][set])
            if v.cost ≥ 0.0
                pop!(base_labels[node][node][set], k)
            end
        end
    end
end

*("3", "3", "4")

for starting_node in graph.N_depots_charging
    vals = [
        length(keys(base_labels[starting_node][end_node]))
        for end_node in graph.N_nodes
    ]
    # vals = [
    #     sum(
    #         [length(base_labels[starting_node][end_node][set])
    #         for set in keys(base_labels[starting_node][end_node])],
    #         init = 0,
    #     )
    #     for end_node in graph.N_nodes
    # ]
    println(*(["$val\t" for val in vals]...))
end


base_labels[28][24]
sum(
    length(base_labels[starting_node][end_node][set])
    for starting_node in graph.N_depots_charging
        for end_node in graph.N_depots_charging
            for set in keys(base_labels[starting_node][end_node])
)
sum(
    length(base_labels[starting_node][end_node])
    for starting_node in graph.N_depots_charging
        for end_node in graph.N_depots_charging
)

plot_instance(data)

### Starting debug of subproblem_iteration_benchmark_incremental_elementarity

κ, μ, ν = CGLP_results["κ"], CGLP_results["μ"], CGLP_results["ν"]
time_windows = false
path_check_customers = true
verbose = true 
warm_start = false 

start_time = time()
S = Int[]
initial_pure_path_labels = nothing

pure_path_labels = find_nondominated_paths_S(
    data, κ, μ, ν,
    ;
    S = S,
    time_windows = time_windows, 
    check_customers = path_check_customers,
    christofides = christofides,
    initial_pure_path_labels = initial_pure_path_labels
)
pure_path_labels_count = sum(
    length(pure_path_labels[starting_node][end_node]) 
    for starting_node in graph.N_depots
        for end_node in keys(pure_path_labels[starting_node]) 
)
# filter for path labels: (i) ending at depot, 
# (ii) w/ negative reduced cost,
# (iii) elementary
d_pure_path_labels = get_depot_pure_path_labels(data, pure_path_labels);
n_d_pure_path_labels = get_negative_pure_path_labels_from_pure_path_labels(graph, d_pure_path_labels);
n_d_pure_path_labels_count = sum(
    length(n_d_pure_path_labels[starting_node][end_node]) 
    for starting_node in graph.N_depots
        for end_node in keys(n_d_pure_path_labels[starting_node]) 
)
if n_d_pure_path_labels_count == 0
    # this activates if the overall CG converges
    return (n_d_pure_path_labels, 0, time() - start_time)
end
(e_n_d_pure_path_labels, ne_n_d_pure_path_labels) = get_elementary_nonelementary_pure_path_labels(
    data, n_d_pure_path_labels; 
    S = graph.N_customers
);
e_n_d_pure_path_labels_count = sum(
    length(e_n_d_pure_path_labels[starting_node][end_node]) 
    for starting_node in graph.N_depots
        for end_node in keys(e_n_d_pure_path_labels[starting_node]) 
)
if e_n_d_pure_path_labels_count > 0
    # this activates if there is some negative path found (elementary w.r.t. all customers)
    return (e_n_d_pure_path_labels, e_n_d_pure_path_labels_count, time() - start_time)
end

# else, expand S
# println(length(S))
if length(S) == graph.n_customers
    error()
end

max_freq = 0
max_custs = Int[]
for starting_node in graph.N_depots
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