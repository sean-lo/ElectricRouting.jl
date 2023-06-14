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

p_b_tw_ie_LP_results, p_b_tw_ie_IP_results, p_b_tw_ie_params, p_b_tw_ie_printlist, p_b_tw_ie_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true,)
p_b_tw_iec_LP_results, p_b_tw_iec_IP_results, p_b_tw_iec_params, p_b_tw_iec_printlist, p_b_tw_iec_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true)
p_b_tw_iecw_LP_results, p_b_tw_iecw_IP_results, p_b_tw_iecw_params, p_b_tw_iecw_printlist, p_b_tw_iecw_some_paths = path_formulation_column_generation(G, data; method = "benchmark", time_windows = true, verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true) # FIXME: debug

p_b_LP_results, p_b_IP_results, p_b_params, p_b_printlist, p_b_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true)
p_b_s_LP_results, p_b_s_IP_results, p_b_s_params, p_b_s_printlist, p_b_s_some_paths = path_formulation_column_generation(G, data; method = "benchmark", path_single_service = true, verbose = true)
p_b_sc_LP_results, p_b_sc_IP_results, p_b_sc_params, p_b_sc_printlist, p_b_sc_some_paths = path_formulation_column_generation(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)

p_b_ie_LP_results, p_b_ie_IP_results, p_b_ie_params, p_b_ie_printlist, p_b_ie_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true,)
p_b_iec_LP_results, p_b_iec_IP_results, p_b_iec_params, p_b_iec_printlist, p_b_iec_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true)
p_b_iecw_LP_results, p_b_iecw_IP_results, p_b_iecw_params, p_b_iecw_printlist, p_b_iecw_some_paths = path_formulation_column_generation(G, data; method = "benchmark", verbose = true, incremental_elementarity = true, path_check_customers = true, warm_start = true) # FIXME: debug this

p_o_LP_results, p_o_IP_results, p_o_params, p_o_printlist, p_o_some_paths = path_formulation_column_generation(G, data; method = "ours", verbose = true)
p_o_s_LP_results, p_o_s_IP_results, p_o_s_params, p_o_s_printlist, p_o_s_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, verbose = true)
p_o_sc_LP_results, p_o_sc_IP_results, p_o_sc_params, p_o_sc_printlist, p_o_sc_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
p_o_ss_LP_results, p_o_ss_IP_results, p_o_ss_params, p_o_ss_printlist, p_o_ss_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
p_o_scsc_LP_results, p_o_scsc_IP_results, p_o_scsc_params, p_o_scsc_printlist, p_o_scsc_some_paths = path_formulation_column_generation(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)

sp_b_tw_LP_results, sp_b_tw_IP_results, sp_b_tw_params, sp_b_tw_printlist, sp_b_tw_some_subpaths, sp_b_tw_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, verbose = true)
sp_b_tw_s_LP_results, sp_b_tw_s_IP_results, sp_b_tw_s_params, sp_b_tw_s_printlist, sp_b_tw_s_some_subpaths, sp_b_tw_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, path_single_service = true, verbose = true)
sp_b_tw_sc_LP_results, sp_b_tw_sc_IP_results, sp_b_tw_sc_params, sp_b_tw_sc_printlist, sp_b_tw_sc_some_subpaths, sp_b_tw_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, path_single_service = true, path_check_customers = true, verbose = true)
sp_b_tw_sca_LP_results, sp_b_tw_sca_IP_results, sp_b_tw_sca_params, sp_b_tw_sca_printlist, sp_b_tw_sca_some_subpaths, sp_b_tw_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

sp_b_tw_ie_LP_results, sp_b_tw_ie_IP_results, sp_b_tw_ie_params, sp_b_tw_ie_printlist, sp_b_tw_ie_some_subpaths, sp_b_tw_ie_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, verbose = true)
sp_b_tw_iec_LP_results, sp_b_tw_iec_IP_results, sp_b_tw_iec_params, sp_b_tw_iec_printlist, sp_b_tw_iec_some_subpaths, sp_b_tw_iec_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", time_windows = true, incremental_elementarity = true, path_check_customers = true, verbose = true)

sp_b_LP_results, sp_b_IP_results, sp_b_params, sp_b_printlist, sp_b_some_subpaths, sp_b_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", verbose = true)
sp_b_s_LP_results, sp_b_s_IP_results, sp_b_s_params, sp_b_s_printlist, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, verbose = true)
sp_b_sc_LP_results, sp_b_sc_IP_results, sp_b_sc_params, sp_b_sc_printlist, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, verbose = true)
sp_b_sca_LP_results, sp_b_sca_IP_results, sp_b_sca_params, sp_b_sca_printlist, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "benchmark", path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

sp_o_LP_results, sp_o_IP_results, sp_o_params, sp_o_printlist, sp_o_some_subpaths, sp_o_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", verbose = true)
sp_o_s_LP_results, sp_o_s_IP_results, sp_o_s_params, sp_o_s_printlist, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, verbose = true)
sp_o_sc_LP_results, sp_o_sc_IP_results, sp_o_sc_params, sp_o_sc_printlist, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, verbose = true)
sp_o_sca_LP_results, sp_o_sca_IP_results, sp_o_sca_params, sp_o_sca_printlist, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, check_customers_accelerated = true, verbose = true)
sp_o_ss_LP_results, sp_o_ss_IP_results, sp_o_ss_params, sp_o_ss_printlist, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, path_single_service = true, verbose = true)
sp_o_scsc_LP_results, sp_o_scsc_IP_results, sp_o_scsc_params, sp_o_scsc_printlist, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, verbose = true)
sp_o_scsca_LP_results, sp_o_scsca_IP_results, sp_o_scsca_params, sp_o_scsca_printlist, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs = subpath_formulation_column_generation_integrated_from_paths(G, data; method = "ours", subpath_single_service = true, subpath_check_customers = true, path_single_service = true, path_check_customers = true, check_customers_accelerated = true, verbose = true)

collect_path_solution_metrics!(p_b_tw_LP_results, data, p_b_tw_some_paths)
collect_path_solution_metrics!(p_b_tw_s_LP_results, data, p_b_tw_s_some_paths)
collect_path_solution_metrics!(p_b_tw_sc_LP_results, data, p_b_tw_sc_some_paths)

collect_path_solution_metrics!(p_b_tw_ie_LP_results, data, p_b_tw_ie_some_paths)
collect_path_solution_metrics!(p_b_tw_iec_LP_results, data, p_b_tw_iec_some_paths)

collect_path_solution_metrics!(p_b_LP_results, data, p_b_some_paths)
collect_path_solution_metrics!(p_b_s_LP_results, data, p_b_s_some_paths)
collect_path_solution_metrics!(p_b_sc_LP_results, data, p_b_sc_some_paths)

collect_path_solution_metrics!(p_b_ie_LP_results, data, p_b_ie_some_paths)
collect_path_solution_metrics!(p_b_iec_LP_results, data, p_b_iec_some_paths)

collect_path_solution_metrics!(p_o_LP_results, data, p_o_some_paths)
collect_path_solution_metrics!(p_o_s_LP_results, data, p_o_s_some_paths)
collect_path_solution_metrics!(p_o_sc_LP_results, data, p_o_sc_some_paths)
collect_path_solution_metrics!(p_o_ss_LP_results, data, p_o_ss_some_paths)
collect_path_solution_metrics!(p_o_scsc_LP_results, data, p_o_scsc_some_paths)

collect_subpath_solution_metrics!(sp_b_tw_LP_results, data, sp_b_tw_some_subpaths, sp_b_tw_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_s_LP_results, data, sp_b_tw_s_some_subpaths, sp_b_tw_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_sc_LP_results, data, sp_b_tw_sc_some_subpaths, sp_b_tw_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_tw_sca_LP_results, data, sp_b_tw_sca_some_subpaths, sp_b_sca_some_charging_arcs)

collect_subpath_solution_metrics!(sp_b_LP_results, data, sp_b_some_subpaths, sp_b_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_s_LP_results, data, sp_b_s_some_subpaths, sp_b_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sc_LP_results, data, sp_b_sc_some_subpaths, sp_b_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_b_sca_LP_results, data, sp_b_sca_some_subpaths, sp_b_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_LP_results, data, sp_o_some_subpaths, sp_o_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_s_LP_results, data, sp_o_s_some_subpaths, sp_o_s_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sc_LP_results, data, sp_o_sc_some_subpaths, sp_o_sc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_sca_LP_results, data, sp_o_sca_some_subpaths, sp_o_sca_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_ss_LP_results, data, sp_o_ss_some_subpaths, sp_o_ss_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsc_LP_results, data, sp_o_scsc_some_subpaths, sp_o_scsc_some_charging_arcs)
collect_subpath_solution_metrics!(sp_o_scsca_LP_results, data, sp_o_scsca_some_subpaths, sp_o_scsca_some_charging_arcs)

p_b_tw_LP_results["objective"]
p_b_tw_s_LP_results["objective"]
p_b_tw_sc_LP_results["objective"]
p_b_tw_ie_LP_results["objective"]
p_b_tw_iec_LP_results["objective"]
p_b_tw_iecw_LP_results["objective"]

p_b_LP_results["objective"]
p_b_s_LP_results["objective"]
p_b_sc_LP_results["objective"]
p_b_ie_LP_results["objective"]
p_b_iec_LP_results["objective"]
p_b_iecw_LP_results["objective"]

p_o_LP_results["objective"]
p_o_s_LP_results["objective"]
p_o_sc_LP_results["objective"]
p_o_ss_LP_results["objective"]
p_o_scsc_LP_results["objective"]

sp_b_tw_LP_results["objective"]
sp_b_tw_s_LP_results["objective"]
sp_b_tw_sc_LP_results["objective"]
sp_b_tw_sca_LP_results["objective"]

sp_b_LP_results["objective"]
sp_b_s_LP_results["objective"]
sp_b_sc_LP_results["objective"]
sp_b_sca_LP_results["objective"]
sp_o_LP_results["objective"]
sp_o_s_LP_results["objective"]
sp_o_sc_LP_results["objective"]
sp_o_sca_LP_results["objective"]
sp_o_ss_LP_results["objective"]
sp_o_scsc_LP_results["objective"]
sp_o_scsca_LP_results["objective"]

@test (
    p_b_LP_results["objective"] 
    ≤ p_b_sc_LP_results["objective"]
)
@test(
    p_o_LP_results["objective"] 
    ≤ p_o_sc_LP_results["objective"]
    ≤ p_o_scsc_LP_results["objective"]
    ≈ p_b_sc_LP_results["objective"]
)
@test p_b_s_LP_results["objective"] ≥ p_b_sc_LP_results["objective"]
@test p_o_s_LP_results["objective"] ≥ p_o_sc_LP_results["objective"]
@test p_o_ss_LP_results["objective"] ≥ p_o_scsc_LP_results["objective"]
@test p_b_LP_results["objective"] ≈ p_o_LP_results["objective"]

@test (
    sp_b_LP_results["objective"] 
    ≤ sp_b_sc_LP_results["objective"]
)
@test(
    sp_o_LP_results["objective"] 
    ≤ sp_o_sc_LP_results["objective"]
    ≤ sp_o_scsc_LP_results["objective"]
    ≈ sp_b_sc_LP_results["objective"]
)
@test (
    sp_b_s_LP_results["objective"] 
    ≥ sp_b_sc_LP_results["objective"] 
    ≈ sp_b_sca_LP_results["objective"]
)
@test (
    sp_o_s_LP_results["objective"] 
    ≥ sp_o_sc_LP_results["objective"] - 1e-8
    ≈ sp_o_sca_LP_results["objective"]
)
@test (
    sp_o_ss_LP_results["objective"] 
    ≥ sp_o_scsc_LP_results["objective"]
    ≈ sp_o_scsca_LP_results["objective"]
)
@test sp_b_LP_results["objective"] ≈ sp_o_LP_results["objective"]
@test p_b_LP_results["objective"] ≈ sp_b_LP_results["objective"]
@test p_b_sc_LP_results["objective"] ≈ sp_b_sc_LP_results["objective"]
@test p_o_LP_results["objective"] ≈ sp_o_LP_results["objective"]
@test p_o_s_LP_results["objective"] ≈ sp_o_s_LP_results["objective"]
@test p_o_sc_LP_results["objective"] ≈ sp_o_sc_LP_results["objective"]
@test p_o_scsc_LP_results["objective"] ≈ sp_o_scsc_LP_results["objective"]


p_b_tw_params["time_taken"]
p_b_tw_s_params["time_taken"]
p_b_tw_sc_params["time_taken"]
p_b_tw_ie_params["time_taken"]
p_b_tw_iec_params["time_taken"]
p_b_tw_iecw_params["time_taken"]

p_b_params["time_taken"]
p_b_s_params["time_taken"]
p_b_sc_params["time_taken"]
p_b_ie_params["time_taken"]
p_b_iec_params["time_taken"]
p_b_iecw_params["time_taken"]

p_o_params["time_taken"]
p_o_s_params["time_taken"]
p_o_sc_params["time_taken"]
p_o_ss_params["time_taken"]
p_o_scsc_params["time_taken"]

sp_b_tw_params["time_taken"]
sp_b_tw_s_params["time_taken"]
sp_b_tw_sc_params["time_taken"]
sp_b_tw_sca_params["time_taken"]

sp_b_params["time_taken"]
sp_b_s_params["time_taken"]
sp_b_sc_params["time_taken"]
sp_b_sca_params["time_taken"]
sp_o_params["time_taken"]
sp_o_s_params["time_taken"]
sp_o_sc_params["time_taken"]
sp_o_sca_params["time_taken"]
sp_o_ss_params["time_taken"]
sp_o_scsc_params["time_taken"]
sp_o_scsca_params["time_taken"]


### Scratch work

some_paths = generate_artificial_paths(data)
path_costs = compute_path_costs(
    data, 
    some_paths,
)
path_service = compute_path_service(data, some_paths)
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
    seed = 1,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 4,
    permissiveness = 0.7,
)
G = construct_graph(data)
plot_instance(data)
(
    CGLP_results, CGIP_results, CG_params, CG_printlist, CG_subpaths, CG_charging_arcs 
) = subpath_formulation_column_generation_integrated_from_paths(G, data);
(
    CGLPc_results, CGIPc_results, CGc_params, CGc_printlist, CGc_subpaths, CGc_charging_arcs 
) = subpath_formulation_column_generation_integrated_from_paths(G, data, check_customers = true);
(
    CGLPb_results, CGIPb_results, CGb_params, CGb_printlist, CGb_subpaths, CGb_charging_arcs 
) = subpath_formulation_column_generation_integrated_from_paths(G, data, method = "benchmark");
(
    CGLPbc_results, CGIPbc_results, CGbc_params, CGbc_printlist, CGbc_subpaths, CGbc_charging_arcs 
) = subpath_formulation_column_generation_integrated_from_paths(G, data, method = "benchmark", check_customers = true);

collect_subpath_solution_metrics!(CGLP_results, data, CG_subpaths, CG_charging_arcs)
collect_subpath_solution_metrics!(CGIP_results, data, CG_subpaths, CG_charging_arcs)
collect_subpath_solution_metrics!(CGLPc_results, data, CGc_subpaths, CGc_charging_arcs)
collect_subpath_solution_metrics!(CGIPc_results, data, CGc_subpaths, CGc_charging_arcs)
collect_subpath_solution_metrics!(CGLPb_results, data, CGb_subpaths, CGb_charging_arcs)
collect_subpath_solution_metrics!(CGIPb_results, data, CGb_subpaths, CGb_charging_arcs)
collect_subpath_solution_metrics!(CGLPbc_results, data, CGbc_subpaths, CGbc_charging_arcs)
collect_subpath_solution_metrics!(CGIPbc_results, data, CGbc_subpaths, CGbc_charging_arcs)

CGLP_results["objective"]
CGLPc_results["objective"]
CGLPb_results["objective"]
CGLPbc_results["objective"]

plot_subpath_solution(CGLP_results, data, CG_subpaths, CG_charging_arcs)

all_results = []
for (
    n_depots, 
    n_customers, 
    n_charging,
    n_vehicles,
    T, 
    B, 
) in [
    (4, 20, 6, 8, 40000, 15000),
    (4, 20, 8, 8, 40000, 15000),
    (4, 20, 10, 8, 40000, 15000),
    (4, 20, 12, 8, 40000, 15000),
]
    for seed in 1:10
    # for seed in [1]
        println("$n_depots, $n_customers, $n_charging, $n_vehicles, $T, $B: $seed")
        data = generate_instance(
            ;
            n_depots = n_depots,
            n_customers = n_customers,
            n_charging = n_charging,
            n_vehicles = n_vehicles,
            depot_pattern = "circular",    
            customer_pattern = "random_box",
            charging_pattern = "grid",
            shrinkage_depots = 1.0,
            shrinkage_charging = 0.7,
            T = T,
            seed = seed,
            B = B,
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
        # plot_instance(data)
        # (
        #     b_LP_results, b_IP_results, b_params, b_printlist, b_subpaths, b_charging_arcs
        # ) = subpath_formulation_column_generation_integrated_from_paths(
        #     G, data, method = "benchmark",
        # )
        (
            b_s_LP_results, b_s_IP_results, b_s_params, b_s_printlist, b_s_subpaths, b_s_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "benchmark", 
            path_single_service = true,
        )
        println("benchmark - single service, NO customer labels\t\t\t\t$(b_s_params["time_taken"])\t$(b_s_LP_results["objective"])")
        (
            b_sc_LP_results, b_sc_IP_results, b_sc_params, b_sc_printlist, b_sc_subpaths, b_sc_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "benchmark", 
            path_single_service = true, path_check_customers = true,
        )
        println("benchmark - single service, customer labels\t\t\t\t$(b_sc_params["time_taken"])\t$(b_sc_LP_results["objective"])")
        (
            b_sca_LP_results, b_sca_IP_results, b_sca_params, b_sca_printlist, b_sca_subpaths, b_sca_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "benchmark", 
            path_single_service = true, path_check_customers = true, check_customers_accelerated = true
        )
        println("benchmark - single service, customer labels (accel)\t\t\t$(b_sca_params["time_taken"])\t$(b_sca_LP_results["objective"])")

        # (
        #     o_LP_results, o_IP_results, o_params, o_printlist, o_subpaths, o_charging_arcs
        # ) = subpath_formulation_column_generation_integrated_from_paths(
        #     G, data, method = "ours",
        # )
        (
            o_s_LP_results, o_s_IP_results, o_s_params, o_s_printlist, o_s_subpaths, o_s_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "ours", 
            subpath_single_service = true,
        )
        println("ours - single service, NO customer labels for subpaths\t\t\t$(o_s_params["time_taken"])\t$(o_s_LP_results["objective"])")
        (
            o_sc_LP_results, o_sc_IP_results, o_sc_params, o_sc_printlist, o_sc_subpaths, o_sc_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "ours", 
            subpath_single_service = true, subpath_check_customers = true,
        )
        println("ours - single service, customer labels for subpaths\t\t\t$(o_sc_params["time_taken"])\t$(o_sc_LP_results["objective"])")

        (
            o_sca_LP_results, o_sca_IP_results, o_sca_params, o_sca_printlist, o_sca_subpaths, o_sca_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "ours", 
            subpath_single_service = true, subpath_check_customers = true, 
            check_customers_accelerated = true,
        )
        println("ours - single service, customer labels for subpaths (accel)\t\t$(o_sca_params["time_taken"])\t$(o_sca_LP_results["objective"])")

        (
            o_ss_LP_results, o_ss_IP_results, o_ss_params, o_ss_printlist, o_ss_subpaths, o_ss_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "ours", 
            subpath_single_service = true,
            path_single_service = true,
        )
        println("ours - single service, NO customer labels for subpaths and paths\t$(o_ss_params["time_taken"])\t$(o_ss_LP_results["objective"])")
        (
            o_scsc_LP_results, o_scsc_IP_results, o_scsc_params, o_scsc_printlist, o_scsc_subpaths, o_scsc_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "ours", 
            subpath_single_service = true, subpath_check_customers = true,
            path_single_service = true, path_check_customers = true,
        )
        println("ours - single service, customer labels for subpaths and paths\t\t$(o_scsc_params["time_taken"])\t$(o_scsc_LP_results["objective"])")        
        (
            o_scsca_LP_results, o_scsca_IP_results, o_scsca_params, o_scsca_printlist, o_scsca_subpaths, o_scsca_charging_arcs
        ) = @suppress subpath_formulation_column_generation_integrated_from_paths(
            G, data, method = "ours", 
            subpath_single_service = true, subpath_check_customers = true, 
            path_single_service = true, path_check_customers = true,
            check_customers_accelerated = true,
        )
        println("ours - single service, customer labels for subpaths and paths (accel)\t$(o_scsca_params["time_taken"])\t$(o_scsca_LP_results["objective"])")

        # collect_subpath_solution_metrics!(b_LP_results, data, b_subpaths, b_charging_arcs)
        # collect_subpath_solution_metrics!(b_IP_results, data, b_subpaths, b_charging_arcs)
        collect_subpath_solution_metrics!(b_s_LP_results, data, b_s_subpaths, b_s_charging_arcs)
        collect_subpath_solution_metrics!(b_s_IP_results, data, b_s_subpaths, b_s_charging_arcs)
        collect_subpath_solution_metrics!(b_sc_LP_results, data, b_sc_subpaths, b_sc_charging_arcs)
        collect_subpath_solution_metrics!(b_sc_IP_results, data, b_sc_subpaths, b_sc_charging_arcs)
        collect_subpath_solution_metrics!(b_sca_LP_results, data, b_sca_subpaths, b_sca_charging_arcs)
        collect_subpath_solution_metrics!(b_sca_IP_results, data, b_sca_subpaths, b_sca_charging_arcs)
        # collect_subpath_solution_metrics!(o_LP_results, data, o_subpaths, o_charging_arcs)
        # collect_subpath_solution_metrics!(o_IP_results, data, o_subpaths, o_charging_arcs)
        collect_subpath_solution_metrics!(o_s_LP_results, data, o_s_subpaths, o_s_charging_arcs)
        collect_subpath_solution_metrics!(o_s_IP_results, data, o_s_subpaths, o_s_charging_arcs)
        collect_subpath_solution_metrics!(o_sc_LP_results, data, o_sc_subpaths, o_sc_charging_arcs)
        collect_subpath_solution_metrics!(o_sc_IP_results, data, o_sc_subpaths, o_sc_charging_arcs)
        collect_subpath_solution_metrics!(o_sca_LP_results, data, o_sca_subpaths, o_sca_charging_arcs)
        collect_subpath_solution_metrics!(o_sca_IP_results, data, o_sca_subpaths, o_sca_charging_arcs)
        collect_subpath_solution_metrics!(o_ss_LP_results, data, o_ss_subpaths, o_ss_charging_arcs)
        collect_subpath_solution_metrics!(o_ss_IP_results, data, o_ss_subpaths, o_ss_charging_arcs)
        collect_subpath_solution_metrics!(o_scsc_LP_results, data, o_scsc_subpaths, o_scsc_charging_arcs)
        collect_subpath_solution_metrics!(o_scsc_IP_results, data, o_scsc_subpaths, o_scsc_charging_arcs)
        collect_subpath_solution_metrics!(o_scsca_LP_results, data, o_scsca_subpaths, o_scsca_charging_arcs)
        collect_subpath_solution_metrics!(o_scsca_IP_results, data, o_scsca_subpaths, o_scsca_charging_arcs)


        # println("benchmark - single service, NO customer labels\t\t\t\t$(b_s_params["time_taken"])\t$(b_s_LP_results["objective"])")
        # println("benchmark - single service, customer labels\t\t\t\t$(b_sc_params["time_taken"])\t$(b_sc_LP_results["objective"])")
        # println("benchmark - single service, customer labels (accel)\t\t\t$(b_sca_params["time_taken"])\t$(b_sca_LP_results["objective"])")
        # println("ours - single service, NO customer labels for subpaths\t\t\t$(o_s_params["time_taken"])\t$(o_s_LP_results["objective"])")
        # println("ours - single service, customer labels for subpaths\t\t\t$(o_sc_params["time_taken"])\t$(o_sc_LP_results["objective"])")
        # println("ours - single service, customer labels for subpaths (accel)\t\t$(o_sca_params["time_taken"])\t$(o_sca_LP_results["objective"])")
        # println("ours - single service, NO customer labels for subpaths and paths\t$(o_ss_params["time_taken"])\t$(o_ss_LP_results["objective"])")
        # println("ours - single service, customer labels for subpaths and paths\t\t$(o_scsc_params["time_taken"])\t$(o_scsc_LP_results["objective"])")        
        # println("ours - single service, customer labels for subpaths and paths (accel)\t$(o_scsca_params["time_taken"])\t$(o_scsca_LP_results["objective"])")

        push!(all_results,
        (
            n_depots = n_depots, 
            n_customers = n_customers, 
            n_charging = n_charging,
            n_vehicles = n_vehicles,
            T = T, 
            B = B, 
            seed = seed,
            # b_LP_objective = b_LP_results["objective"],
            b_s_LP_objective = b_s_LP_results["objective"],
            b_sc_LP_objective = b_sc_LP_results["objective"],
            b_sca_LP_objective = b_sca_LP_results["objective"],
            # o_LP_objective = o_LP_results["objective"],
            o_s_LP_objective = o_s_LP_results["objective"],
            o_sc_LP_objective = o_sc_LP_results["objective"],
            o_sca_LP_objective = o_sca_LP_results["objective"],
            o_ss_LP_objective = o_ss_LP_results["objective"],
            o_scsc_LP_objective = o_scsc_LP_results["objective"],
            o_scsca_LP_objective = o_scsca_LP_results["objective"],

            # b_time_taken = b_params["time_taken"],
            b_s_time_taken = b_s_params["time_taken"],
            b_sc_time_taken = b_sc_params["time_taken"],
            b_sca_time_taken = b_sca_params["time_taken"],
            # o_time_taken = o_params["time_taken"],
            o_s_time_taken = o_s_params["time_taken"],
            o_sc_time_taken = o_sc_params["time_taken"],
            o_sca_time_taken = o_sca_params["time_taken"],
            o_ss_time_taken = o_ss_params["time_taken"],
            o_scsc_time_taken = o_scsc_params["time_taken"],
            o_scsca_time_taken = o_scsca_params["time_taken"],

            # o_SP_base_time_taken = o_params["sp_base_time_taken_total"],
            o_s_SP_base_time_taken = o_s_params["sp_base_time_taken_total"],
            o_sc_SP_base_time_taken = o_sc_params["sp_base_time_taken_total"],
            o_sca_SP_base_time_taken = o_sca_params["sp_base_time_taken_total"],
            o_ss_SP_base_time_taken = o_ss_params["sp_base_time_taken_total"],
            o_scsc_SP_base_time_taken = o_scsc_params["sp_base_time_taken_total"],
            o_scsca_SP_base_time_taken = o_scsca_params["sp_base_time_taken_total"],

            # b_SP_full_time_taken = b_params["sp_full_time_taken_total"],
            b_s_SP_full_time_taken = b_s_params["sp_full_time_taken_total"],
            b_sc_SP_full_time_taken = b_sc_params["sp_full_time_taken_total"],
            b_sca_SP_full_time_taken = b_sca_params["sp_full_time_taken_total"],
            # o_SP_full_time_taken = o_params["sp_full_time_taken_total"],
            o_s_SP_full_time_taken = o_s_params["sp_full_time_taken_total"],
            o_sc_SP_full_time_taken = o_sc_params["sp_full_time_taken_total"],
            o_sca_SP_full_time_taken = o_sca_params["sp_full_time_taken_total"],
            o_ss_SP_full_time_taken = o_ss_params["sp_full_time_taken_total"],
            o_scsc_SP_full_time_taken = o_scsc_params["sp_full_time_taken_total"],
            o_scsca_SP_full_time_taken = o_scsca_params["sp_full_time_taken_total"],

            # b_IP_objective = b_IP_results["objective"],
            b_s_IP_objective = b_s_IP_results["objective"],
            b_sc_IP_objective = b_sc_IP_results["objective"],
            b_sca_IP_objective = b_sca_IP_results["objective"],
            # o_IP_objective = o_IP_results["objective"],
            o_s_IP_objective = o_s_IP_results["objective"],
            o_sc_IP_objective = o_sc_IP_results["objective"],
            o_sca_IP_objective = o_sca_IP_results["objective"],
            o_ss_IP_objective = o_ss_IP_results["objective"],
            o_scsc_IP_objective = o_scsc_IP_results["objective"],
            o_scsca_IP_objective = o_scsca_IP_results["objective"],

            # b_LP_IP_gap = b_params["LP_IP_gap"],
            b_s_LP_IP_gap = b_s_params["LP_IP_gap"],
            b_sc_LP_IP_gap = b_sc_params["LP_IP_gap"],
            b_sca_LP_IP_gap = b_sca_params["LP_IP_gap"],
            # o_LP_IP_gap = o_params["LP_IP_gap"],
            o_s_LP_IP_gap = o_s_params["LP_IP_gap"],
            o_sc_LP_IP_gap = o_sc_params["LP_IP_gap"],
            o_sca_LP_IP_gap = o_sca_params["LP_IP_gap"],
            o_ss_LP_IP_gap = o_ss_params["LP_IP_gap"],
            o_scsc_LP_IP_gap = o_scsc_params["LP_IP_gap"],
            o_scsca_LP_IP_gap = o_scsca_params["LP_IP_gap"],
        ))
    end
end

results_df = DataFrame(all_results)
results_df |> 
    x -> transform!(
        x, 
        [:b_sc_LP_objective, :o_sc_LP_objective, :o_scsc_LP_objective] 
        => ByRow((x, y, z) -> y ≤ z + 1e-6 && z ≤ x + 1e-6) 
        => :sc_LP_objective_check,
    ) |>
    x -> describe(x, :detailed) |>
    x -> show(x, allrows = true)

results_df |>
    x -> filter(r -> !r.sc_LP_objective_check, x)





include("subpath_formulation.jl")

(
    n_depots , 
    n_customers , 
    n_charging ,
    n_vehicles ,
    T , 
    B , 
    seed
) = (4, 10, 9, 6, 40000, 15000, 3)

data = generate_instance(
    ;
    n_depots = n_depots,
    n_customers = n_customers,
    n_charging = n_charging,
    n_vehicles = n_vehicles,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "grid",
    shrinkage_depots = 1.0,
    shrinkage_charging = 0.7,
    T = T,
    seed = seed,
    B = B,
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
# plot_instance(data)
(
    b_LP_results, b_IP_results, b_params, b_printlist, b_subpaths, b_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark",
)
println("benchmark - not single service\t\t\t\t$(b_params["time_taken"])\t$(b_LP_results["objective"])")

(
    b_s_LP_results, b_s_IP_results, b_s_params, b_s_printlist, b_s_subpaths, b_s_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark", 
    path_single_service = true,
)
println("benchmark - single service, NO customer labels\t\t\t\t$(b_s_params["time_taken"])\t$(b_s_LP_results["objective"])")
(
    b_sc_LP_results, b_sc_IP_results, b_sc_params, b_sc_printlist, b_sc_subpaths, b_sc_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark", 
    path_single_service = true, path_check_customers = true,
)
println("benchmark - single service, customer labels\t\t\t\t$(b_sc_params["time_taken"])\t$(b_sc_LP_results["objective"])")
(
    b_sca_LP_results, b_sca_IP_results, b_sca_params, b_sca_printlist, b_sca_subpaths, b_sca_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark", 
    path_single_service = true, path_check_customers = true, check_customers_accelerated = true
)
println("benchmark - single service, customer labels (accel)\t\t\t$(b_sca_params["time_taken"])\t$(b_sca_LP_results["objective"])")

(
    o_LP_results, o_IP_results, o_params, o_printlist, o_subpaths, o_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours",
)
println("ours - not single service\t\t\t\t\t$(o_params["time_taken"])\t$(o_LP_results["objective"])")

(
    o_s_LP_results, o_s_IP_results, o_s_params, o_s_printlist, o_s_subpaths, o_s_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true,
)
println("ours - single service, NO customer labels for subpaths\t\t\t$(o_s_params["time_taken"])\t$(o_s_LP_results["objective"])")
(
    o_sc_LP_results, o_sc_IP_results, o_sc_params, o_sc_printlist, o_sc_subpaths, o_sc_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true,
)
println("ours - single service, customer labels for subpaths\t\t\t$(o_sc_params["time_taken"])\t$(o_sc_LP_results["objective"])")

(
    o_sca_LP_results, o_sca_IP_results, o_sca_params, o_sca_printlist, o_sca_subpaths, o_sca_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true, 
    check_customers_accelerated = true,
)
println("ours - single service, customer labels for subpaths (accel)\t\t$(o_sca_params["time_taken"])\t$(o_sca_LP_results["objective"])")

(
    o_ss_LP_results, o_ss_IP_results, o_ss_params, o_ss_printlist, o_ss_subpaths, o_ss_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true,
    path_single_service = true,
)
println("ours - single service, NO customer labels for subpaths and paths\t$(o_ss_params["time_taken"])\t$(o_ss_LP_results["objective"])")
(
    o_scsc_LP_results, o_scsc_IP_results, o_scsc_params, o_scsc_printlist, o_scsc_subpaths, o_scsc_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true,
    path_single_service = true, path_check_customers = true,
)
println("ours - single service, customer labels for subpaths and paths\t\t$(o_scsc_params["time_taken"])\t$(o_scsc_LP_results["objective"])")        
(
    o_scsca_LP_results, o_scsca_IP_results, o_scsca_params, o_scsca_printlist, o_scsca_subpaths, o_scsca_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true, 
    path_single_service = true, path_check_customers = true,
    check_customers_accelerated = true,
)
println("ours - single service, customer labels for subpaths and paths (accel)\t$(o_scsca_params["time_taken"])\t$(o_scsca_LP_results["objective"])")

collect_subpath_solution_metrics!(b_LP_results, data, b_subpaths, b_charging_arcs)
collect_subpath_solution_metrics!(b_IP_results, data, b_subpaths, b_charging_arcs)
collect_subpath_solution_metrics!(b_s_LP_results, data, b_s_subpaths, b_s_charging_arcs)
collect_subpath_solution_metrics!(b_s_IP_results, data, b_s_subpaths, b_s_charging_arcs)
collect_subpath_solution_metrics!(b_sc_LP_results, data, b_sc_subpaths, b_sc_charging_arcs)
collect_subpath_solution_metrics!(b_sc_IP_results, data, b_sc_subpaths, b_sc_charging_arcs)
collect_subpath_solution_metrics!(b_sca_LP_results, data, b_sca_subpaths, b_sca_charging_arcs)
collect_subpath_solution_metrics!(b_sca_IP_results, data, b_sca_subpaths, b_sca_charging_arcs)
collect_subpath_solution_metrics!(o_LP_results, data, o_subpaths, o_charging_arcs)
collect_subpath_solution_metrics!(o_IP_results, data, o_subpaths, o_charging_arcs)
collect_subpath_solution_metrics!(o_s_LP_results, data, o_s_subpaths, o_s_charging_arcs)
collect_subpath_solution_metrics!(o_s_IP_results, data, o_s_subpaths, o_s_charging_arcs)
collect_subpath_solution_metrics!(o_sc_LP_results, data, o_sc_subpaths, o_sc_charging_arcs)
collect_subpath_solution_metrics!(o_sc_IP_results, data, o_sc_subpaths, o_sc_charging_arcs)
collect_subpath_solution_metrics!(o_sca_LP_results, data, o_sca_subpaths, o_sca_charging_arcs)
collect_subpath_solution_metrics!(o_sca_IP_results, data, o_sca_subpaths, o_sca_charging_arcs)
collect_subpath_solution_metrics!(o_ss_LP_results, data, o_ss_subpaths, o_ss_charging_arcs)
collect_subpath_solution_metrics!(o_ss_IP_results, data, o_ss_subpaths, o_ss_charging_arcs)
collect_subpath_solution_metrics!(o_scsc_LP_results, data, o_scsc_subpaths, o_scsc_charging_arcs)
collect_subpath_solution_metrics!(o_scsc_IP_results, data, o_scsc_subpaths, o_scsc_charging_arcs)
collect_subpath_solution_metrics!(o_scsca_LP_results, data, o_scsca_subpaths, o_scsca_charging_arcs)
collect_subpath_solution_metrics!(o_scsca_IP_results, data, o_scsca_subpaths, o_scsca_charging_arcs)

println("benchmark - not single service\t\t\t\t\t\t$(b_params["time_taken"])\t$(b_LP_results["objective"])")
println("benchmark - single service, NO customer labels\t\t\t\t$(b_s_params["time_taken"])\t$(b_s_LP_results["objective"])")
println("benchmark - single service, customer labels\t\t\t\t$(b_sc_params["time_taken"])\t$(b_sc_LP_results["objective"])")
println("benchmark - single service, customer labels (accel)\t\t\t$(b_sca_params["time_taken"])\t$(b_sca_LP_results["objective"])")
println("ours - not single service\t\t\t\t\t\t$(o_params["time_taken"])\t$(o_LP_results["objective"])")
println("ours - single service, NO customer labels for subpaths\t\t\t$(o_s_params["time_taken"])\t$(o_s_LP_results["objective"])")
println("ours - single service, customer labels for subpaths\t\t\t$(o_sc_params["time_taken"])\t$(o_sc_LP_results["objective"])")
println("ours - single service, customer labels for subpaths (accel)\t\t$(o_sca_params["time_taken"])\t$(o_sca_LP_results["objective"])")
println("ours - single service, NO customer labels for subpaths and paths\t$(o_ss_params["time_taken"])\t$(o_ss_LP_results["objective"])")
println("ours - single service, customer labels for subpaths and paths\t\t$(o_scsc_params["time_taken"])\t$(o_scsc_LP_results["objective"])")        
println("ours - single service, customer labels for subpaths and paths (accel)\t$(o_scsca_params["time_taken"])\t$(o_scsca_LP_results["objective"])")


val, p = b_LP_results["paths"][6]
compute_path_modified_cost(data, p, o_LP_results["κ"], o_LP_results["μ"], o_LP_results["ν"], verbose = true)

    


sum(compute_charging_arc_cost(data, a) for a in p.charging_arcs)




include("subpath_formulation.jl")

# Example with objective not matching - to investigate
(
    n_depots , 
    n_customers , 
    n_charging ,
    n_vehicles ,
    T , 
    B , 
    seed
) = (4, 20, 10, 8, 40000, 15000, 9)
data = generate_instance(
    ;
    n_depots = n_depots,
    n_customers = n_customers,
    n_charging = n_charging,
    n_vehicles = n_vehicles,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "grid",
    shrinkage_depots = 1.0,
    shrinkage_charging = 0.7,
    T = T,
    seed = seed,
    B = B,
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
# plot_instance(data)
(
    b_LP_results, b_IP_results, b_params, b_printlist, b_subpaths, b_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark",
)
println("benchmark - not single service\t\t\t\t$(b_params["time_taken"])\t$(b_LP_results["objective"])")
(
    b_s_LP_results, b_s_IP_results, b_s_params, b_s_printlist, b_s_subpaths, b_s_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark", 
    path_single_service = true,
)
println("benchmark - single service, NO customer labels\t\t\t\t$(b_s_params["time_taken"])\t$(b_s_LP_results["objective"])")
(
    b_sc_LP_results, b_sc_IP_results, b_sc_params, b_sc_printlist, b_sc_subpaths, b_sc_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark", 
    path_single_service = true, path_check_customers = true,
)
println("benchmark - single service, customer labels\t\t\t\t$(b_sc_params["time_taken"])\t$(b_sc_LP_results["objective"])")
(
    b_sca_LP_results, b_sca_IP_results, b_sca_params, b_sca_printlist, b_sca_subpaths, b_sca_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "benchmark", 
    path_single_service = true, path_check_customers = true, check_customers_accelerated = true
)
println("benchmark - single service, customer labels (accel)\t\t\t$(b_sca_params["time_taken"])\t$(b_sca_LP_results["objective"])")

# (
#     o_LP_results, o_IP_results, o_params, o_printlist, o_subpaths, o_charging_arcs
# ) = subpath_formulation_column_generation_integrated_from_paths(
#     G, data, method = "ours",
# )
(
    o_s_LP_results, o_s_IP_results, o_s_params, o_s_printlist, o_s_subpaths, o_s_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true,
)
println("ours - single service, NO customer labels for subpaths\t\t\t$(o_s_params["time_taken"])\t$(o_s_LP_results["objective"])")
(
    o_sc_LP_results, o_sc_IP_results, o_sc_params, o_sc_printlist, o_sc_subpaths, o_sc_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true,
)
println("ours - single service, customer labels for subpaths\t\t\t$(o_sc_params["time_taken"])\t$(o_sc_LP_results["objective"])")

(
    o_sca_LP_results, o_sca_IP_results, o_sca_params, o_sca_printlist, o_sca_subpaths, o_sca_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true, 
    check_customers_accelerated = true,
)
println("ours - single service, customer labels for subpaths (accel)\t\t$(o_sca_params["time_taken"])\t$(o_sca_LP_results["objective"])")

(
    o_ss_LP_results, o_ss_IP_results, o_ss_params, o_ss_printlist, o_ss_subpaths, o_ss_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true,
    path_single_service = true,
)
println("ours - single service, NO customer labels for subpaths and paths\t$(o_ss_params["time_taken"])\t$(o_ss_LP_results["objective"])")
(
    o_scsc_LP_results, o_scsc_IP_results, o_scsc_params, o_scsc_printlist, o_scsc_subpaths, o_scsc_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true,
    path_single_service = true, path_check_customers = true,
)
println("ours - single service, customer labels for subpaths and paths\t\t$(o_scsc_params["time_taken"])\t$(o_scsc_LP_results["objective"])")        
(
    o_scsca_LP_results, o_scsca_IP_results, o_scsca_params, o_scsca_printlist, o_scsca_subpaths, o_scsca_charging_arcs
) = subpath_formulation_column_generation_integrated_from_paths(
    G, data, method = "ours", 
    subpath_single_service = true, subpath_check_customers = true, 
    path_single_service = true, path_check_customers = true,
    check_customers_accelerated = true,
)
println("ours - single service, customer labels for subpaths and paths (accel)\t$(o_scsca_params["time_taken"])\t$(o_scsca_LP_results["objective"])")

# collect_subpath_solution_metrics!(b_LP_results, data, b_subpaths, b_charging_arcs)
# collect_subpath_solution_metrics!(b_IP_results, data, b_subpaths, b_charging_arcs)
collect_subpath_solution_metrics!(b_s_LP_results, data, b_s_subpaths, b_s_charging_arcs)
collect_subpath_solution_metrics!(b_s_IP_results, data, b_s_subpaths, b_s_charging_arcs)
collect_subpath_solution_metrics!(b_sc_LP_results, data, b_sc_subpaths, b_sc_charging_arcs)
collect_subpath_solution_metrics!(b_sc_IP_results, data, b_sc_subpaths, b_sc_charging_arcs)
collect_subpath_solution_metrics!(b_sca_LP_results, data, b_sca_subpaths, b_sca_charging_arcs)
collect_subpath_solution_metrics!(b_sca_IP_results, data, b_sca_subpaths, b_sca_charging_arcs)
# collect_subpath_solution_metrics!(o_LP_results, data, o_subpaths, o_charging_arcs)
# collect_subpath_solution_metrics!(o_IP_results, data, o_subpaths, o_charging_arcs)
collect_subpath_solution_metrics!(o_s_LP_results, data, o_s_subpaths, o_s_charging_arcs)
collect_subpath_solution_metrics!(o_s_IP_results, data, o_s_subpaths, o_s_charging_arcs)
collect_subpath_solution_metrics!(o_sc_LP_results, data, o_sc_subpaths, o_sc_charging_arcs)
collect_subpath_solution_metrics!(o_sc_IP_results, data, o_sc_subpaths, o_sc_charging_arcs)
collect_subpath_solution_metrics!(o_sca_LP_results, data, o_sca_subpaths, o_sca_charging_arcs)
collect_subpath_solution_metrics!(o_sca_IP_results, data, o_sca_subpaths, o_sca_charging_arcs)
collect_subpath_solution_metrics!(o_ss_LP_results, data, o_ss_subpaths, o_ss_charging_arcs)
collect_subpath_solution_metrics!(o_ss_IP_results, data, o_ss_subpaths, o_ss_charging_arcs)
collect_subpath_solution_metrics!(o_scsc_LP_results, data, o_scsc_subpaths, o_scsc_charging_arcs)
collect_subpath_solution_metrics!(o_scsc_IP_results, data, o_scsc_subpaths, o_scsc_charging_arcs)
collect_subpath_solution_metrics!(o_scsca_LP_results, data, o_scsca_subpaths, o_scsca_charging_arcs)
collect_subpath_solution_metrics!(o_scsca_IP_results, data, o_scsca_subpaths, o_scsca_charging_arcs)


println("benchmark - single service, NO customer labels\t\t\t\t$(b_s_params["time_taken"])\t$(b_s_LP_results["objective"])")
println("benchmark - single service, customer labels\t\t\t\t$(b_sc_params["time_taken"])\t$(b_sc_LP_results["objective"])")
println("benchmark - single service, customer labels (accel)\t\t\t$(b_sca_params["time_taken"])\t$(b_sca_LP_results["objective"])")
println("ours - single service, NO customer labels for subpaths\t\t\t$(o_s_params["time_taken"])\t$(o_s_LP_results["objective"])")
println("ours - single service, customer labels for subpaths\t\t\t$(o_sc_params["time_taken"])\t$(o_sc_LP_results["objective"])")
println("ours - single service, customer labels for subpaths (accel)\t\t$(o_sca_params["time_taken"])\t$(o_sca_LP_results["objective"])")
println("ours - single service, NO customer labels for subpaths and paths\t$(o_ss_params["time_taken"])\t$(o_ss_LP_results["objective"])")
println("ours - single service, customer labels for subpaths and paths\t\t$(o_scsc_params["time_taken"])\t$(o_scsc_LP_results["objective"])")        
println("ours - single service, customer labels for subpaths and paths (accel)\t$(o_scsca_params["time_taken"])\t$(o_scsca_LP_results["objective"])")


(val, p) = b_sc_LP_results["paths"][5]
subpath_modified_costs = [
    compute_subpath_modified_cost(
        data, s, b_sca_LP_results["κ"], b_sca_LP_results["μ"], b_sca_LP_results["ν"]
    )
    for s in p.subpaths
]
charging_arc_costs = [
    compute_charging_arc_cost(
        data, a)
    for a in p.charging_arcs
]
sum(subpath_modified_costs)
sum(charging_arc_costs)

















include("subpath_formulation.jl")
subpath_single_service = false
subpath_check_customers = false
path_single_service = false
path_check_customers = false
check_customers_accelerated = false

### Beginning debug of  subpath_formulation_column_generation_integrated_from_paths
some_subpaths = generate_artificial_subpaths(data)
subpath_costs = compute_subpath_costs(
    data, 
    some_subpaths,
)
subpath_service = compute_subpath_service(
    data, 
    some_subpaths,
)
some_charging_arcs = Dict{
    Tuple{
        Tuple{Int, Int, Int}, 
        Tuple{Int, Int, Int}
    }, 
    Vector{ChargingArc}
}()
charging_arc_costs = Dict{
    Tuple{
        Tuple{Int, Int, Int}, 
        Tuple{Int, Int, Int}
    }, 
    Vector{Int}
}()
charging_states_a_s = Set()
charging_states_s_a = Set()
mp_results = Dict()
params = Dict()
params["number_of_subpaths"] = [sum(length(v) for v in values(some_subpaths))]
params["number_of_charging_arcs"] = [0]
params["objective"] = Float64[]
params["κ"] = Dict{Int, Float64}[]
params["μ"] = Dict{Int, Float64}[]
params["ν"] = Vector{Float64}[]
params["lp_relaxation_solution_time_taken"] = Float64[]
params["sp_base_time_taken"] = Float64[]
params["sp_full_time_taken"] = Float64[]
params["sp_total_time_taken"] = Float64[]
params["lp_relaxation_constraint_time_taken"] = Float64[]
params["number_of_new_subpaths"] = Int[]
params["number_of_new_charging_arcs"] = Int[]
params["number_of_new_charging_states"] = Int[]
params["number_of_charging_states"] = Int[]

printlist = String[]
counter = 0
converged = false

mp_model = @suppress Model(Gurobi.Optimizer)
set_attribute(mp_model, "MIPGapAbs", 1e-3)
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
    for key in keys(some_subpaths)
        for p in 1:length(some_subpaths[key])
)
w = Dict{
    Tuple{
        Tuple{
            Tuple{Int, Int, Int}, 
            Tuple{Int, Int, Int}
        }, 
        Int
    }, 
    VariableRef
}()
@constraint(
    mp_model,
    κ[i in data["N_depots"]],
    sum(
        sum(
            z[((i,0,data["B"]),state2),p]
            for p in 1:length(some_subpaths[((i,0,data["B"]),state2)])
        )        
        for (state1, state2) in keys(some_subpaths)
            if state1[1] == i && state1[2] == 0 && state1[3] == data["B"]
    )
    == data["v_start"][findfirst(x -> (x == i), data["N_depots"])]
)

flow_conservation_exprs_s_out = Dict{Tuple{Int, Int, Int}, AffExpr}()
flow_conservation_exprs_s_in = Dict{Tuple{Int, Int, Int}, AffExpr}()
flow_conservation_exprs_a_out = Dict{Tuple{Int, Int, Int}, AffExpr}()
flow_conservation_exprs_a_in = Dict{Tuple{Int, Int, Int}, AffExpr}()
flow_conservation_constrs_a_s = Dict{Tuple{Int, Int, Int}, ConstraintRef}()
flow_conservation_constrs_s_a = Dict{Tuple{Int, Int, Int}, ConstraintRef}()

@constraint(
    mp_model,
    μ[n2 in data["N_depots"]],
    sum(
        sum(
            z[(state1, state2),p]
            for p in 1:length(some_subpaths[(state1, state2)])
        )
        for (state1, state2) in keys(some_subpaths)
            if state2[1] == n2
    ) ≥ data["v_end"][n2]
)
@constraint(
    mp_model,
    ν[j in data["N_customers"]],
    sum(
        sum(
            subpath_service[((state1, state2),j)][p] * z[(state1, state2),p]
            for p in 1:length(some_subpaths[(state1, state2)])
        )
        for (state1, state2) in keys(some_subpaths)
    ) == 1
)
@expression(
    mp_model,
    subpath_costs_expr,
    sum(
        sum(
            subpath_costs[state_pair][p] * z[state_pair,p]
            for p in 1:length(some_subpaths[state_pair])
        )
        for state_pair in keys(some_subpaths)
    )
)
@expression(
    mp_model,
    charging_arc_costs_expr,
    0,
)
@objective(mp_model, Min, subpath_costs_expr + charging_arc_costs_expr)

checkpoint_reached = false



# while (
#     !converged
#     && time_limit ≥ (time() - start_time)
# )
# begin
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
        "w" => Dict(
            (key, p) => value.(w[(key, p)])
            for (key, p) in keys(w)
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





    full_labels_result = @timed find_nondominated_paths(
        G, data, mp_results["κ"], mp_results["μ"], mp_results["ν"],
        ;
        time_windows = false,
        single_service = path_single_service,
        check_customers = path_check_customers,
    )
    full_labels_time = full_labels_result.time
    (generated_subpaths, generated_charging_arcs) = get_subpaths_charging_arcs_from_negative_path_labels(
        data, full_labels_result.value,
    )










    push!(
        params["sp_base_time_taken"],
        0.0
    )
    push!(
        params["sp_full_time_taken"],
        round(full_labels_time, digits=3)
    )
    push!(
        params["sp_total_time_taken"],
        round(full_labels_time, digits=3)
    )



    if length(generated_subpaths) == 0
        push!(params["number_of_new_subpaths"], 0)
        converged = true
    else
        push!(
            params["number_of_new_subpaths"],
            sum(length(v) for v in values(generated_subpaths))
        )
    end

    if length(generated_charging_arcs) == 0
        push!(params["number_of_new_charging_arcs"], 0)
    else
        push!(
            params["number_of_new_charging_arcs"],
            sum(length(v) for v in values(generated_charging_arcs))
        )
    end

    new_charging_states_a_s = Set()
    new_charging_states_s_a = Set()
    mp_constraint_start_time = time()
    for state_pair in keys(generated_subpaths)
        if !(state_pair in keys(some_subpaths))
            some_subpaths[state_pair] = []
            subpath_costs[state_pair] = []
            for i in 1:data["n_customers"]
                subpath_service[(state_pair, i)] = []
            end
            count = 0
        else
            count = length(some_subpaths[state_pair])
        end
        for s_new in generated_subpaths[state_pair]
            if state_pair in keys(some_subpaths)
                add = !any(isequal(s_new, s) for s in some_subpaths[state_pair])
            else
                add = true
            end
            if add
                # 1: include in some_subpaths
                push!(some_subpaths[state_pair], s_new)
                # 2: add subpath cost
                push!(
                    subpath_costs[state_pair], 
                    compute_subpath_cost(data, s_new)
                )
                # 3: add subpath service
                for i in 1:data["n_customers"]
                    push!(subpath_service[(state_pair, i)], s_new.served[i])
                end
                # 4: create variable
                count += 1
                z[(state_pair, count)] = @variable(mp_model, lower_bound = 0)
                (state1, state2) = state_pair
                # 5: modify constraints starting from depot, ending at depot, and flow conservation
                if state1[1] in data["N_depots"] && state1[2] == 0.0 && state1[3] == data["B"]
                    set_normalized_coefficient(κ[state1[1]], z[state_pair,count], 1)
                elseif state1[1] in data["N_charging"]
                    push!(new_charging_states_a_s, state1)
                    if !(state1 in keys(flow_conservation_exprs_s_out))
                        flow_conservation_exprs_s_out[state1] = @expression(mp_model, 0)
                    end
                    add_to_expression!(flow_conservation_exprs_s_out[state1], z[state_pair, count])
                end
                if state2[1] in data["N_depots"]
                    set_normalized_coefficient(μ[state2[1]], z[state_pair,count], 1)
                elseif state2[1] in data["N_charging"]
                    push!(new_charging_states_s_a, state2)
                    if !(state2 in keys(flow_conservation_exprs_s_in))
                        flow_conservation_exprs_s_in[state2] = @expression(mp_model, 0)
                    end
                    add_to_expression!(flow_conservation_exprs_s_in[state2], z[state_pair, count])
                end
                # 6: modify customer service constraints
                for l in data["N_customers"]
                    set_normalized_coefficient(ν[l], z[state_pair, count], s_new.served[i])
                end
                # 7: modify objective
                set_objective_coefficient(mp_model, z[state_pair, count], subpath_costs[state_pair][count])
            end
        end
    end

    for state_pair in keys(generated_charging_arcs)
        if !(state_pair in keys(some_charging_arcs))
            some_charging_arcs[state_pair] = []
            charging_arc_costs[state_pair] = []
            count = 0
        else
            count = length(some_charging_arcs[state_pair])
        end
        for a_new in generated_charging_arcs[state_pair]
            if state_pair in keys(some_charging_arcs)
                add = !any(isequal(a_new, a) for a in some_charging_arcs[state_pair])
            else
                add = true
            end
            if add
                # 1: include in some_charging_arcs
                push!(some_charging_arcs[state_pair], a_new)
                # 2: add charging arc cost
                push!(
                    charging_arc_costs[state_pair], 
                    compute_charging_arc_cost(data, a_new)
                )
                # 4: create variable
                count += 1
                w[(state_pair, count)] = @variable(mp_model, lower_bound = 0)
                (state1, state2) = state_pair
                # 5: modify constraints starting from depot, ending at depot, and flow conservation
                if state1[1] in data["N_charging"]
                    push!(new_charging_states_s_a, state1)
                    if !(state1 in keys(flow_conservation_exprs_a_out))
                        flow_conservation_exprs_a_out[state1] = @expression(mp_model, 0)
                    end
                    add_to_expression!(flow_conservation_exprs_a_out[state1], w[state_pair, count])
                end
                if state2[1] in data["N_charging"]
                    push!(new_charging_states_a_s, state2)
                    if !(state2 in keys(flow_conservation_exprs_a_in))
                        flow_conservation_exprs_a_in[state2] = @expression(mp_model, 0)
                    end
                    add_to_expression!(flow_conservation_exprs_a_in[state2], w[state_pair, count])
                end
                # 7: modify objective
                set_objective_coefficient(
                    mp_model, 
                    w[state_pair, count], 
                    charging_arc_costs[state_pair][count],
                )
            end
        end
    end
    
    for state in new_charging_states_a_s
        if state in charging_states_a_s
            con = pop!(flow_conservation_constrs_a_s, state)
            delete(mp_model, con)
        end
        flow_conservation_constrs_a_s[state] = @constraint(
            mp_model,
            flow_conservation_exprs_s_out[state]
            == flow_conservation_exprs_a_in[state]
        )
    end
    for state in new_charging_states_s_a
        if state in charging_states_s_a
            con = pop!(flow_conservation_constrs_s_a, state)
            delete(mp_model, con)
        end
        flow_conservation_constrs_s_a[state] = @constraint(
            mp_model,
            flow_conservation_exprs_a_out[state]
            == flow_conservation_exprs_s_in[state]
        )
    end
    union!(charging_states_s_a, new_charging_states_s_a)
    union!(charging_states_a_s, new_charging_states_a_s)
    mp_constraint_end_time = time()

    push!(
        params["number_of_new_charging_states"],
        length(new_charging_states_s_a) + length(new_charging_states_a_s)
    )
    push!(
        params["number_of_charging_states"],
        length(charging_states_s_a) + length(charging_states_a_s)
    )
    push!(
        params["number_of_subpaths"], 
        sum(length(v) for v in values(some_subpaths))
    )
    if length(some_charging_arcs) == 0
        push!(params["number_of_charging_arcs"], 0)
    else
        push!(
            params["number_of_charging_arcs"],
            sum(length(v) for v in values(some_charging_arcs))
        )
    end
    push!(
        params["lp_relaxation_constraint_time_taken"],
        round(mp_constraint_end_time - mp_constraint_start_time, digits = 3)
    )
    println(
        @sprintf(
            "Iteration %3d | %.4e | %10d | %10d | %9.3f | %9.3f | %9.3f | %9.3f | %10d | %10d \n", 
            counter,
            params["objective"][counter],
            params["number_of_subpaths"][counter],
            params["number_of_charging_arcs"][counter],
            params["lp_relaxation_solution_time_taken"][counter],
            params["sp_base_time_taken"][counter],
            params["sp_full_time_taken"][counter],
            params["lp_relaxation_constraint_time_taken"][counter],
            params["number_of_new_subpaths"][counter],
            params["number_of_new_charging_arcs"][counter],
        ),
    )
# end


p = full_labels_result.value[14][13][(196452, -3548, 196452)]
full_labels = full_labels_result.value

for starting_node in data["N_depots"]
    for ending_node in data["N_depots"]
        for (key, path) in pairs(full_labels_result.value[starting_node][ending_node])
            # if any(path.served .> 1)
            #     println("$starting_node, $ending_node, $key: $path")
            # end
            if path.cost < 0
                println("$starting_node, $ending_node, $key: \n$path")
            end
        end
    end
end






κ = mp_results["κ"]
μ = mp_results["μ"]
ν = mp_results["ν"]
single_service = false
time_windows = false
check_customers = false
### Beginning debug of find_nondominated_paths
function add_full_path_label_to_collection!(
    collection::Union{
        SortedDict{
            Tuple{Int, Int, Int}, 
            FullPathLabel, 
        },
        SortedDict{
            Tuple{Int, Int, Int, Vararg{Int}}, 
            FullPathLabel, 
        },
    },
    key::Union{
        Tuple{Int, Int, Int},
        Tuple{Int, Int, Int, Vararg{Int}},
    },
    path::FullPathLabel,
    ;
    verbose::Bool = false,
)
    added = true
    for (k, p) in pairs(collection)
        if p.cost ≤ path.cost
            if all(k .≤ key)
                added = false
                if verbose
                    println("$(key), $(path.cost) dominated by $(k), $(p.cost)")
                end
                break
            end
        end
        if path.cost ≤ p.cost
            if all(key .≤ k)
                if verbose
                    println("$(key), $(path.cost) dominates $(k), $(p.cost)")
                end
                pop!(collection, k)
            end
        end
    end
    if added
        if verbose
            println("$(key), $(path.cost) added!")
        end
        insert!(collection, key, path)
    end
    return added
end

modified_costs = data["travel_cost_coeff"] * Float64.(copy(data["c"]))
for j in data["N_customers"]
    for i in data["N_nodes"]
        modified_costs[i,j] -= ν[j]
    end
end

if check_customers
    full_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Int, Int, Int, Vararg{Int}}, 
                FullPathLabel
            }()
            for current_node in data["N_nodes"]
        )
        for starting_node in data["N_depots"]
    )
else
    full_labels = Dict(
        starting_node => Dict(
            current_node => SortedDict{
                Tuple{Int, Int, Int}, 
                FullPathLabel
            }()
            for current_node in data["N_nodes"]
        )
        for starting_node in data["N_depots"]
    )
end

if check_customers
    # label key here has the following fields:
    # 1) current minimum time T_i(min)
    # 2) negative of current max charge -B_i(max)
    # 3) difference between min time and min charge, T_i(min) - B_i(min)
    # 4) if applicable, whether i-th customer served
    key = (0, -data["B"], -data["B"], zeros(Int, data["n_customers"])...)
else
    key = (0, -data["B"], -data["B"])
end
for depot in data["N_depots"]
    full_labels[depot][depot][key] = FullPathLabel(
        0.0,
        [depot],
        Int[],
        Int[],
        0,
        0,
        data["B"],
        data["B"],
        false,
        zeros(Int, data["n_customers"]),
        false,
    )
end
unexplored_states = SortedSet(
    [
        (key..., depot, depot)
        for depot in data["N_depots"]
    ]
)

t = data["t"]
B = data["B"]
q = data["q"]
if time_windows
    α = data["α"]
    β = data["β"]
else
    α = zeros(Int, data["n_nodes"])
    β = repeat([data["T"]], data["n_nodes"])
end





p
arcs = collect(zip(p.nodes[1:end-1], p.nodes[2:end]))
nodes = p.nodes[2:end]
for j in eachindex(nodes)
    node = nodes[j]
    cost = sum(data["t"][a...] for a in arcs[1:j]) + sum(p.slacks[1:j]) + sum(p.excesses[1:j])
    inside = (
        cost
        in [k[1] for k in keys(full_labels[p.nodes[1]][node])]
    )
    println("$node, $cost, $inside")
end
include("subpath_formulation.jl")
compute_path_label_modified_cost(p, data,     κ,
μ,
ν,
;
verbose = true)

p
states = []
subpaths = []
current_subpath = Subpath(
    n_customers = data["n_customers"],
    starting_node = p.nodes[1],
    starting_time = 0.0, 
    starting_charge = data["B"],
)
i = p.nodes[1]
for (j, e, s) in zip(p.nodes[2:end], p.excesses, p.slacks)
    current_subpath.current_node = j
    push!(current_subpath.arcs, (i, j))
    current_subpath.starting_time += (e + s)
    current_subpath.starting_charge += (e + s) 
    current_subpath.current_time += (data["t"][i,j] + e + s)
    current_subpath.current_charge += (- data["q"][i,j] + e + s)
    if j in data["N_charging"]
        push!(
            states, 
            (current_subpath.starting_node, current_subpath.starting_time, current_subpath.starting_charge), 
            (current_subpath.current_node, current_subpath.current_time, current_subpath.current_charge), 
        )
        push!(subpaths, current_subpath)
        current_subpath = Subpath(
            n_customers = data["n_customers"],
            starting_node = j,
            starting_time = current_subpath.current_time, 
            starting_charge = current_subpath.current_charge,
        )
    elseif j in data["N_customers"]
        current_subpath.served[j] += 1
    end
    i = j
end
push!(
    states, 
    (current_subpath.starting_node, current_subpath.starting_time, current_subpath.starting_charge), 
    (current_subpath.current_node, current_subpath.current_time, current_subpath.current_charge), 
)
push!(subpaths, current_subpath)

subpaths

some_subpaths[((14, 0, 15000), (16, 62800, 2440))]

some_subpaths[((16, 74600, 14240), (15, 145800, 0))]

some_subpaths[((15, 154242, 8442), (13, 196452, 0))]


flow_c



first(unexplored_states)
while length(unexplored_states) > 0
# while true
# while first(unexplored_states)[1] <= 176274
# while true
# while true
    # println(length(unexplored_states))
    state = pop!(unexplored_states)
    starting_node = state[end-1]
    i = state[end]
    if !(state[1:end-2] in keys(full_labels[starting_node][i]))
        # println("uh oh")
        continue
    end
    path = full_labels[starting_node][i][state[1:end-2]]
    println("$(path.time_mincharge), $(path.time_maxcharge), $(path.charge_mincharge), $(path.charge_maxcharge)")
    for j in setdiff(outneighbors(G, i), i)
    # for j in 1:12
        # j = 13
        println("$state -> $j")
        if j in data["N_customers"] && single_service && path.served[j] > 0
            println("already served $j")
            continue
        end
        # feasibility checks
        # (1) battery
        excess = max(
            0, 
            q[i,j] - path.charge_mincharge 
        )
        # (2) time windows
        if path.time_mincharge + excess + t[i,j] > β[j]
            # println("$(path.time_mincharge), $excess, $(t[i,j]), $(β[j])")
            println("not time windows feasible")
            continue
        end
        # (3) charge interval 
        if (
            (i in data["N_charging"] && excess > max(B - path.charge_mincharge, 0))
            || 
            (!(i in data["N_charging"]) && excess > max(path.charge_maxcharge - path.charge_mincharge, 0))
        )
            # if i in data["N_charging"]
            #     println("$excess, $(B), $(path.charge_mincharge)")
            # else
            #     println("$excess, $(path.charge_maxcharge), $(path.charge_mincharge)")
            # end
            println("not charge feasible")
            continue
        end
        
        new_path = copy(path)
        push!(new_path.nodes, j)
        if j in data["N_customers"]
            new_path.served[j] += 1
        end

        push!(new_path.excesses, excess)
        new_path.time_mincharge = max(
            α[j],
            path.time_mincharge + t[i,j] + excess
        )
        if i in data["N_charging"]
            slack = max(
                # floating point accuracy
                0, 
                min(
                    new_path.time_mincharge - (path.time_mincharge + t[i,j] + excess),
                    B - (path.charge_mincharge + excess),
                )
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[j],
                max(
                    α[j],
                    path.time_mincharge + (B - path.charge_mincharge) + t[i,j],
                )
            )
        else
            slack = max(
                # floating point accuracy
                0, 
                min(
                    new_path.time_mincharge - (path.time_mincharge + t[i,j] + excess),
                    path.charge_maxcharge - (path.charge_mincharge + excess),
                )
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[j],
                max(
                    α[j],
                    path.time_maxcharge + t[i,j],
                )
            )
        end
        
        new_path.charge_mincharge = (
            path.charge_mincharge 
            + excess 
            + slack
            - q[i,j]
        )
        new_path.charge_maxcharge = (
            new_path.charge_mincharge 
            + new_path.time_maxcharge 
            - new_path.time_mincharge
        )

        new_path.cost += modified_costs[i,j]
        new_path.cost += data["charge_cost_coeff"] * (slack + excess)
        
        # println("$(new_path.time_mincharge), $(new_path.time_maxcharge), $(new_path.charge_mincharge), $(new_path.charge_maxcharge)")
        # add new_path to collection
        if check_customers
            new_key = (
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge, 
                new_path.served...,
            )
        else
            new_key = (
                new_path.time_mincharge, 
                - new_path.charge_maxcharge, 
                new_path.time_mincharge - new_path.charge_mincharge
            )
        end

        added = add_full_path_label_to_collection!(
            full_labels[starting_node][j], 
            new_key, new_path, verbose = true
        )

        if added && !(j in data["N_depots"])
            new_state = (new_key..., starting_node, j)
            if !(new_state in unexplored_states)
                # println("adding state: $(new_state)")
                push!(unexplored_states, new_state)
            end
        end
    end
    # break
end

full_labels[21][13][(80328, -12169, 80328, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1)]