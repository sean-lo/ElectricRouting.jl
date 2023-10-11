include("../../../src/subpath_formulation.jl")
include("../../../src/path_formulation.jl")
include("../../../src/decomposition_heuristic.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

using JuMP, Gurobi

sleep(rand() * 30.0)
const GRB_ENV = Gurobi.Env()

const TIME_LIMIT_SEC = 3600.0

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
        # method
        # path_single_service
        # path_check_customers
        # christofides
        # ngroute
        # ngroute_alt
        ("ours", false, false, false, false, false),
        ("ours", false, false,  true, false, false),
        ("ours", false, false, false,  true, false),
        ("ours", false, false,  true,  true, false),
        ("ours", false, false, false,  true,  true),
        ("ours", false, false,  true,  true,  true),
        ("ours",  true,  true, false, false, false),
    ]
    for method_param in method_params
        (
            method, 
            path_single_service, path_check_customers,
            christofides, ngroute, ngroute_alt,
        ) = method_param
        use_time_windows = false
        CGLP_results, CGIP_results, CG_params, CG_printlist, CG_some_paths = path_formulation_column_generation(
            sample_data,
            ;
            method = method,
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
        )
        DHL_results_paths = path_formulation_decomposition_heuristic(
            sample_data,
            ;
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            use_integer_paths = false,
        )
        DHI_results_paths = path_formulation_decomposition_heuristic(
            sample_data,
            ;
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            use_integer_paths = true,
        )
        DHL_objective = compute_objective_from_path_solution(DHL_results_paths, sample_data)
        DHI_objective = compute_objective_from_path_solution(DHL_results_paths, sample_data)
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
    method = String(args_df[row_index, :method])
    path_single_service = args_df[row_index, :path_single_service]
    path_check_customers = args_df[row_index, :path_check_customers]
    christofides = args_df[row_index, :christofides]
    ngroute = args_df[row_index, :ngroute]
    ngroute_alt = args_df[row_index, :ngroute_alt]

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

    CGLP_results, CGIP_results, CG_params, CG_printlist, CG_some_paths = path_formulation_column_generation(
        data,
        ;
        method = method,
        time_windows = use_time_windows,
        path_single_service = path_single_service,
        path_check_customers = path_check_customers,
        christofides = christofides,
        ngroute = ngroute,
        ngroute_alt = ngroute_alt,
        time_limit = TIME_LIMIT_SEC,
    )
    DHL_objective = missing
    try
        DHL_results_paths = path_formulation_decomposition_heuristic(
            data,
            ;
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            use_integer_paths = false,
            time_limit = TIME_LIMIT_SEC,
        )
        DHL_objective = compute_objective_from_path_solution(DHL_results_paths, data)
    catch
        nothing
    end

    DHI_objective = missing
    try
        DHI_results_paths = path_formulation_decomposition_heuristic(
            data,
            ;
            time_windows = use_time_windows,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            use_integer_paths = true,
            time_limit = TIME_LIMIT_SEC,
        )
        DHI_objective = compute_objective_from_path_solution(DHI_results_paths, data)
    catch
        nothing
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
            method = method,
            path_single_service = path_single_service,
            path_check_customers = path_check_customers,
            christofides = christofides,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            LP_objective = CGLP_results["objective"],
            IP_objective = CGIP_results["objective"],
            LP_IP_gap = CG_params["LP_IP_gap"],
            counter = CG_params["counter"],
            converged = CG_params["converged"],
            time_taken = CG_params["time_taken"],
            time_limit_reached = CG_params["time_limit_reached"],
            DHL_objective = DHL_objective,
            DHI_objective = DHI_objective,
        )
    ]
    CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
    println("Row $row_index processed: $(args_df[row_index, :])")
end