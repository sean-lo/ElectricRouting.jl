using CSV
using DataFrames
using Glob

seed_range = collect(1:20)
data_params = [
    # n_depots
    # n_customers
    # n_charging
    # n_vehicles
    # T
    # B
    (4, 16, 9, 8, 40000, 15000),
    (4, 20, 9, 8, 40000, 15000),
    (4, 24, 9, 8, 40000, 15000),
    (4, 28, 9, 8, 40000, 15000),
]
setting_params = [
    # load, 
    # time windows
    (false, false),
    (false, true),
]
method_params = [
    # formulation
    # method
    # subpath_single_service
    # subpath_check_customers
    # path_single_service
    # path_check_customers
    # check_customers_accelerated
    ("path", "benchmark", false, false, false, false, false),
    ("path", "benchmark", false, false, true, false, false),
    ("path", "benchmark", false, false, true, true, false),
    ("path", "ours", false, false, false, false, false),
    ("path", "ours", true, false, false, false, false),
    ("path", "ours", true, true, false, false, false),
    ("path", "ours", true, false, true, false, false),
    ("path", "ours", true, true, true, true, false),
    ("subpath", "benchmark", false, false, false, false, false),
    ("subpath", "benchmark", false, false, true, false, false),
    ("subpath", "benchmark", false, false, true, true, false),
    ("subpath", "benchmark", false, false, true, true, true),
    ("subpath", "ours", false, false, false, false, false),
    ("subpath", "ours", true, false, false, false, false),
    ("subpath", "ours", true, true, false, false, false),
    ("subpath", "ours", true, true, false, false, true),
    ("subpath", "ours", true, false, true, false, false),
    ("subpath", "ours", true, true, true, true, false),
    ("subpath", "ours", true, true, true, true, true),
]

args_df = DataFrame(
    n_depots = Int[],
    n_customers = Int[],
    n_charging = Int[],
    n_vehicles = Int[],
    depot_pattern = String[], 
    customer_pattern = String[],
    charging_pattern = String[],
    shrinkage_depots = Float64[],
    shrinkage_charging = Float64[],
    T = Int[],
    B = Int[],
    seed = Int[],
    μ = Int[],
    travel_cost_coeff = Int[],
    charge_cost_coeff = Int[],
    load_scale = Float64[],
    load_shape = Float64[],
    load_tolerance = Float64[],
    batch = Int[],
    permissiveness = Float64[],
    use_load = Bool[],
    use_time_windows = Bool[],
    formulation = String[],
    method = String[],
    subpath_single_service = Bool[],
    subpath_check_customers = Bool[],
    path_single_service = Bool[],
    path_check_customers = Bool[],
    check_customers_accelerated = Bool[],
)
for data_param in data_params, seed in seed_range
    for method_param in method_params, setting_param in setting_params
        if setting_param[2] == true && method_param[2] == "ours" 
            # time_windows not compatible with "ours" (for now)
            continue
        end
        push!(args_df, 
            (
                data_param[1], data_param[2], data_param[3], data_param[4],
                "circular", "random_box", "grid",
                1.0, 0.7, 
                data_param[5], data_param[6], 
                seed,
                5, 7, 3, 
                5.0, 20.0, 1.3,
                1, 0.2,
                setting_param...,
                method_param...
            )
        )
    end
end
# results_df = vcat(
#     [
#         CSV.read(filepath, DataFrame)
#         for filepath in glob("$(@__DIR__)/combined_*.csv")
#     ]...
# ) |>
#     x -> select(
#         x, 
#         names(args_df)
#     )
# new_args_df = antijoin(
#     args_df, 
#     results_df, 
#     on = names(results_df)
# )
new_args_df = args_df
CSV.write("$(@__DIR__)/args.csv", new_args_df)