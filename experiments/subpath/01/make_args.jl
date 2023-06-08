using CSV
using DataFrames

seed_range = collect(1:20)
data_params = [
    # n_depots
    # n_customers
    # n_charging
    # n_vehicles
    # T
    # B
    (4, 20, 6, 8, 40000, 15000),
    (4, 20, 8, 8, 40000, 15000),
    (4, 20, 10, 8, 40000, 15000),
    (4, 20, 12, 8, 40000, 15000),
]
setting_params = [
    # load, 
    # time windows
    (false, false),
]
method_params = [
    # method
    # subpath_single_service
    # subpath_check_customers
    # path_single_service
    # path_check_customers
    # check_customers_accelerated
    ("benchmark", false, false, false, false, false),
    ("benchmark", false, false, true, false, false),
    ("benchmark", false, false, true, true, false),
    ("benchmark", false, false, true, true, true),
    ("ours", false, false, false, false, false),
    ("ours", true, false, false, false, false),
    ("ours", true, true, false, false, false),
    ("ours", true, true, false, false, true),
    ("ours", true, false, true, false, false),
    ("ours", true, true, true, true, false),
    ("ours", true, true, true, true, true),
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
    method = String[],
    subpath_single_service = Bool[],
    subpath_check_customers = Bool[],
    path_single_service = Bool[],
    path_check_customers = Bool[],
    check_customers_accelerated = Bool[],
)
for data_param in data_params, seed in seed_range
    for method_param in method_params, setting_param in setting_params
        push!(args_df, 
            (
                data_param[1], data_param[2], data_param[3], data_param[4],
                "circular", "random_box", "grid",
                1.0, 0.7, 
                data_param[5], data_param[6], 
                seed,
                5, 7, 3, 
                5.0, 20.0, 1.3,
                1, 0.7,
                setting_param...,
                method_param...
            )
        )
    end
end

CSV.write("$(@__DIR__)/args.csv", args_df)