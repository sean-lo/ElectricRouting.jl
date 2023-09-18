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
    (4, 16, 7, 6, 40000, 15000),
    (4, 20, 7, 6, 40000, 15000),
    (4, 24, 7, 6, 40000, 15000),
    (4, 28, 7, 6, 40000, 15000),
]
setting_params = [
    # load, 
    # time windows
    (false, false),
]
method_params = [
    # formulation
    # method
    # ngroute
    # ngroute_alt
    # ngroute_neighborhood_depots_size
    # ngroute_neighborhood_charging_size
    # use_lmSR3_cuts
    ("path", "benchmark", true, false, "small", "small", false)
    ("path", "benchmark", true, false, "small", "small", true)
    ("path", "benchmark", true, false, "small", "medium", false)
    ("path", "benchmark", true, false, "small", "medium", true)
    ("path", "benchmark", true, true, "small", "small", false)
    ("path", "benchmark", true, true, "small", "small", true)
    ("path", "benchmark", true, true, "small", "medium", false)
    ("path", "benchmark", true, true, "small", "medium", true)
    ("path", "ours", true, false, "small", "small", false)
    ("path", "ours", true, false, "small", "small", true)
    ("path", "ours", true, false, "small", "medium", false)
    ("path", "ours", true, false, "small", "medium", true)
    ("path", "ours", true, true, "small", "small", false)
    ("path", "ours", true, true, "small", "small", true)
    ("path", "ours", true, true, "small", "medium", false)
    ("path", "ours", true, true, "small", "medium", true)
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
    ngroute = Bool[],
    ngroute_alt = Bool[],
    ngroute_neighborhood_depots_size = String[],
    ngroute_neighborhood_charging_size = String[],
    use_lmSR3_cuts = Bool[],
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
                "circular", "random_box", "circular_packing",
                1.0, 1.0, 
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

test_args_df = filter(
    x -> (x.n_customers == 16 && x.seed == 1),
    args_df
)
CSV.write("$(@__DIR__)/test_args.csv", test_args_df)