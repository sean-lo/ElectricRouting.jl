using CSV
using DataFrames
using Glob

seed_range = collect(1:20)
data_params = collect(Iterators.product(
    [4], # n_depots
    20:4:40, # n_customers
    [7], # n_charging
    [8], # n_vehicles
    40000:10000:80000, # T
    [15000], # B
))
setting_params = [
    # load, 
    # time windows
    (false, false),
]
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
    # ngroute_alt
    # ngroute_neighborhood_charging_depots_size
    ("subpath", "ours", false, false, false, false, false,  true, false, false, "none"),
    ("subpath", "ours", false, false, false, false, false,  true,  true, false, "small"),
    ("subpath", "ours", false, false, false, false, false,  true,  true,  true, "small"),
    ("subpath", "ours", false, false, false, false, false,  true,  true, false, "large"),
    ("subpath", "ours", false, false, false, false, false,  true,  true,  true, "large"),
    ("subpath", "ours",  true,  true, false, false,  true,  true, false, false, "none"),
    ("subpath", "ours",  true,  true, false, false, false,  true, false, false, "none"),
    ("subpath", "ours",  true,  true,  true,  true,  true,  true, false, false, "none"),
    ("subpath", "ours",  true,  true,  true,  true, false,  true, false, false, "none"),
    ("path", "ours", false, false, false, false, false,  true, false, false, "none"),
    ("path", "ours", false, false, false, false, false,  true,  true, false, "small"),
    ("path", "ours", false, false, false, false, false,  true,  true,  true, "small"),
    ("path", "ours", false, false, false, false, false,  true,  true, false, "large"),
    ("path", "ours", false, false, false, false, false,  true,  true,  true, "large"),
    ("path", "ours",  true,  true, false, false, false,  true, false, false, "none"),
    ("path", "ours",  true,  true,  true,  true, false,  true, false, false, "none"),
    ("subpath", "benchmark", false, false, false, false, false,  true, false, false, "none"),
    ("subpath", "benchmark", false, false, false, false, false,  true,  true, false, "small"),
    ("subpath", "benchmark", false, false, false, false, false,  true,  true,  true, "small"),
    ("subpath", "benchmark", false, false, false, false, false,  true,  true, false, "large"),
    ("subpath", "benchmark", false, false, false, false, false,  true,  true,  true, "large"),
    ("path", "benchmark", false, false, false, false, false,  true, false, false, "none"),
    ("path", "benchmark", false, false, false, false, false,  true,  true, false, "small"),
    ("path", "benchmark", false, false, false, false, false,  true,  true,  true, "small"),
    ("path", "benchmark", false, false, false, false, false,  true,  true, false, "large"),
    ("path", "benchmark", false, false, false, false, false,  true,  true,  true, "large"),
    ("subpath", "benchmark", false, false,  true,  true, false, false, false, false, "none"),
    ("subpath", "benchmark", false, false,  true,  true,  true, false, false, false, "none"),
    ("path", "benchmark", false, false,  true,  true, false, false, false, false, "none"),
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
    christofides = Bool[],
    ngroute = Bool[],
    ngroute_alt = Bool[],
    ngroute_neighborhood_charging_depots_size = String[],
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
                method_param...,
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
sort!(
    new_args_df,
    [
        order(:christofides, rev = true),
        order(:method, rev = true),
        order(:formulation, rev = true),
        order(:T),
        order(:n_customers),
        order(:ngroute, rev = true),
        order(:ngroute_neighborhood_charging_depots_size, rev = true),
        order(:ngroute_alt),
        order(:path_check_customers),
        order(:check_customers_accelerated, rev = true),
    ]
)
CSV.write("$(@__DIR__)/args.csv", new_args_df)