using CSV
using DataFrames
using Glob

seed_range = collect(1:10)

# data_params
(
    n_depots,
    depot_pattern,
    customer_pattern,
    charging_pattern,
    customer_spread,
    xmin,
    ymin,
    ymax,
    B,
    μ,
    travel_cost_coeff,
    charge_cost_coeff,
    load_scale,
    load_shape,
    load_tolerance,
    batch,
    permissiveness,
) = (
    4, "grid", "random_box", "grid_clipped", 
    0.1, 0.0, 0.0, 2.0,
    15000, 5,
    7, 3, 
    5.0, 20.0, 1.3,
    1, 0.2,
)

xmax_range = [1.0, 2.0, 3.0]
n_customers_range = [12, 15, 18, 21]
k_range = [
    1.5, 2.0, 2.5, 3.0
]
# n_charging = Int((xmax - xmin + 1)*(ymax - ymin + 1) - 4)
k = 1.5
T = Int(B * k * (μ + 1) / μ)
# n_vehicles = 6

setting_params = [
    # load
    # time windows
    (false, false),
]
method_params = [
    # method
    # elementary
    # ngroute
    # ngroute_alt
    ("benchmark",  true, false, false,),
    ("benchmark", false, false, false,),
    # ("benchmark", false,  true, false,),
    ("benchmark", false,  true,  true,),
    ("ours",  true, false, false,),
    ("ours", false, false, false,),
    # ("ours", false,  true, false,),
    ("ours", false,  true,  true,),
]

args_df = DataFrame(
    n_depots = Int[],
    n_customers = Int[],
    n_charging = Int[],
    depot_pattern = String[], 
    customer_pattern = String[],
    charging_pattern = String[],
    customer_spread = Float64[],
    xmin = Float64[],
    xmax = Float64[],
    ymin = Float64[],
    ymax = Float64[],
    n_vehicles = Int[],
    T = Int[],
    B = Int[],
    μ = Int[],
    seed = Int[],
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
    elementary = Bool[],
    ngroute = Bool[],
    ngroute_alt = Bool[],
)
for xmax in xmax_range,
    n_customers in n_customers_range,
    k in k_range,
    method_param in method_params, 
    setting_param in setting_params,
    seed in seed_range
    if setting_param[2] == true && method_param[1] == "ours" 
        # time_windows not compatible with "ours" (for now)
        continue
    end
    n_charging = Int((xmax - xmin + 1)*(ymax - ymin + 1) - 4)
    T = Int(B * k * (μ + 1) / μ)
    n_vehicles = 6
    push!(args_df, 
        (
            n_depots,
            n_customers,
            n_charging,
            depot_pattern,
            customer_pattern,
            charging_pattern,
            customer_spread,
            xmin,
            xmax,
            ymin,
            ymax,
            n_vehicles,
            T,
            B,
            μ,
            seed,
            travel_cost_coeff,
            charge_cost_coeff,
            load_scale,
            load_shape,
            load_tolerance,
            batch,
            permissiveness,
            setting_param...,
            method_param...,
        )
    )
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
    r -> (
        r.seed == 1
        && r.xmax == 1.0
        && r.T == 27000
        && (
            r.elementary == true
            || r.n_customers == 12
        )
    ),
    args_df
)
CSV.write("$(@__DIR__)/test_args.csv", test_args_df)