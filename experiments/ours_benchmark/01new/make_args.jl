using CSV
using DataFrames
using Glob

seed_range = collect(1:20)

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

n_customers_range = [12, 15, 18, 21]
# k_range = [
    # 2.5, 3.0, 3.5, 4.0,
# ]
# xmax_range = [2.0, 3.0, 4.0,]
xmax_k_range = [
    (2.0, 2.0),
    (2.0, 2.5),
    (2.0, 3.0),
    (3.0, 3.0),
    (2.0, 3.5),
    (3.0, 3.5),
    (3.0, 4.0),
    (4.0, 4.0),
]
density_range = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

setting_params = [
    # load
    # time windows
    (false, false),
]
method_params = [
    # method
    # elementary
    # ngroute
    ("benchmark", false, false,),
    ("benchmark", false,  true,),
    ("benchmark",  true, false,),
    ("ours", false, false,),
    ("ours", false,  true,),
    ("ours",  true, false,),
]

args_df = DataFrame(
    density = Float64[],
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
)
for density in density_range,
    (xmax, k) in xmax_k_range,
    method_param in method_params, 
    setting_param in setting_params,
    seed in seed_range
    if setting_param[2] == true && method_param[1] == "ours" 
        # time_windows not compatible with "ours" (for now)
        continue
    end
    n_customers = Int(density * (xmax - xmin) * (ymax - ymin))
    n_charging = Int((xmax - xmin + 1)*(ymax - ymin + 1) - 4)
    T = Int(B * k * (μ + 1) / μ)
    n_vehicles = 6
    push!(args_df, 
        (
            density,
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

test_args_df = vcat(
    args_df |> 
        x -> filter(
            r -> (
                r.seed == 1
                && r.density == 1.5
                && r.xmax == 2.0
                && r.T == 36000
            ),
            x
        ),
    args_df |> 
        x -> filter(
            r -> (
                r.seed == 1
                && r.elementary == true
                && r.method == "benchmark"
            ),
            x
        ) |> 
        x -> unique(
            x,
            [:n_customers]
        ),
    args_df |> 
        x -> filter(
            r -> (
                r.seed == 1
                && r.elementary == true
                && r.method == "ours"
            ),
            x
        ) |> 
        x -> unique(
            x,
            [:n_customers]
        ),
)
CSV.write("$(@__DIR__)/test_args.csv", test_args_df)

tune_args_df = copy(args_df)
filter!(
    r -> (
        r.seed == 1
        && r.elementary == false
        && r.ngroute == true
        && r.method == "ours"
    ),
    tune_args_df
)
CSV.write("$(@__DIR__)/tune_args.csv", tune_args_df)
