using CSV
using DataFrames
using Glob

seed_range = collect(1:20)

# data_params
(
    n_vehicles,
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
    charge_cost_coeff_increment,
    load_scale,
    load_shape,
    load_tolerance,
    batch,
    permissiveness,
) = (
    10, 4, "grid", "random_box", "grid_clipped", 
    0.1, 0.0, 0.0, 2.0,
    18000, 5,
    7, 3, 1, 
    5.0, 20.0, 1.3,
    1, 0.2,
)

xmax_k_range = [
    # k is proportional to time horizon
    # xmax defines the grid area
    # (ymax fixed at 2)
    (4.0, 4.0),
]
T_range = [86400]
density_range = [
    # number of customers per unit area
    # multiply by xmax * ymax (= 8) to get number of customers
    2.5, 3.0, 3.5, 4.0, 4.5, 5.0,
]
charge_cost_nlevels_range = [
    # charge_cost_nlevels
    1, 2, 3, 4, 5,
]
setting_params = [
    # load
    # time windows
    # sparsify graph
    (false, false, false),
    (false, false, true),
]
method_params = [
    # method
    # ngroute_neighborhood_charging_size
    # use_lmSR3_cuts
    # max_SR3_cuts
    ("ours", "small", true, 5),
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
    charge_cost_coeff_increment = Int[],
    charge_cost_nlevels = Int[],

    load_scale = Float64[],
    load_shape = Float64[],
    load_tolerance = Float64[],
    batch = Int[],
    permissiveness = Float64[],

    use_load = Bool[],
    use_time_windows = Bool[],
    sparse_graph = Bool[],

    method = String[],
    ngroute_neighborhood_charging_size = String[],
    use_lmSR3_cuts = Bool[],
    max_SR3_cuts = Int[],
)
for (xmax, k) in xmax_k_range,
    density in density_range,
    charge_cost_nlevels in charge_cost_nlevels_range,
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
            charge_cost_coeff_increment,
            charge_cost_nlevels,
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
new_args_df = args_df
sort!(
    new_args_df,
    [
        order(:n_customers),
    ],
)
CSV.write("$(@__DIR__)/args.csv", new_args_df)

test_args_df = args_df |> 
    x -> filter(
        r -> (
            r.seed == 1
            && r.density == 3.0
            && r.sparse_graph
            && r.charge_cost_nlevels == maximum(charge_cost_nlevels_range)
            && r.T == minimum(T_range)
        ),
        x
    )
CSV.write("$(@__DIR__)/test_args.csv", test_args_df)