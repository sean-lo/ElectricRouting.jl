using CSV
using DataFrames

seed_range = collect(1:20)
params = [
    (2, 9, 2, 3, 1500.0, 1050.0, 3, 3600.0),
    (2, 12, 2, 3, 1900.0, 1550.0, 4, 3600.0),
    (2, 15, 2, 3, 2250.0, 1850.0, 5, 3600.0),
    (2, 18, 2, 3, 2550.0, 2050.0, 6, 3600.0),
]

args_df = DataFrame(
    n_depots = Int[],
    n_customers = Int[],
    n_charging = Int[],
    charging_repeats = Int[],
    n_vehicles = Int[],
    shrinkage_depots = Float64[],
    shrinkage_charging = Float64[],
    T = Float64[],
    T_step = Float64[],
    B = Float64[],
    B_step = Float64[],
    seed = Int[],
    μ = Float64[],
    batch = Int[],
    time_limit = Float64[],
)
for param in params, seed in seed_range
    push!(args_df, 
        (
            param[1], param[2], param[3], 1, param[4],
            1.4, 0.6, 
            param[5], 50.0, param[6], 50.0, 
            seed,
            5.0, param[7], param[8],
        )
    )
end
CSV.write("$(@__DIR__)/args.csv", args_df)