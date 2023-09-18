using DataFrames, CSV
using DelimitedFiles

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)

# results_fn = readdlm("$(@__DIR__)/out.txt")
# results_fn = vec([
#     parse(Int, x[1:end-4])
#     for x in results_fn
# ])
# results[!, :filename] = results_fn

# args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
# results[!, :ngroute] = args[results_fn, :ngroute]
# results[!, :ngroute_alt] = args[results_fn, :ngroute_alt]

# (
#     results 
#     |> x -> select!(
#         x, 
#         names(results)[1:24],
#         names(results)[65:66],
#         names(results)[25:63],
#     )
# )
# CSV.write("$(@__DIR__)/combined.csv", results)

data_fields = [
    :n_customers, 
    :seed,
]
method_fields = [
    :method, 
    :ngroute_alt, 
    :ngroute_neighborhood_charging_size, 
    :use_lmSR3_cuts,
]

(
    results 
    |> x -> sort!(
        x, 
        [
            data_fields...,
            method_fields...,
        ]
    )
)
results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100
results.LP_IP_gap_last_SR3 .*= 100

SR3_lmSR3_names = [
    "time_taken_total_SR3", 
    "sp_time_taken_mean_last_SR3", 
    "LP_IP_gap_last_SR3", 
    "converged_last_SR3", 
    "n_CG_iterations_SR3", 
]
results_SR3_lmSR3 = innerjoin(
    (
        results
        |> x -> unique(
            x, 
            vcat(data_fields, setdiff(method_fields, [:use_lmSR3_cuts])),
        )
        |> x -> select(
            x, 
            Not(SR3_lmSR3_names)
        )
    ),
    [
        (
            results 
            |> x -> unstack(
                x, 
                vcat(data_fields, setdiff(method_fields, [:use_lmSR3_cuts])),
                :use_lmSR3_cuts,
                name,
            ) 
            |> x -> rename(
                x, 
                "false" => name,
                "true" => replace(name, "SR3" => "lmSR3"),
            )
        )
        for name in SR3_lmSR3_names
    ]...,
    ;
    on = vcat(data_fields, setdiff(method_fields, [:use_lmSR3_cuts])),
)



summary = (
    results_SR3_lmSR3
    |> x -> groupby(
        x, 
        vcat(setdiff(data_fields, [:seed]), setdiff(method_fields, [:use_lmSR3_cuts]))
    )
    |> x -> combine(
        x, 
        :time_taken_first => geomean,
        :time_taken_total => geomean,
        :time_taken_total_SR3 => geomean,
        :time_taken_total_SR3 => geomean ∘ skipmissing,
        :time_taken_total_lmSR3 => geomean,
        :sp_time_taken_mean_first => geomean,
        :sp_time_taken_mean_last => geomean,
        :sp_time_taken_mean_last_SR3 => geomean,
        :sp_time_taken_mean_last_SR3 => geomean ∘ skipmissing,
        :sp_time_taken_mean_last_lmSR3 => geomean,
        :LP_IP_gap_first => mean,
        :LP_IP_gap_last => mean,
        :LP_IP_gap_last_SR3 => mean,
        :LP_IP_gap_last_SR3 => mean ∘ skipmissing,
        :LP_IP_gap_last_lmSR3 => mean,
        :converged_first => mean,
        :converged_last => mean,
        :converged_last_SR3 => mean,
        :converged_last_lmSR3 => mean,
        :n_iterations => mean,
        :n_CG_iterations => mean,
        :n_CG_iterations_SR3 => mean,
        :n_CG_iterations_lmSR3 => mean,
    )
    # |> x -> filter(
    #     r -> (
    #         r.ngroute_alt == true
    #         && r.ngroute_neighborhood_charging_size == "small"
    #     ),
    #     x,
    # )
    # |> x -> select(
    #     x, 
    #     Not([:ngroute_alt, :ngroute_neighborhood_charging_size])
    # )
)
CSV.write("$(@__DIR__)/summary.csv", summary)