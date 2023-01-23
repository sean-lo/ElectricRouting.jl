using CSV
using DataFrames
using Plots
using StatsBase
using StatsPlots

results_filepath = "$(@__DIR__)/combined.csv"
results_df = CSV.read(results_filepath, DataFrame)
sort!(results_df, [
    :n_customers, 
    :seed,
])
results_df[!, :arc_ip_has_solution] .= .!ismissing.(results_df[!, :arc_ip_objective])
results_df[!, :arc_ip_time_limit_reached] .= (
    results_df[!, :time_limit] .< results_df[!, :arc_ip_time_taken]
)
results_df[!, :arc_ip_has_optimal_solution] .= (
    results_df[!, :arc_ip_has_solution]
    .&& .!results_df[!, :arc_ip_time_limit_reached]
)
CSV.write(results_filepath, results_df)

results_df[!, [
    :n_customers, 
    :seed, :time_limit, 
    :arc_ip_objective, :arc_ip_time_taken, 
    :arc_lp_objective, :arc_lp_time_taken,
    :arc_ip_time_limit_reached,
    :arc_ip_has_solution,
    :arc_ip_has_optimal_solution,
]] |>
    x -> filter(r -> (r.n_customers == 9), x)

results_df |>
    x -> groupby(x, [:n_customers]) |> 
    x -> combine(x, 
        nrow,
        :arc_ip_time_taken => mean,
        :arc_ip_time_taken => std,
        :arc_ip_time_limit_reached => mean => :arc_ip_time_limit_reached_prop,
        :arc_ip_has_solution => mean => :arc_ip_has_solution_prop,
        :arc_ip_has_optimal_solution => mean => :arc_ip_has_optimal_solution_prop,
    )


unstack(
    results_df, 
    [:seed],
    :n_customers,
    :arc_ip_time_taken,
) |> 
    dropmissing |> 
    x -> select(x, Not([:seed])) |> 
    x -> describe(x, :detailed)

unstack(
    results_df, 
    [:seed],
    :n_customers,
    :arc_ip_objective,
)