using CSV
using DataFrames
using Glob
using Plots

idnames = [
    :n_depots,
    :n_charging, 
    :n_customers,
    :n_vehicles,
    :T,
    :B,
]
methodnames = [
    :method, 
    :path_single_service, 
    :path_check_customers, 
    :christofides, 
    :ngroute, 
    :ngroute_alt, 
]

(tall_results_df 
    = vcat(
        [CSV.read(filepath, DataFrame)
        for filepath in glob("experiments/heuristic_benchmark/01/rundata/*.csv")]...
    ) 
    |> x -> unique(x, vcat(idnames, methodnames, [:seed])) 
    |> x -> sort(x, vcat(idnames, methodnames, [:seed]))
)

tall_results_df[!, :idname] = join.(eachrow(tall_results_df[!, idnames]), "_")
tall_results_df[!, :methodname] = join.(eachrow(tall_results_df[!, methodnames]), "_")

(tall_results_df 
    |> x -> filter!(
        r -> r.converged,
        x
    )
)

tall_results_df[!, :LP_DHL_gap] = tall_results_df[!, :DHL_objective] ./ tall_results_df[!, :LP_objective] .- 1
tall_results_df[!, :IP_DHI_gap] = tall_results_df[!, :DHI_objective] ./ tall_results_df[!, :IP_objective] .- 1

(cdf = tall_results_df
    |> x -> filter(
        r -> r.methodname == "ours_false_false_true_false_false",
        x
    )
    |> x -> groupby(x, idnames) 
    |> x -> combine(
        x, 
        nrow, 
        :time_taken => mean => :time_taken_mean,
        :LP_IP_gap => mean => :LP_IP_gap_mean,
        :LP_DHL_gap => mean ∘ skipmissing => :LP_DHL_gap_mean,
        :IP_DHI_gap => mean ∘ skipmissing => :IP_DHI_gap_mean,
    )
)

(cdf 
    |> x -> unstack(x, :T, :n_customers, :LP_DHL_gap_mean)
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], digits = 3))
    |> x -> CSV.write("$(@__DIR__)/results/LP_DHL_gap.csv", x)
)
(cdf 
    |> x -> unstack(x, :T, :n_customers, :IP_DHI_gap_mean)
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], digits = 3))
    |> x -> CSV.write("$(@__DIR__)/results/IP_DHI_gap.csv", x)
)
