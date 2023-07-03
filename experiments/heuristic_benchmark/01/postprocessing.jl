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

tall_results_df = vcat(
    [CSV.read(filepath, DataFrame)
    for filepath in glob("experiments/heuristic_benchmark/01/rundata/*.csv")]...
) |>
    x -> sort(x, vcat(idnames, methodnames, [:seed]))

tall_results_df[!, :idname] = join.(eachrow(tall_results_df[!, idnames]), "_")
tall_results_df[!, :methodname] = join.(eachrow(tall_results_df[!, methodnames]), "_")

