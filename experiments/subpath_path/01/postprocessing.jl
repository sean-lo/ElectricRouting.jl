using CSV
using DataFrames
using Glob
using Plots

idnames = [:n_charging, :n_customers]
methodnames = [:method, :subpath_single_service, :path_single_service, :subpath_check_customers, :path_check_customers, :check_customers_accelerated]

tall_results_df = vcat(
    [CSV.read(filepath, DataFrame)
    for filepath in glob("experiments/subpath/01/combined_*.csv")]...
) |>
    x -> sort(x, vcat(idnames, methodnames, [:seed]))

tall_results_df[!, :methodname] = string.(
    tall_results_df[!, :method], "_", 
    tall_results_df[!, :subpath_single_service], "_", 
    tall_results_df[!, :subpath_check_customers], "_", 
    tall_results_df[!, :path_single_service], "_", 
    tall_results_df[!, :path_check_customers], "_", 
    tall_results_df[!, :check_customers_accelerated]
)

time_results_df = unstack(tall_results_df, vcat(idnames, [:seed]), :methodname, :time_taken)
LP_objective_results_df = unstack(tall_results_df, vcat(idnames, [:seed]), :methodname, :LP_objective)
IP_objective_results_df = unstack(tall_results_df, vcat(idnames, [:seed]), :methodname, :IP_objective)


# Objective tests
l = collect(skipmissing(LP_objective_results_df[!, :benchmark_false_false_false_false_false] .- 1e-4 .≤ LP_objective_results_df[!, :benchmark_false_false_true_true_false]))
count(l) / length(l)
l = collect(skipmissing(LP_objective_results_df[!, :benchmark_false_false_true_true_false] .- 1e-4 .≤ LP_objective_results_df[!, :benchmark_false_false_true_false_false]))
count(l) / length(l)

l = collect(skipmissing(LP_objective_results_df[!, :benchmark_false_false_true_true_false] .- 1e-4 .≤ LP_objective_results_df[!, :benchmark_false_false_true_true_true]))
count(l) / length(l)

# TODO: investigate
l = collect(skipmissing(LP_objective_results_df[!, :benchmark_false_false_true_true_false] .≈ LP_objective_results_df[!, :benchmark_false_false_true_true_true]))
count(l) / length(l)

l = collect(skipmissing(LP_objective_results_df[!, :ours_false_false_false_false_false] .- 1e-4 .≤ LP_objective_results_df[!, :ours_true_true_false_false_false]))
count(l) / length(l)

l = collect(skipmissing(LP_objective_results_df[!, :ours_true_true_false_false_false] .- 1e-4 .≤ LP_objective_results_df[!, :ours_true_false_false_false_false]))
count(l) / length(l)

l = collect(skipmissing(LP_objective_results_df[!, :ours_true_true_false_false_false] .≈ LP_objective_results_df[!, :ours_true_true_false_false_true]))
count(l) / length(l)

l = collect(skipmissing(LP_objective_results_df[!, :ours_true_true_false_false_false] .- 1e-4 .≤ LP_objective_results_df[!, :ours_true_false_false_false_false]))
count(l) / length(l)

l = collect(skipmissing(LP_objective_results_df[!, :ours_true_true_false_false_false] .- 1e-4 .≤ LP_objective_results_df[!, :ours_true_true_true_true_false]))
count(l) / length(l)

# TODO: investigate
l = collect(skipmissing(LP_objective_results_df[!, :ours_true_true_true_true_false] .- 1e-4 .≤ LP_objective_results_df[!, :ours_true_false_true_false_false]))
count(l) / length(l)

# TODO: investigate
l = collect(skipmissing(LP_objective_results_df[!, :ours_true_true_true_true_false] .≈ LP_objective_results_df[!, :ours_true_true_true_true_true]))
count(l) / length(l)

LP_objective_results_df[!, :b_b_sc_ratio] .= LP_objective_results_df[!, :benchmark_false_false_false_false_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :b_s_b_sc_ratio] .= LP_objective_results_df[!, :benchmark_false_false_true_false_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :b_sca_b_sc_ratio] .= LP_objective_results_df[!, :benchmark_false_false_true_true_true] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_b_sc_ratio] .= LP_objective_results_df[!, :ours_false_false_false_false_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_s_b_sc_ratio] .= LP_objective_results_df[!, :ours_true_false_false_false_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_sc_b_sc_ratio] .= LP_objective_results_df[!, :ours_true_true_false_false_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_sca_b_sc_ratio] .= LP_objective_results_df[!, :ours_true_true_false_false_true] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_ss_b_sc_ratio] .= LP_objective_results_df[!, :ours_true_false_true_false_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_scsc_b_sc_ratio] .= LP_objective_results_df[!, :ours_true_true_true_true_false] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]
LP_objective_results_df[!, :o_scsca_b_sc_ratio] .= LP_objective_results_df[!, :ours_true_true_true_true_true] ./ LP_objective_results_df[!, :benchmark_false_false_true_true_false]

mgdf = groupby(tall_results_df, vcat(idnames, [:seed]))

for g in mgdf
    valdf = filter(
        r -> (
            r.method == "benchmark" 
            && r.path_single_service 
            && r.path_check_customers
            && !r.check_customers_accelerated
        ),
        g
    )
    if nrow(valdf) > 0
        val = valdf[1, :LP_objective]
        g[!, :b_sc_ratio] = g[!, :LP_objective] ./ val
    else
        g[!, :b_sc_ratio] .= missing
    end
end

# tall_results_df[!, [idnames..., methodnames..., :seed, :b_sc_ratio]] |>
#     x -> dropmissing(x, :b_sc_ratio) |>
#     x -> filter(
#         r -> (
#             r.b_sc_ratio > 1
#             && (r.subpath_check_customers || r.path_check_customers)
#         ), x)

full_gdf = groupby(tall_results_df, vcat(idnames, methodnames))
full_combine_df = combine(
    full_gdf,
    :n_depots => first => :n_depots,
    :n_vehicles => first => :n_vehicles,
    nrow, 
    :time_taken .=> [median, std],
    :b_sc_ratio => median,
    # :sp_base_time_taken_total .=> [median, std],
    # :sp_full_time_taken_total .=> [median, std],
) |>
    x -> select(x, [:n_depots, :n_customers, :n_charging, :n_vehicles], All())

gdf = groupby(full_combine_df, idnames)

gdf[1]
gdf[2]
gdf[3]
gdf[4]

gdf[4]
gdf[5]
gdf[6]
gdf[7]