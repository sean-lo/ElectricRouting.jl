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
    :formulation, 
    :method, 
    :subpath_single_service, 
    :subpath_check_customers, 
    :path_single_service, 
    :path_check_customers, 
    :check_customers_accelerated, 
    :christofides, 
    :ngroute, 
    :ngroute_alt, 
    :ngroute_neighborhood_charging_depots_size, 
]

(tall_results_df 
    = vcat(
        [CSV.read(filepath, DataFrame)
        for filepath in glob("experiments/subpath_path/02/combined_*.csv")]...
    ) 
    |> x -> unique(x, vcat(idnames, methodnames, [:seed])) 
    |> x -> sort(x, vcat(idnames, methodnames, [:seed]))
)

tall_results_df[!, :idname] = join.(eachrow(tall_results_df[!, idnames]), "_")
tall_results_df[!, :methodname] = join.(eachrow(tall_results_df[!, methodnames]), "_")

tall_results_df.LP_objective = ifelse.(tall_results_df.LP_objective .> 1e8, Inf, tall_results_df.LP_objective)
tall_results_df.IP_objective = ifelse.(tall_results_df.IP_objective .> 1e8, Inf, tall_results_df.IP_objective)

(time_results_df 
    = unstack(tall_results_df, vcat(idnames, [:seed]), :methodname, :time_taken) 
    |> x -> sort(x, vcat(idnames, [:seed]))
)
(LP_objective_results_df 
    = unstack(tall_results_df, vcat(idnames, [:seed]), :methodname, :LP_objective) 
    |> x -> sort(x, vcat(idnames, [:seed]))
)
(IP_objective_results_df 
    = unstack(tall_results_df, vcat(idnames, [:seed]), :methodname, :IP_objective) 
    |> x -> sort(x, vcat(idnames, [:seed]))
)

objective_results_df = innerjoin(
    LP_objective_results_df, 
    IP_objective_results_df, 
    on = vcat(idnames, [:seed]), 
    renamecols = "_LP" => "_IP",
)

(objective_results_df[!, "LP_objective_best"] 
    = objective_results_df
    |> x -> select(x, names(x, y -> endswith(y, "_LP")))
    |> x -> select(x, 
        AsTable(:) => ByRow(
            y -> maximum([z for z in skipmissing(y) if z != Inf], init = -Inf)
        ) => "LP_objective_best"
    )[!, 1]
)
(objective_results_df[!, "IP_objective_best"] 
    = objective_results_df 
    |> x -> select(x, names(x, y -> endswith(y, "_LP")))
    |> x -> select(x, 
        AsTable(:) => ByRow(
            y -> minimum([z for z in skipmissing(y)], init = Inf)
        ) => "IP_objective_best"
    )[!, 1]
)

objective_results_df[!, ["LP_objective_best", "IP_objective_best"]]

for name in names(objective_results_df, y -> endswith(y, "_LP"))
    objective_results_df[!, name * "_ratio"] = objective_results_df[!, name] ./ objective_results_df[!, "LP_objective_best"]
end
for name in names(objective_results_df, y -> endswith(y, "_IP"))
    objective_results_df[!, name * "_ratio"] = objective_results_df[!, name] ./ objective_results_df[!, "IP_objective_best"]
end

(LP_ratio_tall_df 
    = objective_results_df 
    |> x -> select(x, idnames, :seed, names(x, y -> endswith(y, "_LP_ratio")))
    |> x -> stack(x, Not(vcat(idnames, [:seed]))) 
    |> x -> select(
        x, 
        idnames,
        :seed,
        :variable => (x -> String.(chopsuffix.(x, "_LP_ratio"))) => :methodname, 
        :value => :LP_ratio
    )
)
(IP_ratio_tall_df 
    = objective_results_df 
    |> x -> select(x, idnames, :seed, names(x, y -> endswith(y, "_IP_ratio")))
    |> x -> stack(x, Not(vcat(idnames, [:seed]))) 
    |> x -> select(
        x, 
        idnames,
        :seed,
        :variable => (x -> String.(chopsuffix.(x, "_IP_ratio"))) => :methodname, 
        :value => :IP_ratio
    )
)

tall_results_df = outerjoin(tall_results_df, LP_ratio_tall_df, IP_ratio_tall_df, on = vcat(idnames, :seed, :methodname))
tall_results_df[!, :idname] = join.(eachrow(tall_results_df[!, idnames]), "_")

(cdf 
    = groupby(tall_results_df, vcat(:idname, :methodname)) 
    |> x -> combine(
        x, 
        nrow, 
        :time_taken => mean ∘ skipmissing => :time_taken_mean,
        :LP_ratio => (x -> mean([y for y in skipmissing(x) if y != Inf])) => :LP_ratio_mean,
        :IP_ratio => (x -> mean([y for y in skipmissing(x) if y != Inf])) => :IP_ratio_mean,
    ) 
    |> x -> select(x, :time_taken_mean, :nrow, :LP_ratio_mean, :IP_ratio_mean, :idname, :methodname) 
    |> x -> sort(x, :idname)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :nrow) 
    |> x -> CSV.write("$(@__DIR__)/nrow.csv", x)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :time_taken_mean) 
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 5)) 
    |> x -> CSV.write("$(@__DIR__)/time_taken.csv", x)
)

(cdf 
    |> x -> unstack(x, :methodname, :idname, :LP_ratio_mean) 
    # |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 5)) 
    |> x -> CSV.write("$(@__DIR__)/LP_ratio.csv", x)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :IP_ratio_mean) 
    # |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 5)) 
    |> x -> CSV.write("$(@__DIR__)/IP_ratio.csv", x)
)



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