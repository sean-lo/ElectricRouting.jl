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
        for filepath in glob("experiments/subpath_path/03/rundata/*.csv")]...
    ) 
    |> x -> sort(x, vcat(idnames, methodnames, [:seed]))
    |> x -> unique(x, vcat(idnames, methodnames, [:seed]))
)

tall_results_df[!, :idname] = join.(eachrow(tall_results_df[!, idnames]), "_")
tall_results_df[!, :methodname] = join.(eachrow(tall_results_df[!, methodnames]), "_")

tall_results_df.LP_objective = ifelse.(tall_results_df.LP_objective .> 1e8, Inf, tall_results_df.LP_objective)
tall_results_df.IP_objective = ifelse.(tall_results_df.IP_objective .> 1e8, Inf, tall_results_df.IP_objective)

n20_results_df = filter(
    r -> r.n_customers == 20,
    tall_results_df, 
)

(cdf = groupby(n20_results_df, vcat(:idname, :methodname)) 
    |> x -> combine(
        x, 
        nrow, 
        :time_taken => mean => :time_taken_mean,
        :LP_IP_gap => mean => :LP_IP_gap_mean,
    ) 
    |> x -> select(x, 
        :time_taken_mean, 
        :nrow, 
        :LP_IP_gap_mean, 
        :idname, 
        :methodname,
    ) 
    |> x -> sort(x, :idname)
) 
(cdf 
    |> x -> unstack(x, :methodname, :idname, :nrow) 
    |> x -> CSV.write("$(@__DIR__)/results/20_nrow.csv", x)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :time_taken_mean) 
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 5)) 
    |> x -> CSV.write("$(@__DIR__)/results/20_time_taken.csv", x)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :LP_IP_gap_mean) 
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 3)) 
    |> x -> CSV.write("$(@__DIR__)/results/20_LP_IP_gap.csv", x)
)

(pslength_metrics = n20_results_df 
    |> x -> filter(
        r -> (r.LP_objective != Inf && !r.time_limit_reached),
        x,
    )
    |> x -> groupby(x, idnames) 
    |> x -> combine(x, 
        nrow,
        :lp_weighted_mean_subpath_length 
        => (x -> mean(filter(!isnan, skipmissing(x))))
        => :lp_weighted_mean_subpath_length_mean,
        :lp_weighted_mean_path_length 
        => (x -> mean(filter(!isnan, skipmissing(x))))
        => :lp_weighted_mean_path_length_mean,
        :lp_weighted_mean_ps_length 
        => (x -> mean(filter(!isnan, skipmissing(x))))
        => :lp_weighted_mean_ps_length_mean,
    )
)
(pslength_metrics 
    |> x -> hcat(x[!, 1:end-5], round.(x[!, end-4:end], sigdigits = 3)) 
    |> x -> CSV.write("$(@__DIR__)/results/20_pslength_metrics.csv", x)
)

t4_results_df = filter(
    r -> (
        r.T == 40000
    ),
    tall_results_df
)

(cdf = groupby(t4_results_df, vcat(:idname, :methodname)) 
    |> x -> combine(
        x, 
        nrow, 
        :time_taken => mean => :time_taken_mean,
        :LP_IP_gap => mean => :LP_IP_gap_mean,
    ) 
    |> x -> select(x, 
        :time_taken_mean, 
        :nrow, 
        :LP_IP_gap_mean, 
        :idname, 
        :methodname,
    ) 
    |> x -> sort(x, :idname)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :nrow) 
    |> x -> CSV.write("$(@__DIR__)/results/4_nrow.csv", x)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :time_taken_mean) 
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 5)) 
    |> x -> CSV.write("$(@__DIR__)/results/4_time_taken.csv", x)
)
(cdf 
    |> x -> unstack(x, :methodname, :idname, :LP_IP_gap_mean) 
    |> x -> hcat(x[!, 1], round.(x[!, 2:end], sigdigits = 3)) 
    |> x -> CSV.write("$(@__DIR__)/results/4_LP_IP_gap.csv", x)
)

(pslength_metrics = t4_results_df 
    |> x -> filter(
        r -> (r.LP_objective != Inf && !r.time_limit_reached),
        x,
    )
    |> x -> groupby(x, idnames) 
    |> x -> combine(x, 
        nrow,
        :lp_weighted_mean_subpath_length 
        => (x -> mean(filter(!isnan, skipmissing(x))))
        => :lp_weighted_mean_subpath_length_mean,
        :lp_weighted_mean_path_length 
        => (x -> mean(filter(!isnan, skipmissing(x))))
        => :lp_weighted_mean_path_length_mean,
        :lp_weighted_mean_ps_length 
        => (x -> mean(filter(!isnan, skipmissing(x))))
        => :lp_weighted_mean_ps_length_mean,
    )
)
(pslength_metrics 
    |> x -> hcat(x[!, 1:end-5], round.(x[!, end-4:end], sigdigits = 3)) 
    |> x -> CSV.write("$(@__DIR__)/results/4_pslength_metrics.csv", x)
)