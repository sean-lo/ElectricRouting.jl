using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using StatsPlots
using ColorSchemes
using Plots
using CairoMakie

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
results.density .= results.n_customers ./ (results.xmax .* 2)
data_fields = [
    :n_customers,
    :xmax,
    :T,
    :density,
]
method_fields = [
    :method,
    :ngroute_neighborhood_size,
    :ngroute_neighborhood_charging_size,
]
(
    results 
    |> x -> select!(
        x, 
        :density, 
        Not(:density),
    )
    |> x -> sort!(
        x, 
        vcat(data_fields, [:seed], method_fields),
    )
)
CSV.write("$(@__DIR__)/combined.csv", results)

# merge individual CSVs in results folder
begin
    args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
    args.index = collect(1:nrow(args))
    all_dfs = DataFrame[]
    for ind in 1:nrow(args)
        fp = "$(@__DIR__)/results/$ind.csv"
        !isfile(fp) && continue
        data = CSV.read(fp, DataFrame)
        data.instance .= ind
        data.iteration .= collect(1:nrow(data))
        push!(all_dfs, data)
    end
    all_data = vcat(all_dfs...)
    select!(all_data, [:instance, :iteration], Not([:instance, :iteration]))
    CSV.write("$(@__DIR__)/all_results.csv", all_data)
end

# Note: here, these LP_IP_gaps are at least zero up to a small floating point tolerance
results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100
results.LP_IP_gap_first .= ifelse.(results.LP_IP_gap_first .≤ 1e-10, 0.0, results.LP_IP_gap_first)
results.LP_IP_gap_last .= ifelse.(results.LP_IP_gap_last .≤ 1e-10, 0.0, results.LP_IP_gap_last)

all_results = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
all_results.CG_LP_IP_gap .*= 100
all_results.CG_LP_IP_gap .= ifelse.(all_results.CG_LP_IP_gap .≤ 1e-10, 0.0, all_results.CG_LP_IP_gap)

args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args.index = collect(1:nrow(args))
leftjoin!(
    all_results, 
    args, 
    on = :instance => :index, 
    makeunique = true,
)
select!(all_results, Not(:method_1))

method_args_fields = [
    :method,
    :ngroute_neighborhood_size_string,
    :ngroute_neighborhood_charging_size,
]

all_summary = (
    all_results
    |> x -> groupby(x, :instance)
    |> x -> combine(
        x, 
        :CG_time_taken => first => :CG_time_taken_first,
        :CG_time_taken => sum => :CG_time_taken_last,
        :CGLP_objective => (x -> first(x) / last(x)) => :LP_objective_first,
        :CGLP_objective => (x -> 1.0) => :LP_objective_last,
        :CGLP_artificial => last,
        :CGIP_artificial => last,
        :CG_LP_IP_gap => first => :LP_IP_gap_first,
        :CG_LP_IP_gap => last => :LP_IP_gap_last,
    )
    |> x -> transform(
        x,
        [:CGLP_artificial_last, :CGLP_artificial_last] => ((x, y) -> x .|| y) => :artificial 
    )
    |> x -> filter(
        r -> !r.artificial,
        x
    )
    |> x -> innerjoin(x, args, on = :instance => :index)
    |> x -> groupby(x, vcat(data_fields, method_args_fields))
    |> x -> combine(
        x, 
        nrow,
        :CG_time_taken_first => geomean => :CG_time_taken_first,
        :CG_time_taken_last => geomean => :CG_time_taken_last,
        :LP_objective_first => geomean => :LP_objective_first,
        :LP_objective_last => geomean => :LP_objective_last,
        :LP_IP_gap_first => (x -> geomean(x .+ 100) - 100) => :LP_IP_gap_first,
        :LP_IP_gap_last => (x -> geomean(x .+ 100) - 100) => :LP_IP_gap_last,
    )
    |> x -> sort!(
        x, :n_customers,
    )
)

# Plotting of pareto frontier
begin
    n_customers_range = 20:4:32
    # n_customers_range = unique(all_summary.n_customers)
    colors = Dict(
        n_customers => c
        for (n_customers, c) in zip(
            n_customers_range[end:-1:1],
            get(ColorSchemes.viridis, collect(0:1:(length(n_customers_range)-1)) ./(length(n_customers_range)-1))
        )
    )
    ngroute_neighborhood_size_string_range = ["small", "medium"]

    for ngroute_neighborhood_charging_size in ngroute_neighborhood_size_string_range
        for varname in ["LP_IP_gap", "LP_objective"]
            if varname == "LP_IP_gap"
                ylims = (-5, 100)
                ylabel = "Optimality gap (%)"
            else
                ylims = (0.82, 1.01)
                ylabel = "Relaxation bound (relative)"
            end
            xtickvalues = [1, 10, 100, 1000, 3600]
            p = Plots.plot(
                figsize = (500, 400),
                xscale = :log10,
                xlabel = "Computational time (s)",
                xticks = 10.0.^(0:4),
                ylim = ylims,
                ylabel = ylabel,
                format = :png,
            )
            vline!([3600], linestyle = :dash, label = false, color = :black)
            Plots.xticks!(p, xtickvalues, string.(xtickvalues))
            Plots.hline!(
                [0.0],
                [varname == "LP_IP_gap" ? 0.0 : 1.0],
                color = :black,
                style = :dash,
                label = false,
            )
            for n_customers in n_customers_range
                data = filter(
                    r -> (
                        r.n_customers == n_customers
                        && (
                            r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
                            || r.ngroute_neighborhood_size_string in ["1", "all"]
                        )
                    ),
                    all_summary
                )
                for ngroute_neighborhood_size_string in ["1", "cbrt", "sqrt", "third"]
                    data1 = filter(r -> r.ngroute_neighborhood_size_string == ngroute_neighborhood_size_string, data)
                    Plots.plot!(
                        [data1[1,:CG_time_taken_first], data1[1,:CG_time_taken_last]],
                        [data1[1, varname * "_first"], data1[1, varname * "_last"]],
                        label = false, 
                        linewidth = 1.6,
                        alpha = 0.8,
                        shape = [:none, :dtriangle],
                        color = colors[n_customers],
                    )
                end
                Plots.plot!(
                    data.CG_time_taken_first,
                    data[:, varname * "_first"],
                    linewidth = 1.6, alpha = 0.8,
                    style = :dash, color = colors[n_customers],
                    shape = :circ, label = "$n_customers tasks",
                )
                # square for P_none
                Plots.scatter!(
                    [data.CG_time_taken_first[1]],
                    [data[1, varname * "_first"]],
                    color = colors[n_customers],
                    shape = :square, markersize = 8, label = false,
                )
                # star for P_elem
                Plots.scatter!(
                    [data.CG_time_taken_first[end]],
                    [data[end, varname * "_first"]],
                    color = colors[n_customers],
                    shape = :star, markersize = 12, label = false,
                )
            end
            Plots.scatter!(
                [4000], [-100],
                color = :grey,
                shape = :square, markersize = 8, label = "No elementarity",
            )
            Plots.scatter!(
                [4000], [-100],
                color = :grey,
                shape = :star, markersize = 12, label = "Full elementarity",
            )
            Plots.plot!(
                labelfontsize = 12,
                legendfontsize = 11,
                tickfontsize = 11,
            )
            savefig(p, "$(@__DIR__)/plots/$(varname)_time_$(ngroute_neighborhood_charging_size)_pareto.png")
            savefig(p, "$(@__DIR__)/plots/$(varname)_time_$(ngroute_neighborhood_charging_size)_pareto.pdf")
        end
    end
end