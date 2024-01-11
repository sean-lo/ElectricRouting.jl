using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using StatsPlots
using ColorSchemes
using CairoMakie

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
names(results)

args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args.index = collect(1:nrow(args))

# # merge individual CSVs in results folder
# begin
#     all_dfs = DataFrame[]
#     for ind in 1:nrow(args)
#         fp = "$(@__DIR__)/results/$ind.csv"
#         !isfile(fp) && continue
#         data = CSV.read(fp, DataFrame)
#         data.instance .= ind
#         data.iteration .= collect(1:nrow(data))
#         push!(all_dfs, data)
#     end
#     all_data = vcat(all_dfs...)
#     select!(all_data, [:instance, :iteration], Not([:instance, :iteration]))
#     CSV.write("$(@__DIR__)/all_results.csv", all_data)
# end

data_fields = [
    :n_customers,
    :xmax,
    :T,
    # :density,
]
method_fields = [
    :method,
    :ngroute_neighborhood_size,
    :ngroute_neighborhood_charging_size,
]


(
    results 
    |> x -> sort!(
        x, 
        vcat(data_fields, [:seed], method_fields),
    )
)

results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100

results.LP_IP_gap_first .= max.(results.LP_IP_gap_first, 0)
results.LP_IP_gap_last .= max.(results.LP_IP_gap_last, 0)

all_results = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
all_results.CG_LP_IP_gap .*= 100
all_results.CG_LP_IP_gap .= ifelse.(all_results.CG_LP_IP_gap .≤ 1e-10, 0.0, all_results.CG_LP_IP_gap)

all_results = outerjoin(all_results, args, on = :instance => :index, makeunique = true)
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
        :CG_LP_IP_gap => first => :LP_IP_gap_first,
        :CG_LP_IP_gap => last => :LP_IP_gap_last,
    )
    |> x -> outerjoin(x, args, on = :instance => :index)
    |> x -> groupby(x, vcat(data_fields, method_args_fields))
    |> x -> combine(
        x, 
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


(
    all_summary
    |> x -> filter!(
        r -> (
            r.ngroute_neighborhood_charging_size == "small"
            # r.ngroute_neighborhood_charging_size == "medium"
            || 
            r.ngroute_neighborhood_size_string == "all"
            # r.ngroute_neighborhood_size_string == "1"
        ),
        x
    )
)



n_customers_range = unique(all_summary.n_customers)
colors = Dict(
    n_customers => c
    for (n_customers, c) in zip(
        n_customers_range[end:-1:1],
        get(ColorSchemes.viridis, collect(0:1:(length(n_customers_range)-1)) ./(length(n_customers_range)-1))
    )
)

for ngroute_neighborhood_charging_size in ["small", "medium"]
    for varname in ["LP_IP_gap", "LP_objective"]
        p = Plots.plot(
            figsize = (500, 400),
            xscale = :log10,
            format = :png,
        )
        Plots.hline!(
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
            for ngroute_neighborhood_size_string in ["1", "cbrt", "sqrt", "third", "all"]
                data1 = filter(r -> r.ngroute_neighborhood_size_string == ngroute_neighborhood_size_string, data)
                Plots.plot!(
                    [data1[1,:CG_time_taken_first], data1[1,:CG_time_taken_last]],
                    [data1[1, varname * "_first"], data1[1, varname * "_last"]],
                    label = false,
                    shape = [:none, :dtriangle],
                    color = colors[n_customers],
                )
            end
            Plots.plot!(
                data.CG_time_taken_first,
                data[:, varname * "_first"],
                shape = :circ,
                style = :dash,
                color = colors[n_customers],
                label = "$n_customers customers"
            )
        end
        savefig(p, "$(@__DIR__)/plots/$(varname)_time_$(ngroute_neighborhood_charging_size)_pareto.png")
        savefig(p, "$(@__DIR__)/plots/$(varname)_time_$(ngroute_neighborhood_charging_size)_pareto.pdf")
    end
end