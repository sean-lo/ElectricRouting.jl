using DataFrames, CSV
using DelimitedFiles
using StatsBase
using ColorSchemes
using Plots

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)

# merge individual CSVs in results folder
begin
    args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
    args.index = collect(1:nrow(args))
    all_dfs = DataFrame[]
    for ind in 1:nrow(args)
        data = CSV.read("$(@__DIR__)/results/$ind.csv", DataFrame)
        data.instance .= ind
        data.iteration .= collect(1:nrow(data))
        push!(all_dfs, data)
    end
    all_data = vcat(all_dfs...)
    select!(all_data, [:instance, :iteration], Not([:instance, :iteration]))
    CSV.write("$(@__DIR__)/all_results.csv", all_data)
end

data_fields = [
    :n_customers,
    :xmax,
    :T,
    # :density,
]
method_fields = [
    :method,
    :ngroute_neighborhood_charging_size,
]

(
    results 
    |> x -> sort!(
        x, 
        vcat(data_fields,  [:seed], method_fields),
    )
)

results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100

summary = (
    results 
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x, 
        nrow,
        :time_taken_first => geomean => :time_taken_first,
        :time_taken_total => geomean => :time_taken_total,
        :sp_time_taken_mean_first => geomean => :sp_time_taken_mean_first,
        :sp_time_taken_mean_last => geomean => :sp_time_taken_mean_last,
        :LP_IP_gap_first => mean => :LP_IP_gap_first,
        :LP_IP_gap_last => mean => :LP_IP_gap_last,
        :neighborhood_size_mean_first => mean => :neighborhood_size_mean_first,
        :neighborhood_size_mean_last => mean => :neighborhood_size_mean_last,
    )
)

# Plotting: aggregate 
begin
    μ = 5
    B = 15000
    density_range = collect(2.5:0.5:4.5)
    n_customers_range = Int.(density_range * 8)
    # TB_range = collect(4.0:0.5:5.0)
    TB_range = collect(4.0:0.5:4.0)
    
    colors = Dict(
        n_customers => c
        for (n_customers, c) in zip(
            n_customers_range[end:-1:1],
            get(ColorSchemes.viridis, collect(0:1:(length(n_customers_range)-1)) ./(length(n_customers_range)-1))
        )
    )

    plot(
        xlabel = "Time taken (s)", 
        ylabel = "Optimality gap (%)",
        xscale = :log10,
        ylim = (0.0, 100.0),
    )
    
    for n_customers in n_customers_range,
        TB in TB_range
        println("$n_customers, $TB")
        data = (
            summary 
            |> x -> filter(
                r -> (
                    r.ngroute_neighborhood_charging_size == "small"
                    && r.n_customers == n_customers
                    && r.T == Int(TB * 15000 * (μ + 1)/μ)
                ),
                x
            )
        )
        plot!(
            (
                data 
                |> x -> collect(filter(
                    r -> (
                        r.method == "benchmark"
                    ),
                    x
                )[1, [:time_taken_first, :time_taken_total]])
            ),
            (
                data 
                |> x -> collect(filter(
                    r -> (
                        r.method == "benchmark"
                    ),
                    x
                )[1, [:LP_IP_gap_first, :LP_IP_gap_last]])
            ),
            shape = :circle,
            label = "$n_customers customers, T/B = $TB",
            color = colors[n_customers],
        )
        plot!(
            (
                data 
                |> x -> collect(filter(
                    r -> (
                        r.method == "ours"
                    ),
                    x
                )[1, [:time_taken_first, :time_taken_total]])
            ),
            (
                data 
                |> x -> collect(filter(
                    r -> (
                        r.method == "ours"
                    ),
                    x
                )[1, [:LP_IP_gap_first, :LP_IP_gap_last]])
            ),
            shape = :square,
            label = "$n_customers customers, T/B = $TB, ours",
            color = colors[n_customers],
        )
    end
    plot!(legend = :topright, legendfontsize = 7)
end

# Plotting: granular


all_results = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
all_results = all_results |>
    x -> groupby(x, :instance) |>
    x -> transform(
        x, 
        :CG_time_taken => cumsum,
        :CG_LP_IP_gap => (x -> x .* 100) => :CG_LP_IP_gap,
    )

begin
    xmax = 4.0
    k = 4.0
    density = 2.5
    method = "ours"
    ngroute_neighborhood_charging_size = "small"

    B = 15000
    μ = 5
    args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
    args.index = collect(1:nrow(args))
    instances = filter(
        r -> (
            r.xmax == xmax
            # && r.k == k
            && r.T == Int(B * k * (μ + 1) / μ)
            && r.density == density
            && r.method == method
            && r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
        ),
        args
    )[!, :index]

    plot(
        xlabel = "Time taken (s)",
        # xscale = :log10,
        ylabel = "Objective gap (%)",
        ylim = (-2.0, 50.0),
        xlim = (0.0, 3.5),
    )
    for ind in instances
        data = all_results |> 
            x -> filter(
                r -> r.instance == ind,
                all_results,
            )
        plot!(
            data.CG_time_taken_cumsum,
            data.CG_LP_IP_gap,
            shape = :circle,
            # color = :black,
            label = false,
        )
    end
    plot!(legend = :topright)
end