using DataFrames, CSV
using DelimitedFiles
using StatsBase
using ColorSchemes
using Plots

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
names(results)

args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args.index = collect(1:nrow(args))
# merge individual CSVs in results folder
begin
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
    :elementary,
    :ngroute,
    :ngroute_alt,
    :ngroute_neighborhood_charging_size,
    :use_adaptive_ngroute,
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
        :converged => mean => :converged,
        :time_limit_reached => mean => :time_limit_reached,
        :time_taken_first => geomean => :time_taken_first,
        :time_taken_total => geomean => :time_taken_total,
        :sp_time_taken_mean_first => geomean => :sp_time_taken_mean_first,
        :sp_time_taken_mean_last => geomean => :sp_time_taken_mean_last,
        :LP_objective_first => mean => :LP_objective_first,
        :LP_objective_last => mean => :LP_objective_last,
        :IP_objective_first => mean => :IP_objective_first,
        :IP_objective_last => mean => :IP_objective_last,
        :LP_IP_gap_first => mean => :LP_IP_gap_first,
        :LP_IP_gap_last => mean => :LP_IP_gap_last,
        :neighborhood_size_mean_first => mean => :neighborhood_size_mean_first,
        :neighborhood_size_mean_last => mean => :neighborhood_size_mean_last,
    )
    |> x -> sort(
        x,
        vcat(data_fields, method_fields),
    )
)


# Building table

begin
    ours_ngroute_datas = DataFrame[]
    for size in ["small", "medium", "large"]
        obj_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.method == "ours"
                    && r.use_adaptive_ngroute
                    && r.ngroute_neighborhood_charging_size == size
                ),
                x,
            )
            |> x -> select(
                x, 
                vcat(data_fields, [:LP_objective_first, :LP_objective_last])
            )
            |> x -> stack(
                x, [:LP_objective_first, :LP_objective_last]
            )
            |> x -> select(
                x,
                :n_customers, 
                :T => (x -> x ./ 15000) => :TB,
                :variable => :adaptive,
                :value => String(String(:LP_objective) * "_" * size),
            )
        )
        replace!(obj_df.adaptive, "LP_objective_first" => "non-adaptive")
        replace!(obj_df.adaptive, "LP_objective_last" => "adaptive")
        push!(ours_ngroute_datas, obj_df)
        time_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.method == "ours"
                    && r.use_adaptive_ngroute
                    && r.ngroute_neighborhood_charging_size == size
                ),
                x,
            )
            |> x -> select(
                x, 
                vcat(data_fields, [:time_taken_first, :time_taken_total])
            )
            |> x -> stack(
                x, [:time_taken_first, :time_taken_total]
            )
            |> x -> select(
                x,
                :n_customers, 
                :T => (x -> x ./ 15000) => :TB,
                :variable => :adaptive,
                :value => String(String(:time_taken) * "_" * size),
            )
        )
        replace!(time_df.adaptive, "time_taken_first" => "non-adaptive")
        replace!(time_df.adaptive, "time_taken_total" => "adaptive")
        push!(ours_ngroute_datas, time_df)
    end

    ours_ngroute_data = outerjoin(
        ours_ngroute_datas...,
        on = [:n_customers, :TB, :adaptive]
    ) |> 
        x -> sort(
            x, 
            [:n_customers, :TB, order(:adaptive, rev = true)]
        )

    ours_elementary_data = (
        summary 
        |> x -> filter(
            r -> (
                r.method == "ours"
                && r.elementary
            ),
            x
        )
        |> x -> select(
            x, 
            :n_customers, 
            :T => (x -> x ./ 15000) => :TB,
            [:time_limit_reached, :LP_objective_last] => ((x, y) -> ifelse.(x .> 0.0, NaN, y)) => :LP_objective_elementary,
            [:time_limit_reached, :time_taken_total] => ((x, y) -> ifelse.(x .> 0.0, NaN, y)) => :time_taken_elementary,
        )
    )

    ours_data_table = outerjoin(ours_ngroute_data, ours_elementary_data, on = [:n_customers, :TB])

    for n_customers in 20:4:36, TB in 4.8:0.6:6.0
        df = @view ours_data_table[(ours_data_table.n_customers .== n_customers) .& (ours_data_table.TB .== TB), :]
        df[!, :LP_objective_elementary] .= df[df.adaptive .== "adaptive", :LP_objective_small]
    end

    ours_data_table.percent_gap_small .= max.(0.0, 100 .- 100 .* ours_data_table.LP_objective_small ./ ours_data_table.LP_objective_elementary)
    ours_data_table.percent_gap_medium .= max.(0.0, 100 .- 100 .* ours_data_table.LP_objective_medium ./ ours_data_table.LP_objective_elementary)
    ours_data_table.percent_gap_large .= max.(0.0, 100 .- 100 .* ours_data_table.LP_objective_large ./ ours_data_table.LP_objective_elementary)
    ours_data_table[
        2:2:nrow(ours_data_table), 
        [
            :LP_objective_small, :LP_objective_medium, :LP_objective_large, 
            :percent_gap_small, :percent_gap_medium, :percent_gap_large,
        ]
    ] .= 0.0

    (
        ours_data_table
        |> x -> transform(
            x,
            [:time_taken_small, :time_taken_medium, :time_taken_large, :time_taken_elementary]
            .=> (x -> round.(round.(x, sigdigits = 4), digits = 3)) 
            .=> [:time_taken_small, :time_taken_medium, :time_taken_large, :time_taken_elementary],
            [:LP_objective_small, :LP_objective_medium, :LP_objective_large, :LP_objective_elementary]
            .=> (x -> round.(x, sigdigits = 4) ./ 1e6) 
            .=> [:LP_objective_small, :LP_objective_medium, :LP_objective_large, :LP_objective_elementary],
            [:percent_gap_small, :percent_gap_medium, :percent_gap_large]
            .=> (x -> round.(round.(x, sigdigits = 4), digits = 3))
            .=> [:percent_gap_small, :percent_gap_medium, :percent_gap_large],
        )
        |> x -> select(
            x, 
            :n_customers, :TB, :adaptive,
            :LP_objective_elementary, :time_taken_elementary,
            :LP_objective_large, :percent_gap_large, :time_taken_large,
            :LP_objective_medium, :percent_gap_medium, :time_taken_medium,
            :LP_objective_small, :percent_gap_small, :time_taken_small,
        )
        |> x -> CSV.write("$(@__DIR__)/combined_summary_v2.csv", x)
    )
end

begin
    ours_summary = (
        summary 
        |> x -> filter(
            r -> r.method == "ours",
            x,
        )
    )
    ours_nonadaptive_summary = innerjoin(
        [
            (
                ours_summary 
                |> x -> filter(
                    r -> (
                        r.use_adaptive_ngroute
                        && r.ngroute_neighborhood_charging_size == size
                    ),
                    x,
                ) 
                |> x -> select(
                    x,
                    vcat(data_fields, [:LP_objective_first, :time_taken_first])
                )
            )
            for size in ["small", "medium", "large"]
        ]..., 
        on = data_fields, 
        makeunique = true,
    ) |>
        x -> select(
            x, 
            data_fields,
            :LP_objective_first   => :LP_objective_small,
            :time_taken_first   => :time_taken_small,
            :LP_objective_first_1 => :LP_objective_medium,
            :time_taken_first_1 => :time_taken_medium,
            :LP_objective_first_2 => :LP_objective_large,
            :time_taken_first_2 => :time_taken_large,
        )
    
    ours_adaptive_summary = innerjoin(
        [
            (
                ours_summary 
                |> x -> filter(
                    r -> (
                        r.use_adaptive_ngroute
                        && r.ngroute_neighborhood_charging_size == size
                    ),
                    x,
                ) 
                |> x -> select(
                    x,
                    vcat(data_fields, [:LP_objective_last, :time_taken_total])
                )
            )
            for size in ["small", "medium", "large"]
        ]..., 
        on = data_fields, 
        makeunique = true,
    ) |>
        x -> select(
            x, 
            data_fields,
            :LP_objective_last => :LP_objective_elementary,
            :time_taken_total   => :time_taken_adaptive_small,
            :time_taken_total_1 => :time_taken_adaptive_medium,
            :time_taken_total_2 => :time_taken_adaptive_large,
        )

    ours_elementary_summary = (
        ours_summary 
        |> x -> filter(
            r -> r.elementary,
            x,
        ) 
        |> x -> select(
            x,
            data_fields,
            [:time_limit_reached, :time_taken_total] => ((x, y) -> ifelse.(x .> 0.0, NaN, y)) => :time_taken_elementary,
        )
    )
    ours_combined_summary = (
        innerjoin(
            ours_nonadaptive_summary,
            ours_adaptive_summary,
            ours_elementary_summary,
            on = data_fields,
        )
        |> x -> transform(
            x,
            [:LP_objective_small, :LP_objective_elementary] => ((x, y) -> 100.0 .- (100.0 .* x ./ y)) => :LP_objective_gap_small,
            [:LP_objective_medium, :LP_objective_elementary] => ((x, y) -> 100.0 .- (100.0 .* x ./ y)) => :LP_objective_gap_medium,
            [:LP_objective_large, :LP_objective_elementary] => ((x, y) -> 100.0 .- (100.0 .* x ./ y)) => :LP_objective_gap_large,
        )
        |> x -> select(
            x, 
            data_fields,
            :LP_objective_small, :LP_objective_gap_small, :time_taken_small,
            :LP_objective_medium, :LP_objective_gap_medium, :time_taken_medium,
            :LP_objective_large, :LP_objective_gap_large, :time_taken_large,
            :LP_objective_elementary,
            :time_taken_adaptive_small, 
            :time_taken_adaptive_medium, 
            :time_taken_adaptive_large, 
            :time_taken_elementary, 
        )
    )
    (
        ours_combined_summary
        |> x -> select(
            x, 
            :n_customers, 
            :T, 
            Not(data_fields) .=> (x -> round.(round.(x, sigdigits = 3), digits = 3)) .=> Not(data_fields),
        )
        |> x -> CSV.write("$(@__DIR__)/combined_summary.csv", x)
    )
end


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

all_results

xmax = 4.0
k = 4.0
density = 2.5
method = "ours"
ngroute_neighborhood_charging_size = "small"
ngroute_instances = filter(
    r -> (
        r.xmax == xmax
        && !r.elementary
        && r.T == Int(B * k * (μ + 1) / μ)
        && r.density == density
        && r.method == method
        && r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
    ),
    args
)[!, :index]


begin
    xmax_k_range = [
        (4.0, 4.0),
        (4.0, 4.5),
        (4.0, 5.0),
    ]
    density_range = [
        2.5, 3.0, 3.5, 4.0, 4.5, 
    ]
    ngroute_neighborhood_charging_size_range = [
        "small", 
        "medium", 
        "large",
    ]
    method = "ours"

    B = 15000
    μ = 5

    for (xmax, k) in xmax_k_range,
        density in density_range,
        ngroute_neighborhood_charging_size in ngroute_neighborhood_charging_size_range

        ngroute_instances = filter(
            r -> (
                r.xmax == xmax 
                && !r.elementary
                && r.T == Int(B * k * (μ + 1) / μ)
                && r.density == density
                && r.method == method
                && r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
            ),
            args
        )[!, :index]
        println(ngroute_instances)
        elem_instances = filter(
            r -> (
                r.xmax == xmax 
                && r.elementary
                && r.T == Int(B * k * (μ + 1) / μ)
                && r.density == density
                && r.method == method
            ),
            args
        )[!, :index]
        (max_times 
            = all_results 
            |> x -> filter(
                r -> r.instance in ngroute_instances,
                x
            )
            |> x -> groupby(
                x, 
                :instance,
            )
            |> x -> combine(
                x, 
                :CG_time_taken => sum,
            )
            |> x -> sort(
                x, 
                :CG_time_taken_sum
            )
        )
        perm = max_times.instance
        println(perm)

        p = plot(
            xlabel = "Time taken (s)",
            # size = (550, 750),
            # xscale = :log10,
            # yscale = :log10,
            ylabel = "Objective value",
            ylim = (0.79, 1.31),
            xlim = (0.0, max_times.CG_time_taken_sum[end] * 1.05),
        )
        colors = Dict(
            ind => c
            for (ind, c) in zip(
                perm[end:-1:1],
                get(ColorSchemes.viridis, collect(0:1:(length(perm)-1)) ./(length(perm)-1))
            )
        )
        hline!(
            [1.0],
            style = :dash,
            color = :black,
            alpha = 0.7,
            label = false,
        )
        ngroute_ind = 0
        for (ngroute_ind, elem_ind) in zip(ngroute_instances, elem_instances)
            ngroute_data = all_results |> 
                x -> filter(
                    r -> r.instance == ngroute_ind,
                    x,
                )
            elem_data = all_results |> 
                x -> filter(
                    r -> r.instance == elem_ind,
                    x,
                )
            plot!(
                ngroute_data.CG_time_taken_cumsum,
                (ngroute_data.CGLP_objective ./ ngroute_data.CGLP_objective[end]),
                shape = :circle,
                markersize = 3,
                color = colors[ngroute_ind],
                label = false,
                alpha = 0.6,
            )
            plot!(
                [ngroute_data.CG_time_taken_cumsum[end], ngroute_data.CG_time_taken_cumsum[end]],
                [ngroute_data.CGLP_objective[end], ngroute_data.CGIP_objective[end]] ./ ngroute_data.CGLP_objective[end],
                style = :dash,
                color = colors[ngroute_ind],
                label = false,
                alpha = 0.9,
            )
            plot!(
                [ngroute_data.CG_time_taken_cumsum[end]],
                [ngroute_data.CGIP_objective[end] / ngroute_data.CGLP_objective[end]],
                shape = :square,
                markersize = 4,
                color = colors[ngroute_ind],
                label = false,
                alpha = 0.9,
            )
        end
        plot!(
            [-1],
            [0],
            shape = :circle,
            markersize = 3,
            alpha = 0.6,
            color = colors[perm[end]],
            label = "CG iterations (adaptive ng-route)",
        )
        plot!(
            [-1],
            [0],
            shape = :square,
            markersize = 4,
            alpha = 0.9,
            color = colors[perm[end]],
            label = "IP heuristic solution",
        )
        display(p)
        savefig(p, "$(@__DIR__)/plots/LP_IP_objective_trajectory_$(xmax)_$(k)_$(density)_$(method)_$(ngroute_neighborhood_charging_size).pdf")
        savefig(p, "$(@__DIR__)/plots/LP_IP_objective_trajectory_$(xmax)_$(k)_$(density)_$(method)_$(ngroute_neighborhood_charging_size).png")
    end
end

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
                x,
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

