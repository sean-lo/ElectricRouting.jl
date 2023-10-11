using DataFrames, CSV
using DelimitedFiles

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
# results.density .= results.n_customers ./ (results.xmax .* 2)
# CSV.write("$(@__DIR__)/combined.csv", results)

data_fields = [
    :n_customers,
    :xmax,
    :T,
    :density,
]
method_fields = [
    :method,
    :elementary,
    :ngroute,
    :ngroute_alt,
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

results.LP_IP_gap .*= 100

summary = (
    results 
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x, 
        :time_taken => geomean => :time_taken,
        :sp_time_taken_mean => geomean => :sp_time_taken_mean,
        :counter => mean => :counter,
        :converged => mean => :converged,
        :time_limit_reached => mean => :time_limit_reached,
        :LP_IP_gap => mean => :LP_IP_gap,
        :mean_subpath_length => mean => :mean_subpath_length,
        :mean_path_length => mean => :mean_path_length,
        :mean_ps_length => mean => :mean_ps_length,
    )
)

(
    results
    |> x -> filter(
        r -> (
            !r.converged
            # && r.xmax == 2.0
            # && r.T == 27000
            # && r.counter > 100
        ),
        x
    )
    |> x -> select(
        x, 
        vcat(data_fields, [:seed], method_fields, :counter, :LP_IP_gap)
    )
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x, 
        nrow, 
        :counter => mean,
        :LP_IP_gap => mean,
    )
    # |> x -> unstack(
    #     x,
    #     [:n_customers, :xmax,],
    #     :T,
    #     :nrow,
    # )
    # |> x -> sort(x, :xmax)
)

filtered_summary = summary |> 
    x -> filter(
        r -> (
            !(r.xmax == 2.0 && r.T == 36000)
            # && !(r.density in [3.5, 4.0] && r.T == 72000)
        ),
        x
    )


filtered_summary.time_taken_filtered = copy(filtered_summary.time_taken)
filtered_summary.time_taken_filtered[filtered_summary.converged .!= 1.0] .= NaN
filtered_summary.method_combined = [
    join(filtered_summary[i, method_fields], "_")
    for i in 1:nrow(filtered_summary)
]

(
    filtered_summary
    |> x -> sort(
        x, 
        vcat(data_fields, [:elementary, :ngroute, :method])
    )
    |> x -> unstack(
        x, 
        data_fields, 
        :method_combined, 
        :time_taken_filtered,
    )
    |> x -> select(
        x, 
        :n_customers,
        :xmax => (x -> x .* 2.0) => :area,
        :T => (x -> x ./ 15000) => :TB,
        :benchmark_false_false_false, :ours_false_false_false,
        [:benchmark_false_false_false, :ours_false_false_false] => ((x, y) -> (1.0 .- y ./ x) .* 100) => :gap_false_false_false,
        :benchmark_false_true_true, :ours_false_true_true,
        [:benchmark_false_true_true, :ours_false_true_true] => ((x, y) -> (1.0 .- y ./ x) .* 100) => :gap_false_true_true,
        :benchmark_true_false_false, :ours_true_false_false,
        [:benchmark_true_false_false, :ours_true_false_false] => ((x, y) -> (1.0 .- y ./ x) .* 100) => :gap_true_false_false,
    )
    |> x -> select(
        x, 
        :n_customers, 
        :area,
        :TB,
        Not([:n_customers, :area, :TB]) .=> (x -> round.(round.(x, sigdigits = 3), digits = 3)) .=> Not([:n_customers, :area, :TB]),
    )
    |> x -> CSV.write("$(@__DIR__)/filtered_summary.csv", x)
)


density_range = collect(1.5:0.5:4.0)
TB_range = collect(2.5:0.5:4.0)

begin
    (elementary, ngroute, ngroute_alt) = (true, false, false)
    xmax = 3.0
    data = (
        filtered_summary 
        |> x -> filter(
            r -> (
                r.elementary == elementary
                && r.ngroute == ngroute
                && r.ngroute_alt == ngroute_alt
                # && r.method == "benchmark"
                && r.xmax == xmax
                #  && r.density == density
            ),
            x
        )
    )
    xlims = data |>
        x -> select(x, :time_taken)[:, 1] |>
        extrema
    xloglims = (floor(log10(xlims[1]))), ceil(log10(xlims[2]))
    xlims = 10.0^xloglims[1], 10.0^xloglims[2]
    xticks = 10.0.^(xloglims[1]:xloglims[2])
    ylims = data |>
        x -> select(x, :LP_IP_gap)[:, 1] |>
        extrema
    ylims = (0, ylims[2] * 1.1)    

    colors = Dict(
        density => c
        for (density, c) in zip(
            density_range[end:-1:1],
            get(ColorSchemes.viridis, collect(0:1:(length(density_range)-1)) ./(length(density_range)-1))
        )
    )
    plot(
        xlabel = "Time taken (s)",
        xscale = :log,
        xticks = xticks,
        xlim = xlims,
        ylabel = "Optimality gap (%)",
        ylim = ylims,
    )
    for density in density_range
        Plots.plot!(
            (
                filtered_summary 
                |> x -> filter(
                    r -> (
                        r.elementary == elementary
                        && r.ngroute == ngroute
                        && r.ngroute_alt == ngroute_alt
                        && r.method == "benchmark"
                        && r.xmax == xmax
                        && r.density == density
                    ),
                    x
                )[!, :time_taken]
            ),
            (
                filtered_summary 
                |> x -> filter(
                    r -> (
                        r.elementary == elementary
                        && r.ngroute == ngroute
                        && r.ngroute_alt == ngroute_alt
                        && r.method == "benchmark"
                        && r.xmax == xmax
                        && r.density == density
                    ),
                    x
                )[!, :LP_IP_gap]
            ),
            label = "density = $density",
            marker = :circle,
            color = colors[density],
        )
        Plots.plot!(
            (
                filtered_summary 
                |> x -> filter(
                    r -> (
                        r.elementary == elementary
                        && r.ngroute == ngroute
                        && r.ngroute_alt == ngroute_alt
                        && r.method == "ours"
                        && r.xmax == xmax
                        && r.density == density
                    ),
                    x
                )[!, :time_taken]
            ),
            (
                filtered_summary 
                |> x -> filter(
                    r -> (
                        r.elementary == elementary
                        && r.ngroute == ngroute
                        && r.ngroute_alt == ngroute_alt
                        && r.method == "ours"
                        && r.xmax == xmax
                        && r.density == density
                    ),
                    x
                )[!, :LP_IP_gap]
            ),
            label = "density = $density, ours",
            marker = :square,
            color = colors[density],
        )
    end
    plot!(legend = :bottomright,)
end

function make_pivottable(
    df::DataFrame,
    elementary::Bool,
    ngroute::Bool,
    ngroute_alt::Bool,
    var::Symbol, # one of :time_taken, :LP_IP_gap_mean
    val::Symbol, # one of :benchmark, :ours, :ratio, :percent_gap
)
    if !((elementary, ngroute, ngroute_alt) in [
        (true, false, false),
        (false, false, false), 
        (false, true, true),
    ])
        error()
    end
    result = (
        df 
        |> x -> filter(
            r -> (
                r.elementary == elementary
                && r.ngroute == ngroute
                && r.ngroute_alt == ngroute_alt
            ),
            x
        )
        |> x -> unstack(
            x, 
            [:n_customers, :xmax, :T],
            :method,
            var,
        )
        |> x -> transform(
            x,
            [:benchmark, :ours] => ((x, y) -> x ./ y) => :ratio,
            [:n_customers, :xmax] => ((x, y) -> x ./ (2 * y)) => :customer_density,
        )
        |> x -> sort(
            x, 
            [:customer_density, :T],
        )
        |> x -> unstack(
            x, 
            :T,
            :customer_density,
            val == :percent_gap ? :ratio : val,
            combine = geomean,
        )
        |> x -> sort(
            x, 
            order(:T, rev = true),
        )
    )
    if val == :percent_gap
        result[:, Not(:T)] .= (1 .- (1 ./ result[:, Not(:T)])) .* 100.0
    end
    return result
end

time_taken_ratio_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :time_taken, :ratio
)
time_taken_percent_gap_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :time_taken, :percent_gap
)
time_taken_benchmark_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :time_taken, :benchmark
)
time_taken_ours_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :time_taken, :ours
)
time_taken_ratio_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :time_taken, :ratio
)
time_taken_percent_gap_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :time_taken, :percent_gap
)
time_taken_benchmark_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :time_taken, :benchmark
)
time_taken_ours_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :time_taken, :ours
)
time_taken_ratio_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :time_taken, :ratio
)
time_taken_percent_gap_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :time_taken, :percent_gap
)
time_taken_benchmark_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :time_taken, :benchmark
)
time_taken_ours_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :time_taken, :ours
)

LP_IP_gap_ratio_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :LP_IP_gap, :ratio
)
LP_IP_gap_benchmark_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :LP_IP_gap, :benchmark
)
LP_IP_gap_ours_none_df = make_pivottable(
    filtered_summary,
    false, false, false, 
    :LP_IP_gap, :ours
)
LP_IP_gap_ratio_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :LP_IP_gap, :ratio
)
LP_IP_gap_benchmark_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :LP_IP_gap, :benchmark
)
LP_IP_gap_ours_ngroute_df = make_pivottable(
    filtered_summary,
    false, true, true, 
    :LP_IP_gap, :ours
)
LP_IP_gap_ratio_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :LP_IP_gap, :ratio
)
LP_IP_gap_benchmark_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :LP_IP_gap, :benchmark
)
LP_IP_gap_ours_elementary_df = make_pivottable(
    filtered_summary,
    true, false, false,
    :LP_IP_gap, :ours
)

using ColorSchemes
using CairoMakie

begin 
    # masked time taken ratio

    density_range = collect(1.5:0.5:4.0)
    TB_range = collect(2.5:0.5:4.0)

    time_taken_none_ratio = Matrix(make_pivottable(
        filtered_summary,
        false, false, false, 
        :time_taken_filtered, :ratio
    )[end:-1:1,2:end])
    time_taken_elementary_ratio = Matrix(make_pivottable(
        filtered_summary,
        true, false, false, 
        :time_taken_filtered, :ratio
    )[end:-1:1,2:end])
    time_taken_ngroute_ratio = Matrix(make_pivottable(
        filtered_summary,
        false, true, true, 
        :time_taken_filtered, :ratio
    )[end:-1:1,2:end])

    datas = [
        time_taken_none_ratio,
        time_taken_ngroute_ratio,
        time_taken_elementary_ratio,
    ]
    titles = ["No elementarity", "ng-route", "Elementary"]

    cmax = extrema(vcat([vec(x)[findall(!isnan, vec(x))] for x in datas]...))[2]
    cmin = 1.0


    fig = CairoMakie.Figure(resolution = (500, 900), fontsize = 15)
    grid = fig[1,1] = GridLayout()
    fig[0,:] = Label(fig, "Ratio of time taken: our method against benchmark", font = :bold, fontsize = 17)
    for (ind, data) in enumerate(datas)
        if ind == length(datas)
            ax = Axis(
                grid[ind,1],
                xlabel = "Customer density",
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        else
            ax = Axis(
                grid[ind,1],
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        end
        density_range_unfolded = repeat(density_range, inner = length(TB_range))
        TB_range_unfolded = repeat(TB_range, outer = length(density_range))
        mask = findall(!isnan, vec(data))
        CairoMakie.heatmap!(
            ax,
            density_range_unfolded[mask],
            TB_range_unfolded[mask],
            vec(data)[mask],
            colorrange = (cmin, cmax),
        )
        for (density, TB, val) in zip(
            density_range_unfolded[mask],
            TB_range_unfolded[mask],
            vec(data)[mask],
        )
            textcolor = val < 6 ? :white : :black
            text!(
                ax, 
                string(round(val, digits = 2)),
                position = (density, TB),
                color = textcolor,
                align = (:center, :center),
            )
        end
    end
    Colorbar(grid[:,end+1], colorrange = (cmin, cmax))
    colgap!(grid, 20)
    rowgap!(grid, 10)
    display(fig)

    save("$(@__DIR__)/plots/time_taken_ratio_heatmap_masked.png", fig)
    save("$(@__DIR__)/plots/time_taken_ratio_heatmap_masked.pdf", fig)
end

begin 
    # time taken ratio

    density_range = collect(1.5:0.5:4.0)
    TB_range = collect(2.5:0.5:4.0)

    time_taken_none_ratio = Matrix(make_pivottable(
        filtered_summary,
        false, false, false, 
        :time_taken, :ratio
    )[end:-1:1,2:end])
    time_taken_elementary_ratio = Matrix(make_pivottable(
        filtered_summary,
        true, false, false, 
        :time_taken, :ratio
    )[end:-1:1,2:end])
    time_taken_ngroute_ratio = Matrix(make_pivottable(
        filtered_summary,
        false, true, true, 
        :time_taken, :ratio
    )[end:-1:1,2:end])

    datas = [
        time_taken_none_ratio,
        time_taken_ngroute_ratio,
        time_taken_elementary_ratio,
    ]
    titles = ["No elementarity", "ng-route", "Elementary"]

    cmin, cmax = extrema(vcat(datas...))


    fig = CairoMakie.Figure(resolution = (500, 900), fontsize = 15)
    grid = fig[1,1] = GridLayout()
    fig[0,:] = Label(fig, "Ratio of time taken: our method against benchmark", font = :bold, fontsize = 17)
    for (ind, data) in enumerate(datas)
        if ind == length(datas)
            ax = Axis(
                grid[ind,1],
                xlabel = "Customer density",
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        else
            ax = Axis(
                grid[ind,1],
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        end
        CairoMakie.heatmap!(
            ax,
            density_range,
            TB_range,
            data',
            colorrange = (cmin, cmax),
        )
        for i in 1:length(density_range), j in 1:length(TB_range)
            textcolor = data'[i,j] < 6 ? :white : :black
            text!(
                ax, 
                string(round(data'[i,j], digits = 2)),
                position = (density_range[i], TB_range[j]),
                color = textcolor,
                align = (:center, :center),
            )
        end
    end
    Colorbar(grid[:,end+1], colorrange = (cmin, cmax))
    colgap!(grid, 20)
    rowgap!(grid, 10)

    save("$(@__DIR__)/plots/time_taken_ratio_heatmap.pdf", fig)
    save("$(@__DIR__)/plots/time_taken_ratio_heatmap.png", fig)
end

begin 
    # masked time taken percent gap

    density_range = collect(1.5:0.5:4.0)
    TB_range = collect(2.5:0.5:4.0)

    time_taken_none_percent_gap = Matrix(make_pivottable(
        filtered_summary,
        false, false, false, 
        :time_taken_filtered, :percent_gap
    )[end:-1:1,2:end])
    time_taken_elementary_percent_gap = Matrix(make_pivottable(
        filtered_summary,
        true, false, false, 
        :time_taken_filtered, :percent_gap
    )[end:-1:1,2:end])
    time_taken_ngroute_percent_gap = Matrix(make_pivottable(
        filtered_summary,
        false, true, true, 
        :time_taken_filtered, :percent_gap
    )[end:-1:1,2:end])

    datas = [
        time_taken_none_percent_gap,
        time_taken_ngroute_percent_gap,
        time_taken_elementary_percent_gap,
    ]
    titles = ["No elementarity", "ng-route", "Elementary"]

    # cmax = extrema(vcat([vec(x)[findall(!isnan, vec(x))] for x in datas]...))[2]
    # cmin = extrema(vcat([vec(x)[findall(!isnan, vec(x))] for x in datas]...))[1]
    cmin, cmax = (0, 100)


    fig = CairoMakie.Figure(resolution = (600, 900), fontsize = 15)
    grid = fig[1,1] = GridLayout()
    fig[0,:] = Label(fig, "% improvement in time taken: our method against benchmark", font = :bold, fontsize = 17)
    for (ind, data) in enumerate(datas)
        if ind == length(datas)
            ax = Axis(
                grid[ind,1],
                xlabel = "Customer density",
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        else
            ax = Axis(
                grid[ind,1],
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        end
        density_range_unfolded = repeat(density_range, inner = length(TB_range))
        TB_range_unfolded = repeat(TB_range, outer = length(density_range))
        mask = findall(!isnan, vec(data))
        CairoMakie.heatmap!(
            ax,
            density_range_unfolded[mask],
            TB_range_unfolded[mask],
            vec(data)[mask],
            colorrange = (cmin, cmax),
        )
        for (density, TB, val) in zip(
            density_range_unfolded[mask],
            TB_range_unfolded[mask],
            vec(data)[mask],
        )
            textcolor = val < 50 ? :white : :black
            text!(
                ax, 
                string(round(val, digits = 2)),
                position = (density, TB),
                color = textcolor,
                align = (:center, :center),
            )
        end
    end
    Colorbar(grid[:,end+1], colorrange = (cmin, cmax))
    colgap!(grid, 20)
    rowgap!(grid, 10)
    display(fig)

    save("$(@__DIR__)/plots/time_taken_percent_gap_heatmap_masked.pdf", fig)
    save("$(@__DIR__)/plots/time_taken_percent_gap_heatmap_masked.png", fig)
end

begin 
    # time taken percent gap

    density_range = collect(1.5:0.5:4.0)
    TB_range = collect(2.5:0.5:4.0)

    time_taken_none_percent_gap = Matrix(make_pivottable(
        filtered_summary,
        false, false, false, 
        :time_taken, :percent_gap
    )[end:-1:1,2:end])
    time_taken_elementary_percent_gap = Matrix(make_pivottable(
        filtered_summary,
        true, false, false, 
        :time_taken, :percent_gap
    )[end:-1:1,2:end])
    time_taken_ngroute_percent_gap = Matrix(make_pivottable(
        filtered_summary,
        false, true, true, 
        :time_taken, :percent_gap
    )[end:-1:1,2:end])

    datas = [
        time_taken_none_percent_gap,
        time_taken_ngroute_percent_gap,
        time_taken_elementary_percent_gap,
    ]
    titles = ["No elementarity", "ng-route", "Elementary"]

    cmin, cmax = (0, 100)


    fig = CairoMakie.Figure(resolution = (600, 900), fontsize = 15)
    grid = fig[1,1] = GridLayout()
    fig[0,:] = Label(fig, "% improvement in time taken: our method against benchmark", font = :bold, fontsize = 17)
    for (ind, data) in enumerate(datas)
        if ind == length(datas)
            ax = Axis(
                grid[ind,1],
                xlabel = "Customer density",
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        else
            ax = Axis(
                grid[ind,1],
                xticks = density_range,
                ylabel = "T / B",
                yticks = TB_range,
                title = titles[ind],
            )
        end
        CairoMakie.heatmap!(
            ax,
            density_range,
            TB_range,
            data',
            colorrange = (cmin, cmax),
        )
        for i in 1:length(density_range), j in 1:length(TB_range)
            textcolor = data'[i,j] < 6 ? :white : :black
            text!(
                ax, 
                string(round(data'[i,j], digits = 2)),
                position = (density_range[i], TB_range[j]),
                color = textcolor,
                align = (:center, :center),
            )
        end
    end
    Colorbar(grid[:,end+1], colorrange = (cmin, cmax))
    colgap!(grid, 20)
    rowgap!(grid, 10)

    save("$(@__DIR__)/plots/time_taken_percent_gap_heatmap.pdf", fig)
    save("$(@__DIR__)/plots/time_taken_percent_gap_heatmap.png", fig)
end