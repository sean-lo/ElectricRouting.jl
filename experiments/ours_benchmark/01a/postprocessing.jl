using DataFrames, CSV
using DelimitedFiles
using StatsBase
using ColorSchemes
using Plots

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
    :elementary,
    :ngroute,
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

(
    results
    |> x -> filter!(
        r -> (
            # Note: LP_IP_gap only ends up being negative for certain large instances,
            # where the time limit is reached and one still has artificial paths in the LP solution
            r.LP_IP_gap ≥ 0 
            || !r.time_limit_reached
        ),
        x, 
    )
)

results.LP_IP_gap .*= 100
results.LP_IP_gap .= ifelse.(results.LP_IP_gap .≤ 1e-10, 0.0, results.LP_IP_gap)

summary_df = (
    results 
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x, 
        nrow,
        :time_taken => geomean => :time_taken,
        :sp_time_taken_mean => geomean => :sp_time_taken_mean,
        :counter => mean => :counter,
        :converged => mean => :converged,
        :time_limit_reached => mean => :time_limit_reached,
        :LP_objective => mean,
        :IP_objective => mean,
        :LP_IP_gap => mean => :LP_IP_gap,
        :mean_subpath_length => mean => :mean_subpath_length,
        :mean_path_length => mean => :mean_path_length,
        :mean_ps_length => mean => :mean_ps_length,
    )
)

filtered_summary_df = summary_df |> 
    x -> filter(
        r -> (
            !(r.xmax == 2.0 && r.T == 36000)
            # && !(r.density in [3.5, 4.0] && r.T == 72000)
        ),
        x
    )


filtered_summary_df.time_taken_filtered = copy(filtered_summary_df.time_taken)
filtered_summary_df.time_taken_filtered[filtered_summary_df.converged .!= 1.0] .= NaN
filtered_summary_df.method_combined = [
    join(filtered_summary_df[i, method_fields], "_")
    for i in 1:nrow(filtered_summary_df)
]

(
    filtered_summary_df
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
        :benchmark_false_false, :ours_false_false,
        [:benchmark_false_false, :ours_false_false] => ((x, y) -> (1.0 .- y ./ x) .* 100) => :gap_false_false_false,
        :benchmark_false_true, :ours_false_true,
        [:benchmark_false_true, :ours_false_true] => ((x, y) -> (1.0 .- y ./ x) .* 100) => :gap_false_true_true,
        :benchmark_true_false, :ours_true_false,
        [:benchmark_true_false, :ours_true_false] => ((x, y) -> (1.0 .- y ./ x) .* 100) => :gap_true_false_false,
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


function make_pivottable(
    df::DataFrame,
    elementary::Bool,
    ngroute::Bool,
    var::Symbol, # one of :time_taken, :LP_IP_gap_mean
    val::Symbol, # one of :benchmark, :ours, :ratio, :percent_gap
)
    if !((elementary, ngroute) in [
        (true, false),
        (false, false), 
        (false, true),
    ])
        error()
    end
    result = (
        df 
        |> x -> filter(
            r -> (
                r.elementary == elementary
                && r.ngroute == ngroute
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
    filtered_summary_df,
    false, false,
    :time_taken, :ratio
)
time_taken_percent_gap_none_df = make_pivottable(
    filtered_summary_df,
    false, false,
    :time_taken, :percent_gap
)
time_taken_benchmark_none_df = make_pivottable(
    filtered_summary_df,
    false, false,
    :time_taken, :benchmark
)
time_taken_ours_none_df = make_pivottable(
    filtered_summary_df,
    false, false,
    :time_taken, :ours
)
time_taken_ratio_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :time_taken, :ratio
)
time_taken_percent_gap_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :time_taken, :percent_gap
)
time_taken_benchmark_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :time_taken, :benchmark
)
time_taken_ours_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :time_taken, :ours
)
time_taken_ratio_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :time_taken, :ratio
)
time_taken_percent_gap_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :time_taken, :percent_gap
)
time_taken_benchmark_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :time_taken, :benchmark
)
time_taken_ours_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :time_taken, :ours
)

LP_IP_gap_ratio_none_df = make_pivottable(
    filtered_summary_df,
    false, false,
    :LP_IP_gap, :ratio
)
LP_IP_gap_benchmark_none_df = make_pivottable(
    filtered_summary_df,
    false, false,
    :LP_IP_gap, :benchmark
)
LP_IP_gap_ours_none_df = make_pivottable(
    filtered_summary_df,
    false, false,
    :LP_IP_gap, :ours
)
LP_IP_gap_ratio_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :LP_IP_gap, :ratio
)
LP_IP_gap_benchmark_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :LP_IP_gap, :benchmark
)
LP_IP_gap_ours_ngroute_df = make_pivottable(
    filtered_summary_df,
    false, true,
    :LP_IP_gap, :ours
)
LP_IP_gap_ratio_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :LP_IP_gap, :ratio
)
LP_IP_gap_benchmark_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :LP_IP_gap, :benchmark
)
LP_IP_gap_ours_elementary_df = make_pivottable(
    filtered_summary_df,
    true, false,
    :LP_IP_gap, :ours
)

using ColorSchemes
using CairoMakie

begin
    density_range = collect(1.5:0.5:4.0)
    TB_range = collect(2.5:0.5:4.0)
    for (elementary, ngroute, name) in [
        (false, false, "none"),
        (false, true, "ngroute"),
        (true, false, "elementary"),
    ]
        data = -Matrix(make_pivottable(
            filtered_summary_df,
            elementary, ngroute,
            :time_taken_filtered, :percent_gap
        )[end:-1:1,2:end])


        fig = CairoMakie.Figure(resolution = (600, 400), fontsize = 18)
        mygrid = fig[1,1] = GridLayout()
        ax = Axis(
            mygrid[1,1],
            xlabel = "Task density",
            xticks = density_range,
            ylabel = "Scaled time horizon",
            yticks = TB_range,
        )
        hidedecorations!(ax, label = false, ticklabels = false, ticks = false, minorticks = false)
        cmin, cmax = (-100, 0)
        density_range_unfolded = repeat(density_range, inner = length(TB_range))
        TB_range_unfolded = repeat(TB_range, outer = length(density_range))
        mask = findall(!isnan, vec(data))
        CairoMakie.heatmap!(
            ax,
            density_range_unfolded[mask],
            TB_range_unfolded[mask],
            vec(data)[mask],
            colorrange = (cmin, cmax),
            colormap = Makie.Reverse(:viridis),
            nan_color = :grey30,
        )
        for (density, TB, val) in zip(
            density_range_unfolded[mask],
            TB_range_unfolded[mask],
            vec(data)[mask],
        )
            textcolor = val > -50 ? :white : :black
            text!(
                ax, 
                string(round(val, digits = 1)) * "%",
                fontsize = 22,
                position = (density, TB),
                color = textcolor,
                align = (:center, :center),
            )
        end
        display(fig)
        save("$(@__DIR__)/plots/time_taken_percent_gap_heatmap_masked_$name.pdf", fig)
        save("$(@__DIR__)/plots/time_taken_percent_gap_heatmap_masked_$name.png", fig)
    end
end