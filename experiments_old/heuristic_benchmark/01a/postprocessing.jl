using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using CairoMakie

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)

data_fields = [
    :n_customers,
    :xmax,
    :T,
    :density,
]
method_fields = [
    :heuristic_use_SR3_cuts,
    :use_SR3_cuts,
]


args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args.index = collect(1:nrow(args))
results = innerjoin(
    results, args, on = setdiff(names(args), ["density", "index"])
)

(
    results 
    |> x -> sort!(
        x, 
        vcat(data_fields,  [:seed], method_fields),
    )
    |> x -> filter!(
        r -> !isinf(r.hc_IP_objective),
        x,
    )
    # |> x -> filter!(
    #     r -> r.hc_IP_objective > r.opt_IP_objective,
    #     x,
    # )
)

results.hc_IP_opt_IP_gap .= 100 .* results.hc_IP_objective ./ results.opt_IP_objective .- 100
results.opt_IP_opt_LP_gap .= 100 .* results.opt_IP_objective ./ results.opt_LP_objective .- 100
results.hc_IP_opt_LP_gap .= 100 .* results.hc_IP_objective ./ results.opt_LP_objective .- 100

results.hc_IP_opt_IP_ratio .= results.hc_IP_objective ./ results.opt_IP_objective
results.opt_IP_opt_LP_ratio .= results.opt_IP_objective ./ results.opt_LP_objective
results.hc_IP_opt_LP_ratio .= results.hc_IP_objective ./ results.opt_LP_objective

results.opt_IP_hc_IP_ratio .= results.opt_IP_objective ./ results.hc_IP_objective
results.opt_LP_opt_IP_ratio .= results.opt_LP_objective ./ results.opt_IP_objective
results.opt_LP_hc_IP_ratio .= results.opt_LP_objective ./ results.hc_IP_objective

summary = (
    results 
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields),
    )
    |> x -> combine(
        x, 
        nrow,
        # :h_time => mean,
        # :opt_time => mean,
        # :hc_IP_opt_IP_gap => mean,
        # :opt_IP_opt_LP_gap => mean,
        # :hc_IP_opt_LP_gap => mean,
        :hc_IP_opt_IP_ratio => geomean,
        :opt_IP_opt_LP_ratio => geomean,
        :hc_IP_opt_LP_ratio => geomean,
        :opt_IP_hc_IP_ratio => geomean,
        :opt_LP_opt_IP_ratio => geomean,
        :opt_LP_hc_IP_ratio => geomean,
    )
)

opt_IP_hc_IP_percent_gap_df = (
    summary
    |> x -> filter(
        r -> r.use_SR3_cuts,
        x,
    )
    |> x -> transform(
        x, 
        :opt_IP_hc_IP_ratio_geomean => (x -> (x .- 1) .* 100) => :opt_IP_hc_IP_ratio_geomean,
    )
    |> x -> unstack(
        x,
        :T,
        :n_customers,
        :opt_IP_hc_IP_ratio_geomean,
    )
)


opt_IP_opt_LP_percent_gap_df = (
    summary
    |> x -> filter(
        r -> r.use_SR3_cuts,
        x,
    )
    |> x -> transform(
        x, 
        :opt_IP_opt_LP_ratio_geomean => (x -> (x .- 1) .* 100) => :opt_IP_opt_LP_ratio_geomean,
    )
    |> x -> unstack(
        x,
        :T,
        :n_customers,
        :opt_IP_opt_LP_ratio_geomean,
    )
)


opt_IP_hc_IP_percent_gap_df

begin
    density_range = collect(2.5:0.5:5.0)
    TB_range = collect(4.0:0.5:5.0)
    data = Matrix(opt_IP_hc_IP_percent_gap_df[:, 2:end])

    (cmin, cmax) = round.(extrema(data), digits = 0)

    fig = CairoMakie.Figure(resolution = (600, 400), fontsize = 17)
    ax = Axis(
        fig[1,1],
        xlabel = "Customer density",
        xticks = density_range,
        ylabel = "Time horizon",
        yticks = TB_range,
        aspect = 1.6,
    )

    CairoMakie.heatmap!(
        ax,
        density_range, 
        TB_range,
        data',
        colorrange = (cmin, cmax),
        colormap = Makie.Reverse(:viridis),
    )
    for i in 1:length(density_range), j in 1:length(TB_range)
        val = data'[i,j]
        textcolor = val > (cmin + cmax) / 2 ? :white : :black
        text!(
            ax, 
            string(round(val, digits = 1)) * "%",
            position = (density_range[i], TB_range[j]),
            color = textcolor,
            align = (:center, :center),
            fontsize = 22,
        )
    end
    display(fig)

    save("$(@__DIR__)/plots/optimal_heuristic_percent_gap_heatmap.pdf", fig)
    save("$(@__DIR__)/plots/optimal_heuristic_percent_gap_heatmap.png", fig)
end

begin
    density_range = collect(2.5:0.5:5.0)
    TB_range = collect(4.0:0.5:5.0)
    data = Matrix(opt_IP_opt_LP_percent_gap_df[:, 2:end])

    (cmin, cmax) = round.(extrema(data), digits = 0)

    fig = CairoMakie.Figure(resolution = (600, 400), fontsize = 17)
    ax = Axis(
        fig[1,1],
        xlabel = "Customer density",
        xticks = density_range,
        ylabel = "Time horizon",
        yticks = TB_range,
        aspect = 1.6,
    )

    CairoMakie.heatmap!(
        ax,
        density_range, 
        TB_range,
        data',
        colorrange = (cmin, cmax),
        colormap = Makie.Reverse(:viridis),
    )
    for i in 1:length(density_range), j in 1:length(TB_range)
        val = data'[i,j]
        textcolor = val > (cmin + cmax) / 2 ? :white : :black
        text!(
            ax, 
            string(round(val, digits = 1)) * "%",
            position = (density_range[i], TB_range[j]),
            color = textcolor,
            align = (:center, :center),
            fontsize = 22,
        )
    end
    display(fig)

    save("$(@__DIR__)/plots/optimal_IP_LP_percent_gap_heatmap.pdf", fig)
    save("$(@__DIR__)/plots/optimal_IP_LP_percent_gap_heatmap.png", fig)
end