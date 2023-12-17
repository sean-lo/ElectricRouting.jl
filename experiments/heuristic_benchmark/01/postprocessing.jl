using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using CairoMakie

results = vcat(
    [
        CSV.read(filepath, DataFrame)
        for filepath in glob("combined_*.csv", "$(@__DIR__)")
    ]...
)
names(results)

data_fields = [
    :n_customers,
    :xmax,
    :T,
    # :density,
]
method_fields = [
    :heuristic_use_adaptive_ngroute,
    :heuristic_use_SR3_cuts,
    :heuristic_use_lmSR3_cuts,
    :use_adaptive_ngroute,
    :use_SR3_cuts,
    :use_lmSR3_cuts,
]

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
)

# results.LP_IP_gap_first .*= 100
# results.LP_IP_gap_last .*= 100

results.hc_opt_IP_gap .= 100 .* results.hc_IP_objective ./ results.opt_IP_objective .- 100
results.opt_IP_LP_gap .= 100 .* results.opt_IP_objective ./ results.opt_LP_objective .- 100

results.hc_opt_IP_ratio .= results.hc_IP_objective ./ results.opt_IP_objective
results.opt_IP_LP_ratio .= results.opt_IP_objective ./ results.opt_LP_objective

results.opt_hc_IP_ratio .= results.opt_IP_objective ./ results.hc_IP_objective
results.opt_LP_IP_ratio .= results.opt_LP_objective ./ results.opt_IP_objective


summary = (
    results 
    |> x -> groupby(
        x, 
        vcat(data_fields, :use_SR3_cuts)
    )
    |> x -> combine(
        x, 
        nrow,
        :h_time => mean,
        :opt_time => mean,
        :h_IP_objective => mean,
        :hc_IP_objective => mean,
        :opt_IP_objective => mean,
        :opt_LP_objective => mean,
        :hc_opt_IP_gap => mean,
        :opt_IP_LP_gap => mean,
        :hc_opt_IP_ratio => geomean,
        :opt_IP_LP_ratio => geomean,
        :opt_hc_IP_ratio => geomean,
        :opt_LP_IP_ratio => geomean,
    )
)

opt_hc_IP_percent_gap_df = (
    summary
    |> x -> filter(
        r -> (
            r.use_SR3_cuts
        ), 
        x
    )
    |> x -> transform(
        x, 
        :opt_hc_IP_ratio_geomean => (x -> (x .- 1) .* 100) => :opt_hc_IP_ratio_geomean,
    )
    |> x -> unstack(
        x,
        :T,
        :n_customers,
        :opt_hc_IP_ratio_geomean,
    )
)

begin
    density_range = collect(2.5:0.5:5.0)
    TB_range = collect(4.0:0.5:5.0)
    data = Matrix(opt_hc_IP_percent_gap_df[end:-1:1, 2:end])

    (cmin, cmax) = round.(extrema(data), digits = 0)

    fig = CairoMakie.Figure(resolution = (600, 400), fontsize = 17)
    mygrid = fig[1,1] = GridLayout()
    fig[0,1] = Label(
        fig, 
        "% reduction in cost: our method against business-as-usual", 
        font = :bold, fontsize = 18
    )

    ax = Axis(
        mygrid[1,1],
        xlabel = "Customer density",
        xticks = density_range,
        ylabel = "Time horizon",
        yticks = TB_range,
        # title = titles[ind],
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
        textcolor = val > -15 ? :white : :black
        text!(
            ax, 
            string(Int(round(val, digits = 0))) * "%",
            position = (density_range[i], TB_range[j]),
            color = textcolor,
            align = (:center, :center),
            fontsize = 22,
        )
    end
    Colorbar(mygrid[:,end+1], colorrange = (cmin, cmax), colormap = Makie.Reverse(:viridis))
    colgap!(mygrid, 20)
    rowgap!(mygrid, 10)
    display(fig)

    save("$(@__DIR__)/plots/optimal_heuristic_percent_gap_heatmap.pdf", fig)
    save("$(@__DIR__)/plots/optimal_heuristic_percent_gap_heatmap.png", fig)
end