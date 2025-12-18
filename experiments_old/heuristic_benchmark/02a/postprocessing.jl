using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using CairoMakie

args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
undone_indexes = setdiff(
    1:nrow(args),
    [
        parse(Int, x[1:end-4])
        for x in readlines("$(@__DIR__)/out.txt")
    ]
)
d = 7
[
    intersect(i:d:nrow(args), undone_indexes)
    for i in 1:d
]

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)

data_fields = [
    :n_customers,
    :xmax,
    :T,
    :charge_cost_nlevels,
    :density,
]
method_fields = [
    :use_adaptive_ngroute,
    :use_SR3_cuts,
    :use_lmSR3_cuts,
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
        r -> (
            # !ismissing(r.time_taken_het) && 3600.0 > r.time_taken_het
            # && !ismissing(r.time_taken_hom) && 3600.0 > r.time_taken_hom
            r.LP_objective_last_SR3_het .< 1e10
            && r.LP_objective_last_SR3_hom .< 1e10
        ),
        x
    )
)



names(results)

results.IP_ratio_SR3_het_hom = results.IP_objective_last_SR3_het ./ results.IP_objective_last_SR3_hom
results.LP_ratio_SR3_het_hom = results.LP_objective_last_SR3_het ./ results.LP_objective_last_SR3_hom
sum(skipmissing(results.IP_ratio_SR3_het_hom) .< 1)
sum(skipmissing(results.IP_ratio_SR3_het_hom) .> 1)
sum(skipmissing(results.LP_ratio_SR3_het_hom) .< 1)
sum(skipmissing(results.LP_ratio_SR3_het_hom) .> 1)

results.IP_chargecost_ratio_SR3_het_hom = results.IP_chargecost_last_SR3_het ./ results.IP_chargecost_last_SR3_hom
results.LP_chargecost_ratio_SR3_het_hom = results.LP_chargecost_last_SR3_het ./ results.LP_chargecost_last_SR3_hom
sum(skipmissing(results.IP_chargecost_ratio_SR3_het_hom) .< 1)
sum(skipmissing(results.IP_chargecost_ratio_SR3_het_hom) .> 1)
sum(skipmissing(results.LP_chargecost_ratio_SR3_het_hom) .< 1)
sum(skipmissing(results.LP_chargecost_ratio_SR3_het_hom) .> 1)

results.IP_LP_ratio_het = results.IP_objective_last_SR3_het ./ results.LP_objective_last_SR3_het
results.IP_LP_ratio_hom = results.IP_objective_last_SR3_hom ./ results.LP_objective_last_SR3_hom


summary = (
    results 
    |> x -> groupby(
        x, 
        data_fields
    )
    |> x -> combine(
        x, 
        nrow,
        :IP_ratio_SR3_het_hom => geomean,
        :LP_ratio_SR3_het_hom => geomean,
        :IP_chargecost_ratio_SR3_het_hom => geomean,
        :LP_chargecost_ratio_SR3_het_hom => geomean,
        :IP_LP_ratio_het => geomean,
        :IP_LP_ratio_hom => geomean,
    )
    |> x -> transform(
        x,
        :IP_ratio_SR3_het_hom_geomean => (x -> (x .- 1) .* 100) => :IP_gap_SR3_het_hom_geomean,
        :LP_ratio_SR3_het_hom_geomean => (x -> (x .- 1) .* 100) => :LP_gap_SR3_het_hom_geomean,
        :IP_chargecost_ratio_SR3_het_hom_geomean => (x -> (x .- 1) .* 100) => :IP_chargecost_gap_SR3_het_hom_geomean,
        :LP_chargecost_ratio_SR3_het_hom_geomean => (x -> (x .- 1) .* 100) => :LP_chargecost_gap_SR3_het_hom_geomean,
        :IP_LP_ratio_het_geomean => (x -> (x .- 1) .* 100) => :IP_LP_gap_het_geomean,
        :IP_LP_ratio_hom_geomean => (x -> (x .- 1) .* 100) => :IP_LP_gap_hom_geomean,
    )
)

begin
    T_range = sort(unique(results.T))
    T_range = T_range[1:1]

    for name in [
        :IP_gap_SR3_het_hom_geomean,
        :IP_chargecost_gap_SR3_het_hom_geomean,
        :IP_LP_gap_het_geomean,
        :IP_LP_gap_hom_geomean,
    ],
        T in T_range

    # begin
        density_range = collect(2.5:0.5:5.0)
        charge_cost_nlevels_range = collect(2:1:5)
        df = (
            summary 
            |> x -> filter(
                r -> r.T == T,
                x
            )
            |> x -> unstack(
                x, 
                :charge_cost_nlevels,
                :n_customers,
                name,
            )
        )
        data = Matrix(df[:, 2:end])

        # (cmin, cmax) = (-25.0, -5.0)
        (cmin, cmax) = (floor(minimum(data), digits = 0), ceil(maximum(data), digits = 0))

        fig = CairoMakie.Figure(resolution = (600, 400), fontsize = 17)
        ax = Axis(
            fig[1,1],
            xlabel = "Customer density",
            xticks = density_range,
            ylabel = "Number of charge cost levels",
            yticks = charge_cost_nlevels_range,
            aspect = 1.5,
        )

        CairoMakie.heatmap!(
            ax,
            density_range, 
            charge_cost_nlevels_range,
            data',
            colorrange = (cmin, cmax),
            colormap = Makie.Reverse(:viridis),
        )
        for i in 1:length(density_range), j in 1:length(charge_cost_nlevels_range)
            val = data'[i,j]
            textcolor = val > (cmin + cmax) / 2 ? :white : :black
            text!(
                ax, 
                string(round(val, digits = 1)) * "%",
                position = (density_range[i], charge_cost_nlevels_range[j]),
                color = textcolor,
                align = (:center, :center),
                fontsize = 22,
            )
        end
        resize_to_layout!(fig)
        display(fig)

        save("$(@__DIR__)/plots/$(string(name))_$(T)_heatmap.pdf", fig)
        save("$(@__DIR__)/plots/$(string(name))_$(T)_heatmap.png", fig)
    end
end