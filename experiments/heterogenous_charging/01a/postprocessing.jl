using DataFrames, CSV
using StatsBase
using StatsPlots

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
results.density .= results.n_customers ./ (results.xmax .* 2)
data_fields = [
    :n_customers,
    :xmax,
    :T,
    :density,
    :charge_cost_nlevels,
]
method_fields = [
    :method,
    :ngroute_neighborhood_charging_size,
    :use_lmSR3_cuts,
    :sparse_graph,
]
(
    results 
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

results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100
results.LP_IP_gap_last_SR3 .*= 100
results.LP_IP_gap_first .= ifelse.(results.LP_IP_gap_first .≤ 1e-10, 0.0, results.LP_IP_gap_first)
results.LP_IP_gap_last .= ifelse.(results.LP_IP_gap_last .≤ 1e-10, 0.0, results.LP_IP_gap_last)
results.LP_IP_gap_last_SR3 .= ifelse.(results.LP_IP_gap_last_SR3 .≤ 1e-10, 0.0, results.LP_IP_gap_last_SR3)

all_results = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args.index = 1:nrow(args)
all_results.CG_LP_IP_gap .*= 100
all_results.CG_LP_IP_gap .= ifelse.(all_results.CG_LP_IP_gap .≤ 1e-10, 0.0, all_results.CG_LP_IP_gap)

all_results = outerjoin(all_results, args, on = :instance => :index, makeunique = true)
select!(all_results, Not(:method_1))

joincols = vcat(data_fields, method_fields, [:seed])
results = innerjoin(
    results, 
    select(args, vcat(joincols, [:index])), 
    on = joincols,
)
remove_inds = results[results.LP_objective_last_SR3 .> 1e10, :index]

filter!(
    r -> !(r.index in remove_inds),
    results
)
filter!(
    r -> !(r.instance in remove_inds),
    all_results
)

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
        :time_taken_total_SR3 => geomean => :time_taken_total_SR3,
        :sp_time_taken_mean_first => geomean => :sp_time_taken_mean_first,
        :sp_time_taken_mean_last => geomean => :sp_time_taken_mean_last,
        :sp_time_taken_mean_last_SR3 => geomean => :sp_time_taken_mean_last_SR3,
        :LP_objective_first => mean => :LP_objective_first,
        :LP_objective_last => mean => :LP_objective_last,
        :LP_objective_last_SR3 => mean => :LP_objective_last_SR3,
        :IP_objective_first => mean => :IP_objective_first,
        :IP_objective_last => mean => :IP_objective_last,
        :IP_objective_last_SR3 => mean => :IP_objective_last_SR3,
        :LP_IP_gap_first => mean => :LP_IP_gap_first,
        :LP_IP_gap_last => mean => :LP_IP_gap_last,
        :LP_IP_gap_last_SR3 => mean => :LP_IP_gap_last_SR3,
        :neighborhood_size_mean_first => mean => :neighborhood_size_mean_first,
        :neighborhood_size_mean_last => mean => :neighborhood_size_mean_last,
        :neighborhood_size_mean_last_SR3 => mean => :neighborhood_size_mean_last_SR3,
        :implemented_SR3_cuts_count_total => mean => :implemented_SR3_cuts_count,
    )
)
data_df = (
    summary 
    |> x -> filter(
        r -> (
            r.charge_cost_nlevels == 1
            && r.sparse_graph == true
        ),
        x
    )
)

# Plotting: aggregate grouped bar plots with StatsPlots
# 1 plot per charge_cost_nlevels
begin
    n_customers_range = sort(unique(results.n_customers))
    charge_cost_nlevels_range = 1:5
    groupnames = ["CG", "CG + ng-relaxation", "CG + ng-relaxation + cuts"]
    sizenames = string.(n_customers_range)
    for charge_cost_nlevels in charge_cost_nlevels_range, 
        sparse_graph in [false, true]
        data_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.charge_cost_nlevels == charge_cost_nlevels
                    && r.sparse_graph == sparse_graph
                ),
                x
            )
        )
        time_taken_df = (
            data_df
            |> x -> select(x, :n_customers, :time_taken_first, :time_taken_total, :time_taken_total_SR3)
        )
        LP_IP_gap_df = (
            data_df
            |> x -> select(x, :n_customers, :LP_IP_gap_first, :LP_IP_gap_last, :LP_IP_gap_last_SR3)
        )
        p1_ytickvalues = [1, 10, 100, 1000, 3600]
        p1 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(time_taken_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.8,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Computational time (s)",
            ylim = 10.0.^(-0.5, 4.0), 
            yticks = (p1_ytickvalues, string.(p1_ytickvalues)),
            grid = :y,
            legend = :bottomright,
        )
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(charge_cost_nlevels)_$(sparse_graph).pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(charge_cost_nlevels)_$(sparse_graph).png")
        p2_ytickvalues = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        p2 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(LP_IP_gap_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.8,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Optimality gap (%)",
            ylim = 10.0.^(-1, 2), 
            yticks = (p2_ytickvalues, string.(p2_ytickvalues)),
            grid = :y,
            legend = :bottomright,
        )
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(charge_cost_nlevels)_$(sparse_graph).pdf")
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(charge_cost_nlevels)_$(sparse_graph).png")
    end
end


# Plotting: aggregate grouped bar plots with StatsPlots
# 1 plot per n_customers
begin
    n_customers_range = sort(unique(results.n_customers))
    charge_cost_nlevels_range = 1:5
    groupnames = ["CG", "CG + ng-relaxation", "CG + ng-relaxation + cuts"]
    sizenames = string.(charge_cost_nlevels_range)
    for n_customers in n_customers_range,
        sparse_graph in [false, true]
        data_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.n_customers == n_customers
                    && r.sparse_graph == sparse_graph
                ),
                x
            )
        )
        time_taken_df = (
            data_df
            |> x -> select(x, :n_customers, :time_taken_first, :time_taken_total, :time_taken_total_SR3)
        )
        LP_IP_gap_df = (
            data_df
            |> x -> select(x, :n_customers, :LP_IP_gap_first, :LP_IP_gap_last, :LP_IP_gap_last_SR3)
        )
        p1_ytickvalues = [1, 10, 100, 1000, 3600]
        p1 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(time_taken_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.8,
            framestyle = :box,
            xlabel = "Number of charge cost levels",
            yscale = :log10,
            ylabel = "Computational time (s)",
            ylim = 10.0.^(-0.5, 4.0), 
            yticks = (p1_ytickvalues, string.(p1_ytickvalues)),
            grid = :y,
            legend = :bottomright,
        )
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(n_customers)_$(sparse_graph).pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(n_customers)_$(sparse_graph).png")
        p2_ytickvalues = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        p2 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(LP_IP_gap_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.8,
            framestyle = :box,
            xlabel = "Number of charge cost levels",
            yscale = :log10,
            ylabel = "Optimality gap (%)",
            ylim = 10.0.^(-1, 2),
            yticks = (p2_ytickvalues, string.(p2_ytickvalues)),
            grid = :y,
            legend = :bottomright,
        )
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(n_customers)_$(sparse_graph).pdf")
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(n_customers)_$(sparse_graph).png")
    end
end

using ColorSchemes
using CairoMakie

begin
    n_customers_range = sort(unique(results.n_customers))
    charge_cost_nlevels_range = sort(unique(results.charge_cost_nlevels))
    n_customers_range_unfolded = repeat(n_customers_range, inner = length(charge_cost_nlevels_range))
    charge_cost_nlevels_range_unfolded = repeat(charge_cost_nlevels_range, outer = length(n_customers_range))
    for sparse_graph in [false, true]
        data = Matrix((
            summary
            |> filter(
                r -> r.sparse_graph == sparse_graph,
            )
            |> x -> unstack(
                x, 
                :charge_cost_nlevels,
                :n_customers,
                :time_taken_total_SR3
            )
            # |> x -> select(x, :charge_cost_nlevels, :time_taken_first, :time_taken_total, :time_taken_total_SR3)
        )[:,2:end])

        fig = CairoMakie.Figure(resolution = (600, 400), fontsize = 18)
        mygrid = fig[1,1] = GridLayout()
        ax = Axis(
            mygrid[1,1],
            xlabel = "# of tasks",
            xticks = n_customers_range,
            ylabel = "# of charge cost levels",
            yticks = charge_cost_nlevels_range,
        )
        CairoMakie.heatmap!(
            ax,
            n_customers_range_unfolded,
            charge_cost_nlevels_range_unfolded,
            vec(data),
            # colorrange = (cmin, cmax),
            colormap = Makie.Reverse(:viridis),
            colorscale = log10,
            # nan_color = :grey30,
        )
        for (n_customers, charge_cost_level, val) in zip(
            n_customers_range_unfolded,
            charge_cost_nlevels_range_unfolded,
            vec(data),
        )
            # textcolor = val > -50 ? :white : :black
            textcolor = val < 250 ? :black : :white
            text!(
                ax, 
                # string(round(val, sigdigits = 3)),
                string(round(val, digits = 1)),
                fontsize = 22,
                position = (n_customers, charge_cost_level),
                color = textcolor,
                align = (:center, :center),
            )
        end
        display(fig)
        save("$(@__DIR__)/plots/heatmap_time_taken_n_customers_charge_cost_levels_$(sparse_graph).pdf", fig)
        save("$(@__DIR__)/plots/heatmap_time_taken_n_customers_charge_cost_levels_$(sparse_graph).png", fig)
    end
end