using DataFrames, CSV
using DelimitedFiles
using StatsBase
using StatsPlots
using ColorSchemes

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
    :ngroute,
    :use_lmSR3_cuts,
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
    |> x -> unique!(
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
        if !isfile("$(@__DIR__)/results/$ind.csv")
            continue
        end
        data = CSV.read("$(@__DIR__)/results/$ind.csv", DataFrame)
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

all_data = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
all_data.CG_LP_IP_gap *= 100
sort!(all_data, [:instance, :iteration])
(
    all_data 
    |> x -> groupby(
        x, 
        :instance,
    )
    |> x -> transform!(
        x, 
        :CG_time_taken => cumsum,
    )
    |> x -> select!(
        x, 
        [:instance, :iteration, :CG_time_taken, :CG_time_taken_cumsum],
        Not([:instance, :iteration, :CG_time_taken, :CG_time_taken_cumsum]),
    )
)
args_df = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args_df.index = 1:nrow(args_df)
leftjoin!(
    all_data, 
    args_df,
    on = :instance => :index,
    makeunique = true,
)

# Trajectory plot
begin
    p = Plots.plot(
        xscale = :log10, 
        xlabel = "Computational time (s)",
        yaxis = (:log10, [0.1, 110.0]),
        ylabel = "Optimality gap (%)",
    )
    ytickvalues = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
    yticks!(p, ytickvalues, string.(ytickvalues))
    xtickvalues = [1, 10, 100, 1000, 3600]
    xticks!(p, xtickvalues, string.(xtickvalues))
    gb = (
        all_data
        |> x -> filter(
            r -> (
                r.xmax == 4.0
                && r.density == 4.0
                && r.method_1 == "benchmark"
            ),
            all_data
        )
        |> x -> groupby(
            x, :instance,
        )
        # |> x -> length(unique(x.instance))
    )
    for (ind, g) in enumerate(gb)
        g.CG_LP_IP_gap[g.CG_LP_IP_gap .< 0.1] .= 0.05
        plot!(
            vec(g.CG_time_taken_cumsum), vec(g.CG_LP_IP_gap), 
            markersize = 2, 
            # shape = :circle, 
            color = :grey, 
            alpha = 0.7,
            label = false,
        )
        plot!(
            [g.CG_time_taken_cumsum[end]], [g.CG_LP_IP_gap[end]], 
            shape = :circle, 
            color = :grey,
            label = (ind == 1 ? "Label-setting" : false),
        )
    end
    gb = (
        all_data
        |> x -> filter(
            r -> (
                r.xmax == 4.0
                && r.density == 4.0
                && r.method_1 == "ours"
            ),
            all_data
        )
        |> x -> groupby(
            x, :instance,
        )
        # |> x -> length(unique(x.instance))
    )
    for (ind, g) in enumerate(gb)
        g.CG_LP_IP_gap[g.CG_LP_IP_gap .< 0.1] .= 0.05
        plot!(
            vec(g.CG_time_taken_cumsum), vec(g.CG_LP_IP_gap), 
            markersize = 2, 
            # shape = :circle, 
            color = :green, 
            alpha = 0.7,
            label = false,
        )
        plot!(
            [g.CG_time_taken_cumsum[end]], [g.CG_LP_IP_gap[end]], 
            shape = :circle, 
            color = :green,
            label = (ind == 1 ? "Two-level label-setting" : false),
        )
    end
    savefig(p, "$(@__DIR__)/plots/LP_IP_gap_time_taken_trajectories.png")
    savefig(p, "$(@__DIR__)/plots/LP_IP_gap_time_taken_trajectories.pdf")
end

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
        :mean_subpath_length => mean => :mean_subpath_length,
        :mean_path_length => mean => :mean_path_length,
        :mean_ps_length => mean => :mean_ps_length,
    ) |>
    x -> filter!(
        r -> (
            !(r.xmax == 2.0 && r.T == 36000)
        ),
        x
    )
)

# Bar plots comparing benchmark to ours
begin
groupnames = ["Label-setting", "Two-level label-setting"]
xmax_T_range = collect(zip(unique(summary[!, [:xmax, :T]]).xmax, unique(summary[!, [:xmax, :T]]).T))
    for (xmax, T) in xmax_T_range
        data_df = filter(
            r -> r.xmax == xmax && r.T == T,
            summary
        )
        sizes = sort(unique(data_df.n_customers))
        sizenames = string.(sizes, pad = 2)
        time_taken_df = hcat(
            data_df
            |> x -> filter(r -> r.method == "benchmark", x) 
            |> x -> select(x, :T, :n_customers, :time_taken_total_SR3 => :time_taken_total_SR3_benchmark),
            data_df
            |> x -> filter(r -> r.method == "ours", x) 
            |> x -> select(x, :time_taken_total_SR3 => :time_taken_total_SR3_ours)
        )
        
        LP_IP_gap_df = hcat(
            data_df
            |> x -> filter(r -> r.method == "benchmark", x) 
            |> x -> select(x, :n_customers, :LP_IP_gap_last_SR3 => :LP_IP_gap_last_SR3_benchmark),
            data_df
            |> x -> filter(r -> r.method == "ours", x) 
            |> x -> select(x, :LP_IP_gap_last_SR3 => :LP_IP_gap_last_SR3_ours)
        )

        p1_ytickvalues = [1, 10, 100, 1000, 3600]
        p1 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(time_taken_df[!, end-1:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.75,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Computational time (s)",
            ylim = 10.0.^(-0.5, 4.0), 
            grid = :y,
            legend = :topleft,
        )
        yticks!(p1, p1_ytickvalues, string.(p1_ytickvalues))
        xticks!(p1, 0.5:1:(length(sizes)-0.5), string.(sizes))
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_SRI_$(xmax)_$(T)_lmSRI.pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_SRI_$(xmax)_$(T)_lmSRI.png")
    end
end
