using DataFrames, CSV
using DelimitedFiles
using StatsBase
using Plots
using StatsPlots
using ColorSchemes

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
results.density .= results.n_customers ./ (results.xmax .* 2)
data_fields = [
    :n_customers,
    :xmax,
    :T,
]
method_fields = [
    :use_lmSR3_cuts,
    :sparse_graph,
]
(
    results 
    |> x -> sort!(
        x, 
        vcat(data_fields,  [:seed], method_fields),
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

# results matching
all_data = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
args = CSV.read("$(@__DIR__)/args.csv", DataFrame)

for i in 1:nrow(args)
    if Bool(args[i, :use_lmSR3_cuts])
        continue
    end
    density = args[i, :density]
    n_customers = args[i, :n_customers]
    n_charging = args[i, :n_charging]
    T = args[i, :T]
    use_lmSR3_cuts = Bool(args[i, :use_lmSR3_cuts])
    seed = args[i, :seed]
    # println("$density, $n_customers, $n_charging, $T, $ngroute_neighborhood_charging_size, $use_lmSR3_cuts")
    cond = (
        (results.seed .== seed)
        .& (results.n_customers .== n_customers)
        .& (results.n_charging .== n_charging)
        .& (results.T .== T)
        .& (results.use_lmSR3_cuts .== use_lmSR3_cuts)
    )
    ind_df = filter(r -> r.instance == i, all_data)
    ind = findfirst(x -> x in ["use_SR3_cuts", "use_lmSR3_cuts"], filter(r -> r.instance == i, all_data).method)
    if isnothing(ind)
        ind = nrow(ind_df)
    end
    sum(ind_df[1:ind, :CG_time_taken])

    results[cond, :time_taken_total] .= sum(ind_df[1:ind, :CG_time_taken])
    results[cond, :sp_time_taken_mean_last] .= ind_df[ind, :CG_sp_time_taken_mean]
    results[cond, :LP_objective_last] .= ind_df[ind, :CGLP_objective]
    results[cond, :IP_objective_last] .= ind_df[ind, :CGIP_objective]
    results[cond, :LP_IP_gap_last] .= ind_df[ind, :CG_LP_IP_gap]
    results[cond, :neighborhood_size_mean_last] .= ind_df[ind, :ngroute_neighborhood_size]
    println("$i, $(findall(cond))")
end

results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100
results.LP_IP_gap_last_SR3 .*= 100
results.LP_IP_gap_first .= ifelse.(results.LP_IP_gap_first .≤ 1e-10, 0.0, results.LP_IP_gap_first)
results.LP_IP_gap_last .= ifelse.(results.LP_IP_gap_last .≤ 1e-10, 0.0, results.LP_IP_gap_last)
results.LP_IP_gap_last_SR3 .= ifelse.(results.LP_IP_gap_last_SR3 .≤ 1e-10, 0.0, results.LP_IP_gap_last_SR3)

summary_df = (
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
        :time_taken_first => mean => :time_taken_first,
        :time_taken_total => mean => :time_taken_total,
        :time_taken_total_SR3 => mean => :time_taken_total_SR3,
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

(
    summary_df
    |> x -> filter!(
        r -> !(r.n_customers in [70,80]),
        x
    )
)
# Plotting: aggregate grouped bar plots with StatsPlots
# For presentation

begin
    groupnames = ["CG", "CG + ng-relaxation", "CG + ng-relaxation + cuts"]
    sizenames = string.(unique(summary_df.n_customers))
    for sparse_graph in [false, true]
        time_taken_df = (
            summary_df
            |> x -> filter(r -> r.sparse_graph, x)
            |> x -> select(x, :n_customers, :time_taken_first, :time_taken_total, :time_taken_total_SR3)
        )
        LP_IP_gap_df = (
            summary_df
            |> x -> filter(r -> r.sparse_graph, x)
            |> x -> select(x, :n_customers, :LP_IP_gap_first, :LP_IP_gap_last, :LP_IP_gap_last_SR3)
        )
        p1_ytickvalues = [1, 10, 100, 1000, 3600]
        p1 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(time_taken_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.75,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Computational time (s)",
            ylim = 10.0.^(-0.5, 4.0), 
            grid = :y,
            legend = :bottomright,
        )
        Plots.yticks!(p1, p1_ytickvalues, string.(p1_ytickvalues))
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(sparse_graph).pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(sparse_graph).png")
        p2_ytickvalues = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        p2 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(LP_IP_gap_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.75,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Optimality gap (%)",
            ylim = 10.0.^(-1.0, 2),
            grid = :y,
            legend = :bottomright,
        )
        Plots.yticks!(p2, p2_ytickvalues, string.(p2_ytickvalues))
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(sparse_graph).pdf")
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(sparse_graph).png")
    end
end
