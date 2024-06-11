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
    :ngroute_neighborhood_charging_size,
    :use_lmSR3_cuts,
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
    ngroute_neighborhood_charging_size = String(args[i, :ngroute_neighborhood_charging_size])
    # use_lmSR3_cuts = Bool(args[i, :use_lmSR3_cuts])
    seed = args[i, :seed]
    # println("$density, $n_customers, $n_charging, $T, $ngroute_neighborhood_charging_size, $use_lmSR3_cuts")
    cond = (
        (results.seed .== seed)
        .& (results.n_customers .== n_customers)
        .& (results.n_charging .== n_charging)
        .& (results.T .== T)
        .& (results.ngroute_neighborhood_charging_size .== ngroute_neighborhood_charging_size)
        # .& (results.use_lmSR3_cuts .== use_lmSR3_cuts)
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

# Plotting: aggregate grouped bar plots with StatsPlots
# For presentation: SRI cuts
begin
    T_range = 72000:9000:90000
    ngroute_neighborhood_charging_size_range = ["small", "medium"]
    groupnames = ["CG", "CG + ng-relaxation", "CG + ng-relaxation + SRI cuts"]
    sizenames = string.(20:4:40)
    for T in T_range, ngroute_neighborhood_charging_size in ngroute_neighborhood_charging_size_range
        data_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.T == T
                    && r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
                ),
                x
            )
        )
        time_taken_df = (
            data_df
            |> x -> filter(r -> r.use_lmSR3_cuts, x) 
            |> x -> select(x, :n_customers, :time_taken_first, :time_taken_total, :time_taken_total_SR3)
        )
        LP_IP_gap_df = (
            data_df
            |> x -> filter(r -> r.use_lmSR3_cuts, x) 
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
        yticks!(p1, p1_ytickvalues, string.(p1_ytickvalues))
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_SRI_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_SRI_$(T)_$ngroute_neighborhood_charging_size.png")
        p2_ytickvalues = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        p2 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(LP_IP_gap_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.75,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Optimality gap (%)",
            ylim = 10.0.^(-0.5, 2),
            grid = :y,
            legend = :bottomright,
        )
        yticks!(p2, p2_ytickvalues, string.(p2_ytickvalues))
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_SRI_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_SRI_$(T)_$ngroute_neighborhood_charging_size.png")
    end
end

# Plotting: aggregate grouped bar plots with StatsPlots
# For presentation: lm-SRI cuts
begin
    T_range = 72000:9000:90000
    ngroute_neighborhood_charging_size_range = ["small", "medium"]
    groupnames = ["CG", "CG + ng-relaxation", "CG + ng-relaxation + lmSRI cuts"]
    sizenames = string.(20:4:40)
    for T in T_range, ngroute_neighborhood_charging_size in ngroute_neighborhood_charging_size_range
        data_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.T == T
                    && r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
                ),
                x
            )
        )
        time_taken_df = (
            data_df
            |> x -> filter(r -> r.use_lmSR3_cuts, x) 
            |> x -> select(x, :n_customers, :time_taken_first, :time_taken_total, :time_taken_total_SR3)
        )
        LP_IP_gap_df = (
            data_df
            |> x -> filter(r -> r.use_lmSR3_cuts, x) 
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
        yticks!(p1, p1_ytickvalues, string.(p1_ytickvalues))
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_lmSRI_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_lmSRI_$(T)_$ngroute_neighborhood_charging_size.png")
        p2_ytickvalues = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        p2 = groupedbar(
            repeat(sizenames, outer = length(groupnames)),
            Matrix(LP_IP_gap_df[!, 2:end]),
            group = repeat(groupnames, inner = length(sizenames)),
            barwidth = 0.75,
            framestyle = :box,
            xlabel = "# tasks",
            yscale = :log10,
            ylabel = "Optimality gap (%)",
            ylim = 10.0.^(-0.5, 2),
            grid = :y,
            legend = :bottomright,
        )
        yticks!(p2, p2_ytickvalues, string.(p2_ytickvalues))
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_lmSRI_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_lmSRI_$(T)_$ngroute_neighborhood_charging_size.png")
    end
end


# Plotting: aggregate grouped bar plots with StatsPlots
begin
    T_range = 72000:9000:90000
    ngroute_neighborhood_charging_size_range = ["small", "medium"]
    groupnames = ["CG", "CG + ng-relaxation", "CG + ng-relaxation + SRI cuts", "CG + ng-relaxation + lmSRI cuts"]
    sizenames = string.(20:4:40)
    for T in T_range, ngroute_neighborhood_charging_size in ngroute_neighborhood_charging_size_range
        data_df = (
            summary 
            |> x -> filter(
                r -> (
                    r.T == T
                    && r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
                ),
                x
            )
        )
        time_taken_df = innerjoin(
            (
                data_df
                |> x -> filter(r -> !r.use_lmSR3_cuts, x) 
                |> x -> select(x, :n_customers, :time_taken_first, :time_taken_total)
            ),
            (
                data_df
                |> x -> unstack(x, :n_customers, :use_lmSR3_cuts, :time_taken_total_SR3)
                |> x -> select(
                    x,
                    :n_customers, 
                    "false" => :time_taken_total_SR3,
                    "true" => :time_taken_total_lmSR3,
                )
            ),
            on = :n_customers,
        )
        LP_IP_gap_df = innerjoin(
            (
                data_df
                |> x -> filter(r -> !r.use_lmSR3_cuts, x) 
                |> x -> select(x, :n_customers, :LP_IP_gap_first, :LP_IP_gap_last)
            ),
            (
                data_df
                |> x -> unstack(x, :n_customers, :use_lmSR3_cuts, :LP_IP_gap_last_SR3)
                |> x -> select(
                    x,
                    :n_customers, 
                    "false" => :LP_IP_gap_last_SR3,
                    "true" => :LP_IP_gap_last_lmSR3,
                )
            ),
            on = :n_customers,
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
            grid = :y,
            legend = :bottomright,
        )
        yticks!(p1, p1_ytickvalues, string.(p1_ytickvalues))
        hline!([3600], linestyle = :dash, label = false, color = :black)
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p1, "$(@__DIR__)/plots/time_taken_groupedbar_$(T)_$ngroute_neighborhood_charging_size.png")
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
            grid = :y,
            legend = :bottomright,
        )
        yticks!(p2, p2_ytickvalues, string.(p2_ytickvalues))
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p2, "$(@__DIR__)/plots/LP_IP_gap_groupedbar_$(T)_$ngroute_neighborhood_charging_size.png")
    end
end

# Plotting: aggregate 
begin

    μ = 5
    B = 15000
    T_range = 72000:9000:90000
    ngroute_neighborhood_charging_size_range = ["small", "medium"]
    
    density_range = collect(2.5:0.5:5.0)
    n_customers_range = Int.(density_range * 8)
    
    colors = Dict(
        n_customers => c
        for (n_customers, c) in zip(
            n_customers_range[end:-1:1],
            get(ColorSchemes.viridis, collect(0:1:(length(n_customers_range)-1)) ./(length(n_customers_range)-1))
        )
    )

    for T in T_range, 
        ngroute_neighborhood_charging_size in ngroute_neighborhood_charging_size_range
        xtickvalues = [1, 10, 100, 1000, 3600]
        ytickvalues = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]
        p = plot(
            xlabel = "Computational time (s)", 
            ylabel = "Optimality gap (%)",
            xscale = :log10,
            yscale = :log10,
            ylim = 10.0.^(-0.5, 2),
            xlim = 10.0.^(-0.5, 4.0),
        )
        vline!([3600], linestyle = :dash, label = false, color = :black)
        xticks!(p, xtickvalues, string.(xtickvalues))
        yticks!(p, ytickvalues, string.(ytickvalues))
        for n_customers in n_customers_range
            println("$n_customers, $T")
            data = (
                summary 
                |> x -> filter(
                    r -> (
                        r.ngroute_neighborhood_charging_size == ngroute_neighborhood_charging_size
                        && r.n_customers == n_customers
                        && r.T == T
                    ),
                    x
                )
            )
            Plots.scatter!(
                [data[1, :time_taken_first]],
                [data[1, :LP_IP_gap_first]],
                color = colors[n_customers],
                label = false,
                shape = :circle,
            )
            Plots.scatter!(
                [data[1, :time_taken_total]],
                [data[1, :LP_IP_gap_last]],
                color = colors[n_customers],
                label = false,
                shape = :circle,
            )
            Plots.scatter!(
                [data[1, :time_taken_total_SR3]],
                [data[1, :LP_IP_gap_last_SR3]],
                color = colors[n_customers],
                label = false,
                shape = :circle,
            )
            Plots.scatter!(
                [data[2, :time_taken_total_SR3]],
                [data[2, :LP_IP_gap_last_SR3]],
                color = colors[n_customers],
                label = false,
                shape = :square,
            )
            plot!(
                collect(data[1, [:time_taken_first, :time_taken_total]]),
                collect(data[1, [:LP_IP_gap_first, :LP_IP_gap_last]]),
                label = "$n_customers customers",
                color = colors[n_customers],
            )
            plot!(
                collect(data[1, [:time_taken_total, :time_taken_total_SR3]]),
                collect(data[1, [:LP_IP_gap_last, :LP_IP_gap_last_SR3]]),
                label = false,
                color = colors[n_customers],
            )
            plot!(
                [data[1, :time_taken_total], data[2, :time_taken_total_SR3]],
                [data[1, :LP_IP_gap_last], data[2, :LP_IP_gap_last_SR3]],
                label = false,
                color = colors[n_customers],
            )
        end
        plot!(legend = :topright)        
        savefig(p, "$(@__DIR__)/plots/time_taken_LP_IP_gap_scatter_$(T)_$ngroute_neighborhood_charging_size.pdf")
        savefig(p, "$(@__DIR__)/plots/time_taken_LP_IP_gap_scatter_$(T)_$ngroute_neighborhood_charging_size.png")
    end
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
        xlabel = "Computational time (s)",
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