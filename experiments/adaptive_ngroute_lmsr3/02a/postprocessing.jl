using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using StatsPlots
using ColorSchemes
using CairoMakie

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

ecdf_data = Dict(
    (tolerance, n_customers) => (
        all_results
        |> x -> filter(r -> r.n_customers == n_customers, x)
        |> x -> groupby(
            x, 
            :seed,
        )
        |> x -> select(
            x, 
            :seed,
            :CG_LP_IP_gap => :LP_IP_gap,
            :CG_time_taken => cumsum => :time_taken,
        )
        |> x -> groupby(
            x, 
            :seed,
        )
        |> x -> combine(
            x, 
            [:LP_IP_gap, :time_taken] 
            => ((x, y) -> y[something(findfirst(x .≤ tolerance), lastindex(x))])
            => :time_taken,
            :LP_IP_gap 
            => (x -> x[something(findfirst(x .≤ tolerance), lastindex(x))]) 
            => :LP_IP_gap,
        )
        |> x -> select(
            x, 
            :seed,
            :time_taken,
            :LP_IP_gap 
            => (x -> ifelse.(x .≤ tolerance, 0.0, x)) 
            => :LP_IP_gap,
        )
        |> x -> sort(
            x, 
            [:LP_IP_gap, :time_taken]
        )
    )
    for tolerance in [0.1, 1.0, 2.0, 5.0],
        n_customers in 20:4:48
)


function draw_ecdf(
    ecdf_data::Dict{Tuple{Float64, Int}, DataFrame},
    n_customers_range::Vector{Int},
    tolerance::Float64 = 0.0,
    filename::String = "ecdf",
)

    f = Figure(fontsize = 24)
    xtickvalues = [1, 10, 100, 1000]
    ax1 = Axis(
        f[1,1],
        xscale = log10,
        yticks = 0.0:0.2:1.0,
        limits = ((1, 3600), (0, 1)),
        xlabel = "Time taken (s)",
        xticks = (xtickvalues, string.(xtickvalues)),
    )
    ax2 = Axis(
        f[1,2],
        yticks = 0.0:0.2:1.0,
        yaxisposition = :right,
        limits = (0, nothing, 0, 1),
        xlabel = "Optimality gap (%)"
    )
    hideydecorations!(ax2, grid = false)
    linkyaxes!(ax1, ax2)
    colgap!(f.layout, 1, Relative(0.0))
    time_ecdfs = Combined{Makie.stairs, Tuple{Vector{Point{2, Float32}}}}[]
    gap_ecdfs = Combined{Makie.stairs, Tuple{Vector{Point{2, Float32}}}}[]
    for n_customers in n_customers_range
        data = copy(ecdf_data[(tolerance, n_customers)])
        data.time_taken[data.LP_IP_gap .> 0.0] .= 4000
        ind = nrow(data[data.LP_IP_gap .== 0.0, :]) + 1
        println(ind)
        data.time_taken[ind] = 3600
        push!(time_ecdfs, CairoMakie.ecdfplot!(ax1, data.time_taken, xscale = log10, overdraw = true,))
        push!(gap_ecdfs, CairoMakie.ecdfplot!(ax2, data.LP_IP_gap))
    end
    axislegend(ax2, gap_ecdfs, ["$n_customers customers" for n_customers in n_customers_range], position = :rb)
    save("$(@__DIR__)/plots/$filename.pdf", f)
    save("$(@__DIR__)/plots/$filename.png", f)
end

for tolerance in [0.1, 1.0, 2.0, 5.0]
    draw_ecdf(
        ecdf_data, 
        collect(skipmissing(unique(all_results.n_customers))),
        tolerance, 
        "ecdf_$tolerance",
    )
end

p = Plots.plot(
    xscale = :log10,
)
for n_customers in [24, 32, 40]
    Plots.scatter!(
        ecdf_data[(0.1, n_customers)][!, :time_taken],
        ecdf_data[(0.1, n_customers)][!, :LP_IP_gap],
        label = "$n_customers customers",
    )
end
savefig(p, "$(@__DIR__)/plots/scatter.pdf")