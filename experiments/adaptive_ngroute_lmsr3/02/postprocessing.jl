using DataFrames, CSV
using DelimitedFiles
using Glob
using StatsBase
using StatsPlots
using ColorSchemes
using CairoMakie

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
names(results)

args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
args.index = collect(1:nrow(args))

# # merge individual CSVs in results folder
# begin
#     all_dfs = DataFrame[]
#     for ind in 1:nrow(args)
#         fp = "$(@__DIR__)/results/$ind.csv"
#         !isfile(fp) && continue
#         data = CSV.read(fp, DataFrame)
#         data.instance .= ind
#         data.iteration .= collect(1:nrow(data))
#         push!(all_dfs, data)
#     end
#     all_data = vcat(all_dfs...)
#     select!(all_data, [:instance, :iteration], Not([:instance, :iteration]))
#     CSV.write("$(@__DIR__)/all_results.csv", all_data)
# end

data_fields = [
    :n_customers,
    :xmax,
    :T,
    # :density,
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

results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100
results.LP_IP_gap_last_SR3 .*= 100

results.LP_IP_gap_first .= max.(results.LP_IP_gap_first, 0)
results.LP_IP_gap_last .= max.(results.LP_IP_gap_last, 0)
results.LP_IP_gap_last_SR3 .= max.(results.LP_IP_gap_last_SR3, 0)

all_results = CSV.read("$(@__DIR__)/all_results.csv", DataFrame)
all_results.CG_LP_IP_gap .*= 100
all_results.CG_LP_IP_gap .= ifelse.(all_results.CG_LP_IP_gap .≤ 1e-10, 0.0, all_results.CG_LP_IP_gap)

all_results = outerjoin(all_results, args, on = :instance => :index, makeunique = true)
select!(all_results, Not(:method_1))



processed_data |> x -> filter(r -> r.n_customers == 20, x)

function draw_ecdf(
    all_results::DataFrame,
    n_customers_range::Vector{Int},
    tolerance::Float64 = 0.0,
    filename::String = "ecdf",
)

    processed_data = (
        all_results
        |> x -> groupby(
            x, 
            [:seed, :n_customers],
        )
        |> x -> select(
            x, 
            :seed,
            :n_customers,
            :CG_LP_IP_gap => :LP_IP_gap,
            :CG_time_taken => cumsum => :time_taken,
        )
        |> x -> groupby(
            x, [:seed, :n_customers]
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
            :n_customers,
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

    f = Figure()
    ax1 = Axis(
        f[1,1],
        xscale = log10,
        limits = ((1, 3600), (0, 1)),
        xlabel = "Time taken (s)"
    )
    ax2 = Axis(
        f[1,2],
        yaxisposition = :right,
        limits = (0, nothing, 0, 1),
        xlabel = "Optimality gap (%)"
    )
    linkyaxes!(ax1, ax2)
    colgap!(f.layout, 1, Relative(0.0))
    time_ecdfs = Combined{Makie.stairs, Tuple{Vector{Point{2, Float32}}}}[]
    gap_ecdfs = Combined{Makie.stairs, Tuple{Vector{Point{2, Float32}}}}[]
    for n_customers in n_customers_range
        data = (
            processed_data 
            |> x -> filter(
                r -> r.n_customers == n_customers, 
                x
            )
        ) 
        data.time_taken[data.LP_IP_gap .> 0.0] .= 4000.0
        println(data)
        push!(time_ecdfs, CairoMakie.ecdfplot!(ax1, data.time_taken, xscale = log10))
        push!(gap_ecdfs, CairoMakie.ecdfplot!(ax2, data.LP_IP_gap))
    end
    axislegend(ax2, gap_ecdfs, ["$n_customers customers" for n_customers in n_customers_range], position = :rb)
    save("$(@__DIR__)/plots/$filename.pdf", f)
    save("$(@__DIR__)/plots/$filename.png", f)
end

for tolerance in [0.1, 1.0, 2.0, 5.0]
    draw_ecdf(
        all_results, 
        collect(skipmissing(unique(all_results.n_customers))),
        tolerance, 
        "ecdf_$tolerance",
    )
end