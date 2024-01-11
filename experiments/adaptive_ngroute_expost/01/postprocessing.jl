using DataFrames, CSV
# using DelimitedFiles
using StatsBase
# using StatsPlots
# using ColorSchemes

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
names(results)

# merge individual CSVs in results folder
begin
    args = CSV.read("$(@__DIR__)/args.csv", DataFrame)
    args.index = collect(1:nrow(args))
    all_dfs = DataFrame[]
    for ind in 1:nrow(args)
        data = CSV.read("$(@__DIR__)/results/$ind.csv", DataFrame)
        data.expost .= false
        data.instance .= ind
        data.iteration .= collect(1:nrow(data))
        push!(all_dfs, data)
        e_data = CSV.read("$(@__DIR__)/results/expost_$ind.csv", DataFrame)
        e_data.expost .= true
        e_data.instance .= ind
        e_data.iteration .= collect(1:nrow(e_data))
        push!(all_dfs, e_data)
    end
    all_data = vcat(all_dfs...)
    select!(all_data, [:expost, :instance, :iteration], Not([:expost, :instance, :iteration]))
    CSV.write("$(@__DIR__)/all_results.csv", all_data)
end

data_fields = [
    :n_customers,
    :xmax,
    :T,
    # :density,
]
method_fields = [
    :method,
    :ngroute_neighborhood_charging_size,
]

(
    results 
    |> x -> sort!(
        x, 
        vcat(data_fields,  [:seed], method_fields),
    )
    |> x -> transform!(
        x, 
        [:time_taken_first, :time_taken_first_expost] 
        => ((x, y) -> x ./ y) => :time_taken_first_ratio,
        [:time_taken_total, :time_taken_total_expost]
        => ((x, y) -> x ./ y) => :time_taken_total_ratio,
    )
)

results.LP_IP_gap_first .*= 100
results.LP_IP_gap_last .*= 100
results.LP_IP_gap_first_expost .*= 100
results.LP_IP_gap_last_expost .*= 100

(
    results
    |> x -> groupby(
        x,
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x,
        nrow,
        # :converged => mean => :converged,
        # :time_taken_first => geomean => :time_taken_first,
        # :time_taken_total => geomean => :time_taken_total,
        # :time_taken_first_expost => geomean => :time_taken_first_expost,
        # :time_taken_total_expost => geomean => :time_taken_total_expost,
        :time_taken_first_ratio => geomean => :time_taken_first_ratio,
        :time_taken_total_ratio => geomean => :time_taken_total_ratio,
        # :LP_IP_gap_first => mean => :LP_IP_gap_first,
        # :LP_IP_gap_last => mean => :LP_IP_gap_last,
        # :LP_IP_gap_first_expost => mean => :LP_IP_gap_first_expost,
        # :LP_IP_gap_last_expost => mean => :LP_IP_gap_last_expost,
        # :neighborhood_size_mean_first => mean => :neighborhood_size_mean_first,
        # :neighborhood_size_mean_last => mean => :neighborhood_size_mean_last,
        # :neighborhood_size_mean_first_expost => mean => :neighborhood_size_mean_first_expost,
        # :neighborhood_size_mean_last_expost => mean => :neighborhood_size_mean_last_expost,
        
    )
    |> x -> show(
        x,
        allrows = true,
        allcols = true,
    )
)