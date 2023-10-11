using DataFrames, CSV
using DelimitedFiles

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)
names(results)

data_fields = [
    :n_customers,
    :xmax,
    :T,
]
method_fields = [
    :method,
    :elementary,
    :ngroute,
    :ngroute_alt,
]

(
    results 
    |> x -> sort!(
        x, 
        vcat(data_fields, [:seed], method_fields),
    )
)

results.LP_IP_gap .*= 100

summary = (
    results 
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x, 
        :time_taken => geomean => :time_taken,
        :sp_time_taken_mean => geomean => :sp_time_taken_mean,
        :counter => mean => :counter,
        :converged => mean => :converged,
        :time_limit_reached => mean => :time_limit_reached,
        :LP_IP_gap => mean => :LP_IP_gap,
        :mean_subpath_length => mean => :mean_subpath_length,
        :mean_path_length => mean => :mean_path_length,
        :mean_ps_length => mean => :mean_ps_length,
    )
)

(
    summary
    |> x -> filter(
        r -> (
            r.xmax == 1.0
            && r.elementary == true
        ),
        x,
    )
    # |> unstack(
    #     x, 
    #     [:n_customers, :T],
    #     :method,
    #     :time_taken,
    # )
)
unstack(
    summary,
    ,
    ,
    :time_taken,
)

(
    results
    |> x -> filter(
        r -> (
            !r.converged
            # && r.xmax == 2.0
            # && r.T == 27000
            # && r.counter > 100
        ),
        x
    )
    |> x -> select(
        x, 
        vcat(data_fields, [:seed], method_fields, :counter)
    )
    |> x -> groupby(
        x, 
        data_fields
    )
    |> x -> combine(x, nrow)
    |> x -> unstack(
        x,
        [:n_customers, :xmax,],
        :T,
        :nrow,
    )
    |> x -> sort(x, :xmax)
)