using DataFrames, CSV
using DelimitedFiles

results = CSV.read("$(@__DIR__)/combined.csv", DataFrame)

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
        vcat(data_fields, [:seed], method_fields, :counter, :LP_IP_gap)
    )
    |> x -> groupby(
        x, 
        vcat(data_fields, method_fields)
    )
    |> x -> combine(
        x, 
        nrow, 
        :counter => mean,
        :LP_IP_gap => mean,
    )
    # |> x -> unstack(
    #     x,
    #     [:n_customers, :xmax,],
    #     :T,
    #     :nrow,
    # )
    # |> x -> sort(x, :xmax)
)

(
    summary
    |> x -> filter!(
        r -> !(r.xmax == 2.0 && r.T == 36000),
        x
    )
)

(
    summary
    |> x -> filter(
        r -> (
            r.elementary == false
            && r.ngroute == false
            # && r.ngroute_alt == true
        ),
        x
    )
    |> x -> unstack(
        x, 
        [:n_customers, :xmax, :T],
        :method,
        :time_taken,
    )
    |> x -> transform(
        x,
        [:benchmark, :ours] => ((x, y) -> x ./ y) => :ratio,
        [:n_customers, :xmax] => ((x, y) -> x ./ (2 * y)) => :customer_density,
    )
    |> x -> sort(
        x, 
        [:customer_density, :T],
    )
    |> x -> unstack(
        x, 
        :T,
        :customer_density,
        :ratio,
        combine = geomean,
    )
    |> x -> sort(
        x, 
        order(:T, rev = true),
    )
)