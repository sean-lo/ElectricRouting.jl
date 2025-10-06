include("utils.jl")

using CSV
using DataFrames

function transform_floats(
    x,
    pr::Float64 = 100.0,
    eps::Float64 = 1e-9,
)
    if (x * pr - floor.(x) * pr) < eps
        return Int.(floor.(x * pr))
    else
        return Int.(ceil.(x * pr))
    end
end

function read_evrptw_instance(
    fp::String,
    n_vehicles::Int,
    travel_cost_coeff::Int,
    charge_cost_coeff::Int,
    ;
    n_charging::Int = 0, # If 0, read from file
    n_customers::Int = 0, # If 0, read from file
    scale_time_horizon::Float64 = 1.0,
    scale_charge_capacity::Float64 = 1.0,
    scale_load_capacity::Float64 = 1.0,
)
    lines = readlines(fp)
    sep = findfirst(x -> occursin(r"^\s*$", x), lines)
    table_str = join(lines[1:sep-1], "\n")
    df = CSV.read(IOBuffer(table_str), DataFrame, delim=' ', ignorerepeated=true)

    (
        B_, # Unscaled battery capacity
        C,  # Load capacity
        _,
        inverse_refueling_rate,  # Inverse battery refueling rate
        _
    ) = [parse(Float64, match(r"/([\d.]+)/", x)[1]) for x in lines[sep+1:end]]

    # Scale load capacity
    C = Int(round(C * scale_load_capacity))

    # Scaled battery capacity
    B = transform_floats(B_ * scale_charge_capacity * inverse_refueling_rate)

    # DataFrames for depots, customers, charging stations
    Random.seed!(0)
    df = df[shuffle(1:nrow(df)), :]
    depots_df = filter(r -> r.Type == "d", df)
    customers_df = filter(r -> r.Type == "c", df)
    charging_df = filter(r -> r.Type == "f", df)

    if n_customers > 0
        customers_df = customers_df[1:n_customers, :]
    end
    if n_charging > 0
        charging_df = charging_df[1:n_charging, :]
    end

    # Index sets
    n_depots = nrow(depots_df)
    n_customers = nrow(customers_df)
    n_charging = nrow(charging_df)
    n_nodes = n_depots + n_customers + n_charging

    N_customers = 1:n_customers
    N_depots = n_customers+1:n_customers+n_depots
    N_vehicles = 1:n_vehicles
    N_charging = n_customers+n_depots+1:n_customers+n_depots+n_charging
    N_depots_charging = n_customers+1:n_customers+n_depots+n_charging
    N_nodes = 1:n_customers+n_depots+n_charging

    node_labels = merge(Dict(
        i => "Depot $ind" for (ind, i) in enumerate(N_depots)
    ), Dict(
        i => "Customer $ind" for (ind, i) in enumerate(N_customers)
    ), Dict(
        i => "Charging $ind" for (ind, i) in enumerate(N_charging)
    ))

    # Locations 
    depot_coords = depots_df[:, [:x, :y]] |> Matrix{Float64} |> transpose
    customer_coords = customers_df[:, [:x, :y]] |> Matrix{Float64} |> transpose
    charging_coords = charging_df[:, [:x, :y]] |> Matrix{Float64} |> transpose

    coords = hcat(customer_coords, depot_coords, charging_coords)
    (xmin, xmax) = extrema(coords[1, :])
    (ymin, ymax) = extrema(coords[2, :])

    distances = Distances.pairwise(Euclidean(), coords; dims=2)
    
    # Starting locations of vehicles
    start_depots = StatsBase.sample(N_depots, n_vehicles, replace = true)
    V = Dict(i => findall(x -> x==i, start_depots) for i in N_depots)
    v_start = v_end = Dict(i => length(V[i]) for i in N_depots)

    # Travel cost
    c = transform_floats.(distances)
    # (Scaled) travel time
    t = transform_floats.(distances * inverse_refueling_rate)
    # Battery consumption
    q = transform_floats.(distances)

    # Demand
    d = vcat(
        customers_df[!, :demand],
        depots_df[!, :demand],
        charging_df[!, :demand],
    )

    # Time windows
    α = vcat(
        customers_df[!, :ReadyTime],
        depots_df[!, :ReadyTime],
        charging_df[!, :ReadyTime],
    )
    α = transform_floats.(α * scale_time_horizon)
    β = vcat(
        customers_df[!, :DueDateServiceTime],
        depots_df[!, :DueDateServiceTime],
        charging_df[!, :DueDateServiceTime],
    )
    β = transform_floats.(β * scale_time_horizon)
    T = Int(round(maximum(β)))

    # Charging costs (no heterogenous charging)
    charge_cost_coeffs = Dict(
        i => charge_cost_coeff
        for i in N_charging
    )
    charge_cost_levelslist = [charge_cost_coeff]
    charge_cost_levels = Dict(
        i => 1
        for i in N_charging
    )
    charge_cost_nlevels = 1

    data = EVRPData(
        n_depots,
        n_customers,
        n_vehicles,
        n_charging,
        n_depots + n_charging,
        n_nodes,
        N_customers,
        N_depots,
        N_vehicles,
        N_charging,
        N_depots_charging,
        N_nodes,
        node_labels,
        "evrptw", # depot_pattern
        "evrptw", # customer_pattern
        "evrptw", # charging_pattern
        xmin, 
        xmax, 
        ymin,
        ymax, 
        customer_coords,
        depot_coords,
        charging_coords,
        coords,
        distances,
        V,
        v_start,
        v_end,
        c,
        t,
        q,
        d,
        C,
        T,
        α,
        β,
        inverse_refueling_rate,
        B,
        travel_cost_coeff,
        charge_cost_coeffs,
        charge_cost_levels,
        charge_cost_levelslist,
        charge_cost_nlevels,
    )

    return data
end
