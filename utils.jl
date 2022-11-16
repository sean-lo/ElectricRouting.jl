using Random
using Suppressor 
using Combinatorics
using StatsBase
using StatsPlots
using Distributions
using Distances
using Printf

function generate_times(
    T,
    n_customers,
    seed
)
    times_dist = Uniform(0, T)
    α = zeros(2 * n_customers)
    β = zeros(2 * n_customers)

    Random.seed!(seed)
    for i in 1:n_customers
        pickup_i = i
        dropoff_i = i + n_customers
        while true
            s1, s2, e1, e2 = round.(sort(rand(times_dist, 4)))
            if (e1 - s1) / T > 0.4 && (e2 - s2) / T > 0.4
                α[pickup_i] = s1
                α[dropoff_i] = s2
                β[pickup_i] = e1
                β[dropoff_i] = e2
                break
            end
        end
    end
    return α, β
end

function generate_instance_pair(
    ;
    n_depots::Int, 
    n_customers::Int,
    n_charging::Int,
    charging_repeats::Int,
    n_vehicles::Int,
    shrinkage_depots::Float64,
    shrinkage_charging::Float64,
    T::Float64,
    seed::Int,
    B::Float64,
    μ::Float64,
    C::Int = 1,
)
    n_nodes_charge = n_depots + 2 * n_customers + n_charging * charging_repeats
    n_nodes_nocharge = n_depots + 2 * n_customers

    seeds = abs.(rand(MersenneTwister(seed), Int, 6))

    N_pickups = collect(1:n_customers)
    N_dropoffs = collect(n_customers+1:2*n_customers)
    N_depots = collect(2*n_customers+1:2*n_customers+n_depots)
    N_vehicles = collect(1:n_vehicles)
    N_charging = collect(2*n_customers+n_depots+1:2*n_customers+n_depots+n_charging*charging_repeats)
    N_nodes_nocharge = vcat(N_pickups, N_dropoffs, N_depots)
    N_nodes_charge = vcat(N_pickups, N_dropoffs, N_depots, N_charging)

    node_labels_nocharge = merge(Dict(
        i => "Depot $ind" for (ind, i) in enumerate(N_depots)
    ), Dict(
        i => "Pickup $ind" for (ind, i) in enumerate(N_pickups)
    ), Dict(
        i => "Dropoff $ind" for (ind, i) in enumerate(N_dropoffs)
    ))
    node_labels_charge = merge(
        node_labels_nocharge, 
        Dict(
            i => "Charging $((ind+charging_repeats-1) ÷ charging_repeats)" for (ind, i) in enumerate(N_charging)
        )
    )

    pickup_coords = Random.randn(MersenneTwister(seeds[1]), Float64, (2, n_customers))
    dropoff_coords = Random.randn(MersenneTwister(seeds[2]), Float64, (2, n_customers))
    depot_coords = shrinkage_depots * Random.randn(MersenneTwister(seeds[3]), Float64, (2, n_depots))
    charging_coords = shrinkage_charging * Random.randn(MersenneTwister(seeds[4]), Float64, (2, n_charging))
    coords_nocharge = hcat(
        pickup_coords, 
        dropoff_coords, 
        depot_coords,  
    )
    coords_charge = hcat(
        coords_nocharge,
        repeat(
            charging_coords, 
            inner = (1, charging_repeats),
        )
    )

    distances_nocharge = Distances.pairwise(Euclidean(), coords_nocharge, dims=2)
    distances_charge = Distances.pairwise(Euclidean(), coords_charge, dims=2)

    start_depots = StatsBase.sample(MersenneTwister(seeds[5]), N_depots, n_vehicles, replace = true)
    V = Dict(i => findall(x -> x==i, start_depots) for i in N_depots)
    v_start = [length(V[i]) for i in N_depots]
    v_end_vec = repeat([1], n_depots)
    v_end  = Dict(i => v_end_vec[ind] for (ind, i) in enumerate(N_depots))
    
    c_nocharge = Int.(round.(100 .* distances_nocharge))
    c_charge = Int.(round.(100 .* distances_charge))
    t_nocharge = Int.(round.(100 .* distances_nocharge)) # travel times are integer
    t_charge = Int.(round.(100 .* distances_charge)) # travel times are integer
    d_nocharge = vcat(
        repeat([1], n_customers), 
        repeat([-1], n_customers), 
        repeat([0], n_depots),
    )
    d_charge = vcat(
        d_nocharge,
        repeat([0], length(N_charging)),
    )
    q = Int.(round.(100 .* distances_charge)) # charge costs are integer

    A_nocharge = merge(
        Dict(
            (i,j) => t_nocharge[i,j] for (i,j) in combinations(N_nodes_nocharge, 2)
        ), 
        Dict(
            (j,i) => t_nocharge[i,j] for (i,j) in combinations(N_nodes_nocharge, 2)
        ), 
        Dict(
            (i,i) => 0 for i in N_depots # allow self-loops at depots
        ),
    )
    A_charge = merge(
        Dict(
            (i,j) => t_charge[i,j] for (i,j) in combinations(vcat(N_pickups, N_dropoffs, N_depots), 2)
        ), 
        Dict(
            (j,i) => t_charge[i,j] for (i,j) in combinations(vcat(N_pickups, N_dropoffs, N_depots), 2)
        ), 
        # prevent self-loops at charging stations and between its duplicates
        Dict(
            (i,j) => t_charge[i,j] for i in vcat(N_pickups, N_dropoffs, N_depots), j in N_charging
        ), 
        Dict(
            (j,i) => t_charge[i,j] for i in vcat(N_pickups, N_dropoffs, N_depots), j in N_charging
        ), 
        Dict(
            (i,i) => 0 for i in N_depots # allow self-loops at depots
        ),
    )
    
    α, β = generate_times(T, n_customers, seeds[6])
    α_nocharge = vcat(α, repeat([0.0], n_depots))
    β_nocharge = vcat(β, repeat([T], n_depots))
    α_charge = vcat(α_nocharge, repeat([0.0], length(N_charging)))
    β_charge = vcat(β_nocharge, repeat([T], length(N_charging)))

    nocharge_data = Dict(
        "n_depots" => n_depots,
        "n_customers" => n_customers,
        "n_vehicles" => n_vehicles,
        "n_nodes" => n_nodes_nocharge,
        "N_pickups" => N_pickups, 
        "N_dropoffs" => N_dropoffs,
        "N_depots" => N_depots,
        "N_vehicles" => N_vehicles,
        "N_nodes" => N_nodes_nocharge,
        "node_labels" => node_labels_nocharge,
        "shrinkage_depots" => shrinkage_depots,
        "pickup_coords" => pickup_coords,
        "dropoff_coords" => dropoff_coords,
        "depot_coords" => depot_coords,
        "coords" => coords_nocharge,
        "distances" => distances_nocharge,
        "V" => V,
        "v_start" => v_start,
        "v_end" => v_end,
        "c" => c_nocharge,
        "t" => t_nocharge,
        "C" => C,
        "d" => d_nocharge,
        "T" => T,
        "A" => A_nocharge,
        "α" => α_nocharge,
        "β" => β_nocharge,
    )
    charge_data = merge(nocharge_data, Dict(
        "n_charging" => n_charging,
        "charging_repeats" => charging_repeats,
        "n_nodes" => n_nodes_charge,
        "N_charging" => N_charging,
        "N_nodes" => N_nodes_charge,
        "node_labels" => node_labels_charge,
        "shrinkage_charging" => shrinkage_charging,
        "charging_coords" => charging_coords,
        "coords" => coords_charge,
        "distances" => distances_charge,
        "c" => c_charge,
        "t" => t_charge,
        "d" => d_charge,
        "A" => A_charge,
        "α" => α_charge,
        "β" => β_charge,
        "q" => q,
        "B" => B,
        "μ" => μ,
    ))
    return nocharge_data, charge_data
end

function preprocess_arcs(
    in_data, 
    with_charging::Bool = false,
    allow_nonempty_at_charging::Bool = true,
) 
    data = copy(in_data)
    for (i,j) in keys(data["A"])
        # remove arcs due to time window infeasibility
        if data["α"][i] + data["t"][i,j] > data["β"][j]
            delete!(data["A"], (i,j))
        end
        if with_charging
            # remove arcs due to charge infeasibility
            if data["q"][i,j] > data["B"]
                delete!(data["A"], (i,j))
            end
        end
        # (depot, dropoff)-pair infeasible
        if i in data["N_depots"] && j in data["N_dropoffs"]
            delete!(data["A"], (i,j))
        end
        # (pickup, depot)-pair infeasible
        if i in data["N_pickups"] && j in data["N_depots"]
            delete!(data["A"], (i,j))
        end
        if with_charging 
            if !allow_nonempty_at_charging
                # (pickup, charging)-pair infeasible
                if i in data["N_pickups"] && j in data["N_charging"]
                    delete!(data["A"], (i,j))
                end
                # (pickup, charging)-pair infeasible
                if i in data["N_charging"] && j in data["N_dropoffs"]
                    delete!(data["A"], (i,j))
                end
            end
        end
        if with_charging
            # if C == 1: (pickup[i], j) infeasible if j != dropoff[i]
            if (i in data["N_pickups"]) && (j != i + data["n_customers"]) && !(j in data["N_charging"])
                delete!(data["A"], (i,j))
            end
            # if C == 1: (i, dropoff[j]) infeasible if i != pickup[j]
            if (j in data["N_dropoffs"]) && (i != j - data["n_customers"]) && !(i in data["N_charging"])
                delete!(data["A"], (i,j))
            end
        else
            # if C == 1: (pickup[i], j) infeasible if j != dropoff[i]
            if (i in data["N_pickups"]) && (j != i + data["n_customers"])
                delete!(data["A"], (i,j))
            end
            # if C == 1: (i, dropoff[j]) infeasible if i != pickup[j]
            if (j in data["N_dropoffs"]) && (i != j - data["n_customers"])
                delete!(data["A"], (i,j))
            end
        end
    end
    return data
end

function plot_instance(data, with_charging::Bool = true)
    p = plot(
        # xlim = (0, 1), ylim = (0, 1),
        aspect_ratio = :equal, 
        fmt = :png, 
    )
    plot!(
        data["pickup_coords"][1,:], data["pickup_coords"][2,:],
        seriestype = :scatter, 
        label = "Pickup",
        color = :green
    )
    annotate!.(
        data["pickup_coords"][1,:] .+ 0.01, data["pickup_coords"][2,:], 
        text.(
            collect(string(i) for i in 1:data["n_customers"]), 
            :green, :left, 11
        )
    )
    plot!(
        data["dropoff_coords"][1,:], data["dropoff_coords"][2,:],
        seriestype = :scatter, 
        label = "Dropoff",
        color = :red
    )
    annotate!.(
        data["dropoff_coords"][1,:] .+ 0.01, data["dropoff_coords"][2,:], 
        text.(
            collect(string(i) for i in 1:data["n_customers"]), 
            :red, :left, 11
        )
    )
    plot!(
        data["depot_coords"][1,:], data["depot_coords"][2,:],
        seriestype = :scatter, 
        label = "Depots",
        color = :black
    )
    annotate!.(
        data["depot_coords"][1,:] .+ 0.01, data["depot_coords"][2,:], 
        text.(
            collect("M" * string(i) for i in 1:data["n_depots"]), 
            :black, :left, 11
        )
    )
    if with_charging
        plot!(
            data["charging_coords"][1,:], data["charging_coords"][2,:],
            seriestype = :scatter, 
            label = "Recharging stations",
            color = :grey
        )
        annotate!.(
            data["charging_coords"][1,:] .+ 0.01, data["charging_coords"][2,:], 
            text.(
                collect("R" * string(i) for i in 1:data["n_charging"]), 
                :grey, :left, 11
            )
        )
    end

    plot!(legend = :outerright)
end