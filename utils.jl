using Distances
using Random
using Suppressor 
using Combinatorics
using StatsBase
using StatsPlots
using Distributions
using Printf

function generate_instance(
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
    with_charging::Bool,
    B::Union{Float64, Nothing} = nothing,
    μ::Union{Float64, Nothing} = nothing,
    C::Int = 1,
)
    if with_charging
        if isnothing(B) || isnothing(μ)
            error()
        end
    end

    if with_charging
        n_nodes = n_depots + 2 * n_customers + n_charging * charging_repeats
    else
        n_nodes = n_depots + 2 * n_customers
    end

    N_pickups = collect(1:n_customers)
    N_dropoffs = collect(n_customers+1:2*n_customers)
    N_depots = collect(2*n_customers+1:2*n_customers+n_depots)
    N_vehicles = collect(1:n_vehicles)
    if with_charging
        N_charging = collect(2*n_customers+n_depots+1:2*n_customers+n_depots+n_charging*charging_repeats)
        N_nodes = vcat(N_pickups, N_dropoffs, N_depots, N_charging)
    else
        N_nodes = vcat(N_pickups, N_dropoffs, N_depots)
    end

    node_labels = merge(Dict(
        i => "Depot $ind" for (ind, i) in enumerate(N_depots)
    ), Dict(
        i => "Pickup $ind" for (ind, i) in enumerate(N_pickups)
    ), Dict(
        i => "Dropoff $ind" for (ind, i) in enumerate(N_dropoffs)
    ))
    if with_charging
        merge!(node_labels, Dict(
            i => "Charging $((ind+charging_repeats-1) ÷ charging_repeats)" for (ind, i) in enumerate(N_charging)
        ))
    end

    Random.seed!(seed)
    pickup_coords = Random.randn(Float64, (2, n_customers))
    dropoff_coords = Random.randn(Float64, (2, n_customers))
    depot_coords = shrinkage_depots * Random.randn(Float64, (2, n_depots))
    if with_charging
        charging_coords = shrinkage_charging * Random.randn(Float64, (2, n_charging))
        coords = hcat(
            pickup_coords, 
            dropoff_coords, 
            depot_coords, 
            repeat(
                charging_coords, 
                inner = (1, charging_repeats),
            )
        )
    else
        coords = hcat(
            pickup_coords, 
            dropoff_coords, 
            depot_coords, 
        )
    end

    distances = pairwise(Euclidean(), coords, dims=2)

    start_depots = sample(N_depots, n_vehicles, replace = true)
    V = Dict(i => findall(x -> x==i, start_depots) for i in N_depots)
    v_start = [length(V[i]) for i in N_depots]
    v_end_vec = repeat([1], n_depots)
    v_end  = Dict(i => v_end_vec[ind] for (ind, i) in enumerate(N_depots))
    
    c = Int.(round.(100 .* distances))
    t = Int.(round.(100 .* distances)) # travel times are integer
    if with_charging
        d = vcat(
            repeat([1], n_customers), 
            repeat([-1], n_customers), 
            repeat([0], length(N_depots) + length(N_charging)),
        )
        q = Int.(round.(100 .* distances)) # charge costs are integer
    else
        d = vcat(
            repeat([1], n_customers), 
            repeat([-1], n_customers), 
            repeat([0], length(N_depots)),
        )
    end

    if with_charging
        A = merge(
            Dict(
                (i,j) => t[i,j] for (i,j) in combinations(vcat(N_pickups, N_dropoffs, N_depots), 2)
            ), 
            Dict(
                (j,i) => t[i,j] for (i,j) in combinations(vcat(N_pickups, N_dropoffs, N_depots), 2)
            ), 
            # prevent self-loops at charging stations and between its duplicates
            Dict(
                (i,j) => t[i,j] for i in vcat(N_pickups, N_dropoffs, N_depots), j in N_charging
            ), 
            Dict(
                (j,i) => t[i,j] for i in vcat(N_pickups, N_dropoffs, N_depots), j in N_charging
            ), 
            Dict(
                (i,i) => 0 for i in N_depots # allow self-loops at depots
            ),
        )
    else
        A = merge(
            Dict(
                (i,j) => t[i,j] for (i,j) in combinations(N_nodes, 2)
            ), 
            Dict(
                (j,i) => t[i,j] for (i,j) in combinations(N_nodes, 2)
            ), 
            Dict(
                (i,i) => 0 for i in N_depots # allow self-loops at depots
            ),
        )
    end
    
    times_dist = Uniform(0, T)
    Random.seed!(seed)
    println(n_nodes)
    α = zeros(n_nodes)
    β = zeros(n_nodes)
    
    for i in 1:n_customers
        pickup_i = N_pickups[i]
        dropoff_i = N_dropoffs[i]
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
    α[N_depots] .= 0.0
    β[N_depots] .= T
    if with_charging
        α[N_charging] .= 0.0
        β[N_charging] .= T
    end

    data = Dict(
        "n_depots" => n_depots,
        "n_customers" => n_customers,
        "n_vehicles" => n_vehicles,
        "n_nodes" => n_nodes,
        "N_pickups" => N_pickups, 
        "N_dropoffs" => N_dropoffs,
        "N_depots" => N_depots,
        "N_vehicles" => N_vehicles,
        "N_nodes" => N_nodes,
        "node_labels" => node_labels,
        "shrinkage_depots" => shrinkage_depots,
        "pickup_coords" => pickup_coords,
        "dropoff_coords" => dropoff_coords,
        "depot_coords" => depot_coords,
        "coords" => coords,
        "distances" => distances,
        "V" => V,
        "v_start" => v_start,
        "v_end" => v_end,
        "c" => c,
        "t" => t,
        "C" => C,
        "d" => d,
        "T" => T,
        "A" => A,
        "α" => α,
        "β" => β,
    )
    if with_charging
        merge!(data, Dict(
            "n_charging" => n_charging,
            "charging_repeats" => charging_repeats,
            "N_charging" => N_charging,
            "shrinkage_charging" => shrinkage_charging,
            "charging_coords" => charging_coords,
            "q" => q,
            "B" => B,
            "μ" => μ,
        ))
    end
    return data
end

function plot_instance(data)
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

    plot!(legend = :outerright)
end

function results_printout(
    results, 
    params,
    data,
    with_charging::Bool = false,
)
    @printf("Objective:                 %7.1f\n", results["objective"])
    @printf("Total travel cost:         %7.1f\n", results["total_travel_cost"])
    @printf("Total vehicle time:        %7.1f\n", results["total_vehicle_time"])
    if with_charging
        @printf("Total charge required:     %7.1f\n", results["total_charge_required"])
    end
    @printf("Time taken:                %7.1f s\n", params["time_taken"])
    println("")

    if with_charging
        println("Vehicles         From      time   charge             To      time   charge  |  arc_time |  arc_charge ")
        println("----------------------------------------------------------------------------|-----------|-------------")
    else
        println("Vehicles         From      time             To      time  |  arc_time ")
        println("----------------------------------------------------------|-----------")
    end
    for k in data["N_vehicles"]
        arcs = [
            (i,j) for (i,j) in keys(data["A"]) 
            if results["x"][(i,j),k] > 0.5
        ]
        i = findfirst(x -> (x[1] ∈ data["N_depots"]), arcs)
        while true
            a = popat!(arcs, i)
            current_node = a[2]
            if with_charging
                if current_node in data["N_charging"]
                    @printf(
                        "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | -  %6.1f | -    %6.1f \n", 
                        k, 
                        data["node_labels"][a[1]], 
                        results["τ_leave"][a[1],k],
                        results["b_end"][a[1],k],
                        data["node_labels"][a[2]],
                        results["τ_reach"][a[2],k],
                        results["b_start"][a[2],k],
                        data["t"][a[1],a[2]],
                        data["q"][a[1],a[2]],
                    )
                    @printf(
                        "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | -  %6.1f | +    %6.1f \n", 
                        k, 
                        data["node_labels"][a[2]], 
                        results["τ_reach"][a[2],k],
                        results["b_start"][current_node,k],
                        data["node_labels"][a[2]], 
                        results["τ_leave"][a[2],k],
                        results["b_end"][current_node,k],
                        results["δ"][a[2],k],
                        results["δ"][a[2],k] * data["μ"],
                    )
                else
                    @printf(
                        "Vehicle %s: %12s (%6.1f,  %6.1f) -> %12s (%6.1f,  %6.1f) | -  %6.1f | -    %6.1f \n", 
                        k, 
                        data["node_labels"][a[1]], 
                        results["τ_leave"][a[1],k],
                        results["b_end"][a[1],k],
                        data["node_labels"][a[2]],
                        results["τ_reach"][a[2],k],
                        results["b_start"][a[2],k],
                        data["t"][a[1],a[2]],
                        data["q"][a[1],a[2]],
                    )
                end
            else
                @printf(
                    "Vehicle %s: %12s (%6.1f) -> %12s (%6.1f) | -  %6.1f\n", 
                    k, 
                    data["node_labels"][a[1]], 
                    results["τ_leave"][a[1],k],
                    data["node_labels"][a[2]],
                    results["τ_reach"][a[2],k],
                    data["t"][a[1],a[2]],
                )
            end
            if length(arcs) == 0
                break
            end
            i = findfirst(x -> (x[1] == current_node), arcs)
        end
        println("")
    end
end
