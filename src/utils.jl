using CompositeStructs
using Random
using Suppressor 
using Combinatorics
using StatsBase
using StatsPlots
using Distributions
using Distances
using Printf

using Graphs

Base.@kwdef mutable struct Subpath
    n_customers::Int
    starting_node::Int
    starting_time::Float64
    starting_charge::Float64
    current_node::Int = starting_node
    arcs::Vector{Tuple} = []
    time::Float64 = starting_time
    charge::Float64 = starting_charge
    served::BitVector = falses(n_customers)
    serve_times::Vector{Float64} = zeros(n_customers)
    delta_time::Float64 = 0.0
    delta_charge::Float64 = 0.0
    end_time::Float64 = time
    end_charge::Float64 = charge
    round_time::Float64 = end_time
    round_charge::Float64 = end_charge
    artificial::Bool = false
end

Base.copy(s::Subpath) = Subpath(
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    time = s.time,
    charge = s.charge,
    served = copy(s.served),
    serve_times = copy(s.serve_times),
    delta_time = s.delta_time,
    delta_charge = s.delta_charge,
    end_time = s.end_time,
    end_charge = s.end_charge,
    round_time = s.round_time,
    round_charge = s.round_charge,
    artificial = s.artificial
)

Base.show(io::IO, s::Subpath) = begin
    if s.artificial
        print(io, """Subpath (artificial):
        ($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.round_time), $(s.round_charge))
        """)
    else
        print(io, """Subpath:
        ($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.round_time), $(s.round_charge))
        arcs:           $(s.arcs)
        served:         $(s.served)
        serve_times:    $(s.serve_times)
        now:            ($(s.time), $(s.charge))
        delta:          ($(s.delta_time), $(s.delta_charge))
        end:            ($(s.end_time), $(s.end_charge))
        round:          ($(s.round_time), $(s.round_charge))
        """)
    end
end

Base.isequal(s1::Subpath, s2::Subpath) = begin 
    (
        s1.n_customers == s2.n_customers
        && s1.starting_node == s2.starting_node
        && s1.starting_time == s2.starting_time
        && s1.starting_charge == s2.starting_charge
        && s1.current_node == s2.current_node
        && s1.arcs == s2.arcs
        && s1.time == s2.time
        && s1.charge == s2.charge
        && s1.served == s2.served
        && s1.serve_times == s2.serve_times
        && s1.delta_time == s2.delta_time
        && s1.delta_charge == s2.delta_charge
        && s1.end_time == s2.end_time
        && s1.end_charge == s2.end_charge
        && s1.round_time == s2.round_time
        && s1.round_charge == s2.round_charge
        && s1.artificial == s2.artificial
    )
end

@composite Base.@kwdef mutable struct SubpathWithCost 
    Subpath...
    cost::Float64 = 0.0
    explored::Bool = false
end

Base.show(io::IO, s::SubpathWithCost) = begin 
    if s.artificial
        print(io, """SubpathWithCost (artificial):
        ($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.round_time), $(s.round_charge))
        """)
    else
        print(io, """SubpathWithCost:
        ($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.round_time), $(s.round_charge))
        cost:           $(s.cost)
        arcs:           $(s.arcs)
        served:         $(s.served)
        serve_times:    $(s.serve_times)
        now:            ($(s.time), $(s.charge))
        delta:          ($(s.delta_time), $(s.delta_charge))
        end:            ($(s.end_time), $(s.end_charge))
        round:          ($(s.round_time), $(s.round_charge))
        """)
    end
end

Base.copy(s::SubpathWithCost) = SubpathWithCost(
    cost = s.cost,
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    time = s.time,
    charge = s.charge,
    served = copy(s.served),
    serve_times = copy(s.serve_times),
    delta_time = s.delta_time,
    delta_charge = s.delta_charge,
    end_time = s.end_time,
    end_charge = s.end_charge,
    round_time = s.round_time,
    round_charge = s.round_charge,
    artificial = s.artificial,
)

Base.isequal(s1::SubpathWithCost, s2::SubpathWithCost) = begin 
    (
        s1.cost == s2.cost
        && s1.n_customers == s2.n_customers
        && s1.starting_node == s2.starting_node
        && s1.starting_time == s2.starting_time
        && s1.starting_charge == s2.starting_charge
        && s1.current_node == s2.current_node
        && s1.arcs == s2.arcs
        && s1.time == s2.time
        && s1.charge == s2.charge
        && s1.served == s2.served
        && s1.serve_times == s2.serve_times
        && s1.delta_time == s2.delta_time
        && s1.delta_charge == s2.delta_charge
        && s1.end_time == s2.end_time
        && s1.end_charge == s2.end_charge
        && s1.round_time == s2.round_time
        && s1.round_charge == s2.round_charge
        && s1.artificial == s2.artificial
    )
end

Subpath(s::SubpathWithCost) = Subpath(
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    time = s.time,
    charge = s.charge,
    served = copy(s.served),
    serve_times = copy(s.serve_times),
    delta_time = s.delta_time,
    delta_charge = s.delta_charge,
    end_time = s.end_time,
    end_charge = s.end_charge,
    round_time = s.round_time,
    round_charge = s.round_charge,
    artificial = s.artificial,
)

Base.@kwdef mutable struct Path
    subpaths::Vector{Subpath}
end

Base.isequal(p1::Path, p2::Path) = all(isequal(s1, s2) for (s1, s2) in zip(p1.subpaths, p2.subpaths))

Base.@kwdef mutable struct PathWithCost
    subpaths::Vector{SubpathWithCost}
    cost::Float64
    explored::Bool = false
end

Path(p::PathWithCost) = Path(
    subpaths = [Subpath(s) for s in p.subpaths],
)

function dceil(
    x::Float64,
    points,
)
    return points[searchsortedfirst(points, x)]
end

function dfloor(
    x::Float64,
    points,
)
    return points[searchsortedlast(points, x)]
end

function dceilall(
    x::Float64,
    points,
)
    return points[searchsortedfirst(points, x):end]
end

function dfloorall(
    x::Float64,
    points,
)
    return points[1:searchsortedlast(points, x)]
end

function generate_times(
    T::Float64,
    n_customers::Int,
    seed::Int,
    batch::Int,
    permissiveness::Float64 = 0.4,
)
    if n_customers % batch != 0
        error()
    end
    times_dist = Uniform(0.0, T)
    α = zeros(2 * n_customers)
    β = zeros(2 * n_customers)

    Random.seed!(seed)
    for batch_ind in 1:(n_customers ÷ batch)
        pickup_inds = collect((batch_ind-1)*batch+1:batch_ind*batch)
        dropoff_inds = collect(n_customers+(batch_ind-1)*batch+1:n_customers+batch_ind*batch)
        while true
            times = round.(sort(rand(times_dist, 4 * batch)))
            start_times = times[1:end÷2]
            end_times = times[end÷2+1:end]
            if all(
                (end_times .- start_times) ./ T .> permissiveness
            )
                α[pickup_inds] = start_times[1:end÷2]
                α[dropoff_inds] = start_times[end÷2+1:end]
                β[pickup_inds] = end_times[1:end÷2]
                β[dropoff_inds] = end_times[end÷2+1:end]
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
    batch::Int = 1,
    permissiveness::Float64 = 0.4,
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
            (i,j) => t_charge[i,j] for (i,j) in combinations(N_charging, 2)
        ),
        Dict(
            (j,i) => t_charge[i,j] for (i,j) in combinations(N_charging, 2)
        ), # FIXME: what happens to duplicates of charging stations?
        Dict(
            (i,i) => 0 for i in N_depots # allow self-loops at depots
        ),
    )
    
    α, β = generate_times(T, n_customers, seeds[6], batch, permissiveness)
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

function construct_graph(data)
    G = SimpleDiGraph(data["n_nodes"])
    for (i, j) in keys(data["A"])
        add_edge!(G, i, j)
    end
    return G
end

function preprocess_arcs(
    in_data, 
    with_charging::Bool = false,
    allow_nonempty_at_charging::Bool = true,
) 
    data = deepcopy(in_data)
    # Shrinking of time windows (Dumas 1991)
    for j in data["N_dropoffs"]
        data["β"][j] = min(
            data["β"][j],
            maximum(data["β"][data["N_depots"]] .- data["t"][j,data["N_depots"]])
        )
    end
    for i in data["N_pickups"]
        j = i + data["n_customers"]
        data["β"][i] = min(
            data["β"][i], 
            data["β"][j] - data["t"][i,j],
        )
    end
    for i in data["N_pickups"]
        data["α"][i] = max(
            data["α"][i],
            minimum(data["α"][data["N_depots"]] + data["t"][data["N_depots"],i])
        )
    end
    for j in data["N_dropoffs"]
        i = j - data["n_customers"]
        data["α"][j] = max(
            data["α"][j],
            data["α"][i] + data["t"][i,j],
        )
    end

    for (i,j) in keys(data["A"])
        # Remove going from a dropoff to its pickup
        if i in data["N_dropoffs"] && j == i - data["n_customers"]
            delete!(data["A"], (i,j))
        end
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
                # (charging, dropoff)-pair infeasible
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