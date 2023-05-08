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
    current_time::Float64 = starting_time
    current_charge::Float64 = starting_charge
    served::BitVector = falses(n_customers)
    artificial::Bool = false
end

Base.copy(s::Subpath) = Subpath(
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    current_time = s.current_time,
    current_charge = s.current_charge,
    served = copy(s.served),
    artificial = s.artificial
)

Base.show(io::IO, s::Subpath) = begin
    if s.artificial
        message = """Subpath (artificial):
        """
    else
        message = """Subpath:
        """
    end
    message = message * """
    ($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.current_time), $(s.current_charge))
    arcs:           $(s.arcs)
    served:         $(s.served)
    """
    print(io, message)
end

Base.isequal(s1::Subpath, s2::Subpath) = begin 
    (
        s1.n_customers == s2.n_customers
        && s1.starting_node == s2.starting_node
        && s1.starting_time == s2.starting_time
        && s1.starting_charge == s2.starting_charge
        && s1.current_node == s2.current_node
        && s1.arcs == s2.arcs
        && s1.current_time == s2.current_time
        && s1.current_charge == s2.current_charge
        && s1.served == s2.served
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
        message = """SubpathWithCost (artificial):
        """
    else
        message = """SubpathWithCost:
        """
    end
    message = message * """
    ($(s.starting_node), $(s.starting_time), $(s.starting_charge)) -> ($(s.current_node), $(s.current_time), $(s.current_charge))
    cost:           $(s.cost)
    arcs:           $(s.arcs)
    served:         $(s.served)
    """
    print(io, message)
end

Base.copy(s::SubpathWithCost) = SubpathWithCost(
    cost = s.cost,
    n_customers = s.n_customers,
    starting_node = s.starting_node,
    starting_time = s.starting_time,
    starting_charge = s.starting_charge,
    current_node = s.current_node,
    arcs = copy(s.arcs),
    current_time = s.current_time,
    current_charge = s.current_charge,
    served = copy(s.served),
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
        && s1.current_time == s2.current_time
        && s1.current_charge == s2.current_charge
        && s1.served == s2.served
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
    current_time = s.current_time,
    current_charge = s.current_charge,
    served = copy(s.served),
    artificial = s.artificial,
)

Base.@kwdef mutable struct ChargingArc
    starting_node::Int
    starting_time::Float64
    starting_charge::Float64
    delta_time::Float64 = 0.0
    delta_charge::Float64 = 0.0
    current_time::Float64 = starting_time
    current_charge::Float64 = starting_charge
end

Base.isequal(a1::ChargingArc, a2::ChargingArc) = (
    a1.starting_node == a2.starting_node
    && a1.starting_time == a2.starting_time
    && a1.starting_charge == a2.starting_charge
    && a1.delta_time == a2.delta_time
    && a1.delta_charge == a2.delta_charge
    && a1.current_time == a2.current_time
    && a1.current_charge == a2.current_charge
)

Base.copy(a::ChargingArc) = ChargingArc(
    starting_node = a.starting_node,
    starting_time = a.starting_time,
    starting_charge = a.starting_charge,
    delta_time = a.delta_time,
    delta_charge = a.delta_charge,
    current_time = a.current_time,
    current_charge = a.current_charge,
)

Base.@kwdef mutable struct Path
    subpaths::Vector{Subpath}
    charging_arcs::Vector{ChargingArc}
    served::Vector{Int} = sum(s.served for s in subpaths)
end

Base.isequal(p1::Path, p2::Path) = (
    all(isequal(s1, s2) for (s1, s2) in zip(p1.subpaths, p2.subpaths))
    && all(isequal(a1, a2) for (a1, a2) in zip(p1.charging_arcs, p2.charging_arcs))
    && p1.served == p2.served
)

Base.copy(p::Path) = Path(
    subpaths = [copy(s) for s in p.subpaths],
    charging_arcs = [copy(a) for a in p.charging_arcs],
    served = copy(p.served),
)

Base.@kwdef mutable struct PathWithCost
    subpaths::Vector{SubpathWithCost}
    charging_arcs::Vector{ChargingArc}
    cost::Float64
    served::Vector{Int} = sum(s.served for s in subpaths)
    explored::Bool = false
end

Base.isequal(p1::PathWithCost, p2::PathWithCost) = (
    all(isequal(s1, s2) for (s1, s2) in zip(p1.subpaths, p2.subpaths))
    && all(isequal(a1, a2) for (a1, a2) in zip(p1.charging_arcs, p2.charging_arcs))
    && p1.served == p2.served
    && p1.cost == p2.cost
    && p1.explored == p2.explored
)

Base.copy(p::PathWithCost) = PathWithCost(
    subpaths = [copy(s) for s in p.subpaths],
    charging_arcs = [copy(a) for a in p.charging_arcs],
    served = copy(p.served),
    cost = p.cost,
    explored = p.explored,
)

Path(p::PathWithCost) = Path(
    subpaths = [Subpath(s) for s in p.subpaths],
    charging_arcs = p.charging_arcs,
    served = sum(s.served for s in p.subpaths)
)

function dceil(
    x::Float64,
    points,
)
    return points[searchsortedfirst(points, x)]
end

function generate_instance(
    ;
    n_depots::Int, 
    n_customers::Int,
    n_charging::Int,
    n_vehicles::Int,
    shrinkage_depots::Float64,
    shrinkage_charging::Float64,
    T::Float64,
    seed::Int,
    B::Float64,
    μ::Float64,
    travel_cost_coeff::Int,
    charge_cost_coeff::Int,
    load_scale::Float64,
    load_shape::Float64,
    load_tolerance::Float64,
)
    n_nodes = n_depots + n_customers + n_charging

    seeds = abs.(rand(MersenneTwister(seed), Int, 6))

    N_customers = collect(1:n_customers)
    N_depots = collect(n_customers+1:n_customers+n_depots)
    N_vehicles = collect(1:n_vehicles)
    N_charging = collect(n_customers+n_depots+1:n_customers+n_depots+n_charging)
    N_nodes = vcat(N_customers, N_depots, N_charging)

    node_labels = merge(Dict(
        i => "Depot $ind" for (ind, i) in enumerate(N_depots)
    ), Dict(
        i => "Customer $ind" for (ind, i) in enumerate(N_customers)
    ), Dict(
        i => "Charging $ind" for (ind, i) in enumerate(N_charging)
    ))

    customer_coords = Random.randn(MersenneTwister(seeds[1]), Float64, (2, n_customers))
    depot_coords = shrinkage_depots * Random.randn(MersenneTwister(seeds[2]), Float64, (2, n_depots))
    charging_coords = shrinkage_charging * Random.randn(MersenneTwister(seeds[3]), Float64, (2, n_charging))
    coords = hcat(
        customer_coords,
        depot_coords,
        charging_coords,
    )

    distances = Distances.pairwise(Euclidean(), coords, dims=2)

    start_depots = StatsBase.sample(MersenneTwister(seeds[4]), N_depots, n_vehicles, replace = true)
    V = Dict(i => findall(x -> x==i, start_depots) for i in N_depots)
    v_start = [length(V[i]) for i in N_depots]
    v_end_vec = repeat([1], n_depots)
    v_end  = Dict(i => v_end_vec[ind] for (ind, i) in enumerate(N_depots))
    
    c = Int.(round.(100 .* distances))
    t = Int.(round.(100 .* distances)) # travel times are integer
    q = Int.(round.(100 .* distances)) # charge costs are integer
    d = vcat(
        floor.(rand(Gamma(load_scale, load_shape), n_customers) ./ n_customers),
        repeat([0], n_depots + n_charging),
    )
    C = ceil(sum(d) * load_tolerance / n_vehicles)

    A = merge(
        Dict(
            (i,i) => 0 for i in N_depots # allow self-loops at depots
        ),
        Dict(
            (i,j) => t[i,j] for (i,j) in permutations(vcat(N_customers, N_charging, N_depots), 2)
        ),
    )

    data = Dict(
        "n_depots" => n_depots,
        "n_customers" => n_customers,
        "n_vehicles" => n_vehicles,
        "n_charging" => n_charging,
        "n_nodes" => n_nodes,

        "N_customers" => N_customers,
        "N_depots" => N_depots,
        "N_vehicles" => N_vehicles,
        "N_charging" => N_charging,
        "N_nodes" => N_nodes,

        "node_labels" => node_labels,
        "shrinkage_depots" => shrinkage_depots,
        "shrinkage_charging" => shrinkage_charging,

        "customer_coords" => customer_coords,
        "depot_coords" => depot_coords,
        "charging_coords" => charging_coords,
        "coords" => coords,
        "distances" => distances,

        "V" => V,
        "v_start" => v_start,
        "v_end" => v_end,
        "c" => c,
        "t" => t,
        "q" => q,
        "d" => d,
        "C" => C,
        "T" => T,
        "A" => A,
        "μ" => μ,
        "B" => B,
        "travel_cost_coeff" => travel_cost_coeff,
        "charge_cost_coeff" => charge_cost_coeff,
    )
    return data
end

function construct_graph(data)
    G = SimpleDiGraph(data["n_nodes"])
    for (i, j) in keys(data["A"])
        add_edge!(G, i, j)
    end
    return G
end

function plot_instance(data)
    p = plot(
        # xlim = (0, 1), ylim = (0, 1),
        aspect_ratio = :equal, 
        fmt = :png, 
    )
    plot!(
        data["customer_coords"][1,:], data["customer_coords"][2,:],
        seriestype = :scatter, 
        label = "Customer",
        color = :green
    )
    annotate!.(
        data["customer_coords"][1,:] .+ 0.15, data["customer_coords"][2,:], 
        text.(
            collect(string(i) for i in 1:data["n_customers"]), 
            :green, :left, 11
        )
    )
    plot!(
        data["depot_coords"][1,:], data["depot_coords"][2,:],
        seriestype = :scatter, 
        label = "Depots",
        color = :black
    )
    annotate!.(
        data["depot_coords"][1,:] .+ 0.15, data["depot_coords"][2,:], 
        text.(
            collect("M" * string(i) for i in 1:data["n_depots"]), 
            :black, :left, 11
        )
    )
    plot!(
        data["charging_coords"][1,:], data["charging_coords"][2,:],
        seriestype = :scatter, 
        label = "Charging stations",
        color = :grey
    )
    annotate!.(
        data["charging_coords"][1,:] .+ 0.15, data["charging_coords"][2,:], 
        text.(
            collect("R" * string(i) for i in 1:data["n_charging"]), 
            :grey, :left, 11
        )
    )

    plot!(legend = :outerright)
    return p
end