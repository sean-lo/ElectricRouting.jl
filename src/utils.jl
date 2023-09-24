using Random
using Combinatorics
using StatsBase
using Distributions
using Distances
using Printf
using Graphs

Base.@kwdef mutable struct Subpath
    n_customers::Int
    starting_node::Int
    starting_time::Int
    starting_charge::Int
    current_node::Int = starting_node
    arcs::Vector{NTuple{2, Int}} = []
    current_time::Int = starting_time
    current_charge::Int = starting_charge
    served::Vector{Int} = zeros(Int, n_customers)
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
    artificial = s.artificial,
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

Base.@kwdef mutable struct ChargingArc
    starting_node::Int
    starting_time::Int
    starting_charge::Int
    delta::Int = 0
    current_time::Int = starting_time
    current_charge::Int = starting_charge
end

Base.isequal(a1::ChargingArc, a2::ChargingArc) = (
    a1.starting_node == a2.starting_node
    && a1.starting_time == a2.starting_time
    && a1.starting_charge == a2.starting_charge
    && a1.delta == a2.delta
    && a1.current_time == a2.current_time
    && a1.current_charge == a2.current_charge
)

Base.copy(a::ChargingArc) = ChargingArc(
    starting_node = a.starting_node,
    starting_time = a.starting_time,
    starting_charge = a.starting_charge,
    delta = a.delta,
    current_time = a.current_time,
    current_charge = a.current_charge,
)

Base.@kwdef mutable struct Path
    subpaths::Vector{Subpath}
    charging_arcs::Vector{ChargingArc}
    served::Vector{Int} = sum(s.served for s in subpaths)
    arcs::Vector{NTuple{2, Int}} = vcat([s.arcs for s in subpaths]...)
    customer_arcs::Vector{NTuple{2, Int}} = NTuple{2, Int}[]
end

Base.isequal(p1::Path, p2::Path) = (
    all(isequal(s1, s2) for (s1, s2) in zip(p1.subpaths, p2.subpaths))
    && all(isequal(a1, a2) for (a1, a2) in zip(p1.charging_arcs, p2.charging_arcs))
    && p1.served == p2.served
    && p1.arcs == p2.arcs
    && p1.customer_arcs == p2.customer_arcs
)

Base.copy(p::Path) = Path(
    subpaths = [copy(s) for s in p.subpaths],
    charging_arcs = [copy(a) for a in p.charging_arcs],
    served = copy(p.served),
    arcs = copy(p.arcs),
    customer_arcs = copy(p.customer_arcs),
)

function get_nodes(
    p::Path,
)
    return vcat(
        [a[1] for a in p.arcs],
        [p.arcs[end][2]]
    )
end

struct EVRPData
    n_depots::Int
    n_customers::Int
    n_vehicles::Int
    n_charging::Int
    n_depots_charging::Int
    n_nodes::Int
    N_customers::UnitRange{Int}
    N_depots::UnitRange{Int}
    N_vehicles::UnitRange{Int}
    N_charging::UnitRange{Int}
    N_depots_charging::UnitRange{Int}
    N_nodes::UnitRange{Int}
    node_labels::Dict{Int, String}
    shrinkage_depots::Float64
    shrinkage_charging::Float64
    customer_coords::Array{Float64, 2}
    depot_coords::Array{Float64, 2}
    charging_coords::Array{Float64, 2}
    coords::Array{Float64, 2}
    distances::Array{Float64, 2}
    V::Dict{Int, Vector{Int}}
    v_start::Vector{Int}
    v_end::Dict{Int, Int}
    c::Array{Int, 2}
    t::Array{Int, 2}
    q::Array{Int, 2}
    d::Vector{Int}
    C::Int
    T::Int
    α::Vector{Int}
    β::Vector{Int}
    μ::Int
    B::Int
    travel_cost_coeff::Int
    charge_cost_coeff::Int
end

struct EVRPGraph
    G::SimpleDiGraph{Int}
    node_labels::Dict{Int, String}
    c::Array{Int, 2}
    t::Array{Int, 2}
    q::Array{Int, 2}
    N_customers::UnitRange{Int}
    n_customers::Int
    N_depots::UnitRange{Int}
    n_depots::Int
    N_charging::UnitRange{Int}
    n_charging::Int
    N_depots_charging::UnitRange{Int}
    n_depots_charging::Int
    N_nodes::UnitRange{Int}
    n_nodes::Int
    A::Set{Tuple{Int, Int}}
    T::Int
    B::Int
    μ::Int
    α::Vector{Int}
    β::Vector{Int}
    min_t::Vector{Int}
    min_q::Vector{Int}
end

struct TimeLimitException <: Exception end

function add_message!(
    printlist::Vector{String}, 
    message::String, 
    verbose::Bool,
)
    push!(printlist, message)
    if verbose
        print(message)
    end
end

function compute_subpath_cost(
    data::EVRPData,
    graph::EVRPGraph,
    s::Subpath,
    M::Int = Int(1e10),
    ;
)
    if s.artificial 
        return M
    elseif length(s.arcs) == 0
        return 0
    end

    travel_cost = data.travel_cost_coeff * sum(
        graph.c[a...] for a in s.arcs
    )
    return travel_cost
end

function compute_subpath_modified_cost(
    data::EVRPData,
    graph::EVRPGraph,
    s::Subpath,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    verbose = false,
)
    reduced_cost = compute_subpath_cost(data, graph, s)
    verbose && @printf("Subpath cost: \t\t%11.3f\n", reduced_cost)

    service_cost = 0.0
    for (j, c) in enumerate(s.served)
        service_cost += (c * -ν[j])
    end
    verbose && @printf("Service cost: \t\t%11.3f\n", service_cost)
    reduced_cost += service_cost

    if s.starting_node in graph.N_depots
        if s.starting_time == 0.0 && s.starting_charge == graph.B
            verbose && @printf("Starting depot cost: \t%11.3f\n", (- κ[s.starting_node]))
            reduced_cost = reduced_cost - κ[s.starting_node]
        end
    end

    if s.current_node in graph.N_depots
        verbose && @printf("Ending depot cost: \t%11.3f\n", ( - μ[s.current_node]))
        reduced_cost = reduced_cost - μ[s.current_node]
    end

    verbose && @printf("Total modified cost: \t%11.3f\n\n", reduced_cost)

    return reduced_cost
end

function compute_subpath_costs(
    data::EVRPData,
    graph::EVRPGraph,
    all_subpaths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Subpath},
    },
    M::Int = Int(1e10),
    ;
)
    subpath_costs = Dict(
        key => Int[
            compute_subpath_cost(data, graph, s, M;)
            for s in all_subpaths[key]
        ]
        for key in keys(all_subpaths)
    )
    return subpath_costs
end

function compute_subpath_service(
    graph::EVRPGraph,
    all_subpaths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Subpath},
    },
)
    subpath_service = Dict(
        (key, i) => Int[
            s.served[i]
            for s in all_subpaths[key]
        ]
        for key in keys(all_subpaths), i in 1:graph.n_customers
    )
    return subpath_service
end

function compute_charging_arc_cost(
    data::EVRPData,
    a::ChargingArc,
)
    return data.charge_cost_coeff * a.delta
end

function compute_path_cost(
    data::EVRPData,
    graph::EVRPGraph, 
    p::Path,
    M::Int = Int(1e10),
    ;
    verbose = false,
)
    subpath_costs = length(p.subpaths) > 0 ? sum(compute_subpath_cost(data, graph, s, M) for s in p.subpaths) : 0
    verbose && @printf("Subpath costs: \t\t%11.3f\n", subpath_costs)

    charging_arc_costs = length(p.charging_arcs) > 0 ? sum(compute_charging_arc_cost(data, a) for a in p.charging_arcs) : 0
    verbose && @printf("Charging arc costs: \t\t%11d\n", charging_arc_costs)
    
    return subpath_costs + charging_arc_costs
end

function compute_path_modified_cost(
    data::EVRPData,
    graph::EVRPGraph,
    p::Path,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    ;
    verbose = false,
)
    reduced_cost = 0.0
    for s in p.subpaths
        reduced_cost += compute_subpath_modified_cost(data, graph, s, κ, μ, ν, verbose = verbose)
    end
    charging_costs = length(p.charging_arcs) > 0 ? sum(compute_charging_arc_cost(data, a) for a in p.charging_arcs) : 0
    verbose && @printf("Charging arc costs: \t%11d\n", charging_costs)

    reduced_cost += charging_costs

    verbose && @printf("Total modified cost: \t%11.3f\n\n", reduced_cost)

    return reduced_cost
end

function compute_path_costs(
    data::EVRPData,
    graph::EVRPGraph,
    all_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Path},
    },
    M::Int = Int(1e10),
    ;
)
    path_costs = Dict(
        key => Int[
            compute_path_cost(data, graph, p, M;)
            for p in all_paths[key]
        ]
        for key in keys(all_paths)
    )
    return path_costs
end

function compute_path_service(
    graph::EVRPGraph,
    all_paths::Dict{
        Tuple{NTuple{3, Int}, NTuple{3, Int}}, 
        Vector{Path},
    },
)
    path_service = Dict(
        (key, i) => Int[
            p.served[i]
            for p in all_paths[key]
        ]
        for key in keys(all_paths), i in 1:graph.n_customers
    )
    return path_service
end

function generate_instance(
    ;
    n_depots::Int, 
    n_customers::Int,
    n_charging::Int,
    n_vehicles::Int,
    depot_pattern::String,
    customer_pattern::String,
    charging_pattern::String,
    shrinkage_depots::Float64,
    shrinkage_charging::Float64,
    T::Int,
    seed::Int,
    B::Int,
    μ::Int,
    travel_cost_coeff::Int,
    charge_cost_coeff::Int,
    load_scale::Float64,
    load_shape::Float64,
    load_tolerance::Float64,
    batch::Int,
    permissiveness::Float64,
    data_dir::String = "data/",
)
    function complex_coords(
        n::Int,
    )
        return hcat(
            [round.(collect(reim(exp(2 * pi * j * im / n))), digits = 5)
            for j in 1:n]...
        )
    end

    function get_rectangle(
        n::Int,
    )
        a = Int(ceil(sqrt(n)))
        b = n ÷ a
        while b * a != n
            a -= 1
            b = n ÷ a
        end
        return a, b
    end

    function grid_coords(
        a::Int,
        b::Int,
        xmin::Float64 = -1.0,
        xmax::Float64 = 1.0,
        ymin::Float64 = -1.0,
        ymax::Float64 = 1.0,
    )
        return hcat(
            [
                [
                    xmin + ((xmax - xmin) * i) / (a - 1), 
                    ymin + ((ymax - ymin) * j) / (b - 1),
                ]
                for i in 0:a-1, j in 0:b-1
            ]...
        )
    end

    function circle_packing_coords(
        n::Int,
        ;
        data_dir::String = "data/",
    )
        lines = readlines(joinpath(data_dir, "cci_coords/cci$n.txt"))
        return hcat(
            [
                [parse(Float64, x[2]), parse(Float64, x[3])]
                for x in split.(lines)
            ]...
        )
    end

    function complex_coords_random(
        n::Int, 
        seed::Int,
    )
        Random.seed!(seed)
        deg = rand(n) * 2 * pi
        scale = rand(n)
        return hcat(
            [round.(collect(reim(scale[j] * exp(deg[j] * im))), digits = 5)
            for j in 1:n]...
        )
    end

    function uniform_coords_random(        
        n::Int, 
        seed::Int,
    )
        Random.seed!(seed)
        return -1 .+ rand(Float64, 2, n) .* 2
    end

    n_nodes = n_depots + n_customers + n_charging

    seeds = abs.(rand(MersenneTwister(seed), Int, 6))

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

    if customer_pattern == "random_box"
        customer_coords = uniform_coords_random(n_customers, seeds[1])
    elseif customer_pattern == "random_uniform_polar"
        customer_coords = complex_coords_random(n_customers, seeds[1])
    end
    if depot_pattern == "circular"
        depot_coords = shrinkage_depots * complex_coords(n_depots)
    end
    if charging_pattern == "circular"
        charging_coords = shrinkage_charging * complex_coords(n_charging)
    elseif charging_pattern == "grid"
        (a, b) = get_rectangle(n_charging)
        charging_coords = shrinkage_charging * grid_coords(a, b)
    elseif charging_pattern == "circular_packing"
        charging_coords = shrinkage_charging * circle_packing_coords(n_charging, data_dir = data_dir)
    end
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
    
    # c = Int.(round.(100 .* distances))
    # t = Int.(round.(100 .* distances)) # travel times are integer
    # q = Int.(round.(100 .* distances)) # charge costs are integer
    c = Int.(round.(distances .* 10000))
    t = Int.(round.(distances .* 10000) .* μ)
    q = Int.(round.(distances .* 10000))
    d = vcat(
        Int.(floor.(rand(Gamma(load_scale, load_shape), n_customers) ./ n_customers)),
        repeat([0], n_depots + n_charging),
    )
    C = Int(ceil(sum(d) * load_tolerance / n_vehicles))

    α, β = generate_times(T, n_customers, seeds[5], batch, permissiveness)
    α_charge = vcat(α, repeat([0], n_depots + n_charging))
    β_charge = vcat(β, repeat([T], n_depots + n_charging))

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
        shrinkage_depots,
        shrinkage_charging,
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
        T * μ,
        α_charge * μ,
        β_charge * μ,
        μ,
        B,
        travel_cost_coeff,
        charge_cost_coeff,
    )
    return data
end

function generate_graph_from_data(
    data::EVRPData,
)
    A = union(
        Set{Tuple{Int, Int}}(
            (i,i) for i in data.N_depots # allow self-loops at depots
        ),
        Set{Tuple{Int, Int}}(
            (i,j) for (i,j) in permutations(data.N_nodes, 2)
        ),
    )

    G = SimpleDiGraph{Int}(data.n_nodes)
    for (i, j) in A
        add_edge!(G, i, j)
    end

    t_ds = dijkstra_shortest_paths(G, collect(data.N_depots), data.t)
    q_ds = dijkstra_shortest_paths(G, collect(data.N_depots_charging), data.q)

    node_labels = merge(Dict(
        i => "Depot $ind" for (ind, i) in enumerate(data.N_depots)
    ), Dict(
        i => "Customer $ind" for (ind, i) in enumerate(data.N_customers)
    ), Dict(
        i => "Charging $ind" for (ind, i) in enumerate(data.N_charging)
    ))

    return EVRPGraph(
        G,
        node_labels,
        copy(data.c),
        copy(data.t),
        copy(data.q),
        data.N_customers,
        data.n_customers,
        data.N_depots,
        data.n_depots,
        data.N_charging,
        data.n_charging,
        data.N_depots_charging,
        data.n_depots_charging,
        data.N_nodes,
        data.n_nodes,
        A,
        data.T,
        data.B,
        data.μ,
        copy(data.α), 
        copy(data.β),
        t_ds.dists,
        q_ds.dists,
    )

end

function prune_graph(
    graph::EVRPGraph,
    k::Int,
)
    A = Set{Tuple{Int, Int}}()
    # 1. include self-loops in depots
    union!(A, [(i, i) for i in graph.N_depots])
    # 2. include all pairs between depots and charging stations
    union!(A, [(i, j) for (i, j) in permutations(graph.N_depots_charging, 2)])
    # 3. for each depot or charging station, include edges to and from nearest k customers
    for node in graph.N_depots_charging
        closest_customers = sortperm(graph.t[node, graph.N_customers])[1:k]
        union!(A, [(node, i) for i in closest_customers])
        union!(A, [(i, node) for i in closest_customers])
    end
    # 4. for each customer, include edges to and from nearest k customers (excluding themselves)
    for node in graph.N_customers
        closest_customers = setdiff(sortperm(graph.t[node, graph.N_customers])[1:k], node)
        union!(A, [(node, i) for i in closest_customers])
        union!(A, [(i, node) for i in closest_customers])
    end

    G = SimpleDiGraph{Int}(graph.n_nodes)
    for (i, j) in A
        add_edge!(G, i, j)
    end
    t_ds = dijkstra_shortest_paths(G, collect(graph.N_depots), graph.t)
    q_ds = dijkstra_shortest_paths(G, collect(graph.N_depots_charging), data.q)

    node_labels = merge(Dict(
        i => "Depot $ind" for (ind, i) in enumerate(graph.N_depots)
    ), Dict(
        i => "Customer $ind" for (ind, i) in enumerate(graph.N_customers)
    ), Dict(
        i => "Charging $ind" for (ind, i) in enumerate(graph.N_charging)
    ))

    return EVRPGraph(
        G,
        node_labels,
        copy(graph.c),
        copy(graph.t),
        copy(graph.q),
        graph.N_customers,
        graph.n_customers,
        graph.N_depots,
        graph.n_depots,
        graph.N_charging,
        graph.n_charging,
        graph.N_depots_charging,
        graph.n_depots_charging,
        graph.N_nodes,
        graph.n_nodes,
        A,
        graph.T,
        graph.B,
        graph.μ,
        copy(graph.α), 
        copy(graph.β),
        t_ds.dists,
        q_ds.dists,
    )
end

function generate_time_windows(
    T::Int,
    n_customers::Int,
    seed::Int,
    time_window_min_width::Float64,
    time_window_max_width::Float64,
)
    if !(
        0 < time_window_min_width ≤ time_window_max_width < 1
    )
        error("`time_window_min_width` and `time_window_max_width` out of bounds!")
    end
    Random.seed!(seed)
    time_window_dist = Uniform(time_window_min_width * T, time_window_max_width * T)
    time_window_widths = Int.(round.(rand(time_window_dist, n_customers)))
    time_window_posdist = Uniform(0.0, 1.0)
    time_window_pos = rand(time_window_posdist, n_customers)
    α = Int.(round.(time_window_pos .* (T .- time_window_widths)))
    β = α .+ time_window_widths

    return (α, β)
end

function generate_times(
    T::Int,
    n_customers::Int,
    seed::Int,
    batch::Int,
    permissiveness::Float64 = 0.4,
)
    if n_customers % batch != 0
        error()
    end
    times_dist = Uniform(0.0, T)
    α = zeros(Int, n_customers)
    β = zeros(Int, n_customers)

    Random.seed!(seed)
    for batch_ind in 1:(n_customers ÷ batch)
        inds = collect((batch_ind-1)*batch+1:batch_ind*batch)
        while true
            times = Int.(round.(sort(rand(times_dist, 2 * batch))))
            start_times = times[1:end÷2]
            end_times = times[end÷2+1:end]
            if all(
                (end_times .- start_times) ./ T .> permissiveness
            )
                α[inds] = start_times
                β[inds] = end_times
                break
            end
        end
    end
    return α, β
end

function charge_to_specified_level(
    starting_charge::Int, 
    desired_end_charge::Int, 
    starting_time::Int, 
)
    if desired_end_charge ≤ starting_charge
        return (0, starting_time, starting_charge)
    end
    delta = desired_end_charge - starting_charge
    end_time = starting_time + delta
    return (delta, end_time, desired_end_charge)
end

function compute_ngroute_neighborhoods(
    graph::EVRPGraph,
    k::Int, 
    ;
    depots_size::String = "small",
    charging_size::String = "small",
)
    if !(1 ≤ k ≤ graph.n_customers)
        error()
    end
    neighborhoods = falses(graph.n_nodes, graph.n_nodes)
    for i in graph.N_customers
        neighborhoods[sortperm(graph.c[i, graph.N_customers])[1:k], i] .= true
        # do not include any charging stations / depots in the neighborhoods of customers,
        # since there is no limit on repeat visits to charging stations / depots
    end
    if depots_size == "small"
        for i in graph.N_charging
            neighborhoods[i,i] = true
        end
    elseif depots_size == "medium"
        for i in graph.N_charging
            neighborhoods[i,i] = true
            neighborhoods[sortperm(graph.c[i, graph.N_customers])[1:k], i] .= true
        end
    elseif depots_size == "large"
        for i in graph.N_charging
            neighborhoods[i,i] = true
            neighborhoods[graph.N_customers, i] .= true
        end
    else
        error("`depots_size` argument not recognized.")
    end
    if charging_size == "small"
        for i in graph.N_charging
            neighborhoods[i,i] = true
        end
    elseif charging_size == "medium"
        for i in graph.N_charging
            neighborhoods[i,i] = true
            neighborhoods[sortperm(graph.c[i, graph.N_customers])[1:k], i] .= true
        end
    elseif charging_size == "large"
        for i in graph.N_charging
            neighborhoods[i,i] = true
            neighborhoods[graph.N_customers, i] .= true
        end
    else
        error("`charging_size` argument not recognized.")
    end
    return neighborhoods
end

function ngroute_check_create_fset(
    neighborhoods::BitMatrix,
    set::BitVector,
    next_node::Int,
)
    if set[next_node]
        # if next_node is a customer not yet visited, proceed
        # only if one can extend current_subpath along next_node according to ng-route rules
        return (false, set)
    end
    new_set = copy(set)
    new_set .&= neighborhoods[:, next_node]
    new_set[next_node] = true
    return (true, new_set)
end

function ngroute_create_bset(
    neighborhoods::BitMatrix,
    nodes::Vector{Int},
    set::BitVector,
)
    new_set = copy(set)
    val = neighborhoods[nodes[end], nodes[1]]
    for i in nodes[2:end-1]
        val &= neighborhoods[nodes[end], i]
    end
    new_set[nodes[end]] |= val
    return new_set 
end

function compute_arc_modified_costs(
    graph::EVRPGraph,
    data::EVRPData,
    ν::Vector{Float64}, 
    ;
)
    modified_costs = data.travel_cost_coeff * Float64.(copy(graph.c))
    for j in graph.N_customers
        for i in graph.N_nodes
            modified_costs[i,j] -= ν[j]
        end
    end
    return modified_costs
end


function prepare_lambda(
    λ::Dict{NTuple{3, Int}, Float64},
    n_nodes::Int,
)
    λvals = collect(values(λ))
    λcust = falses(length(λ), n_nodes)
    for (i, k) in enumerate(keys(λ))
        λcust[i, collect(k)] .= true
    end
    λmemory = falses(length(λ))
    return λvals, λcust
end

function prepare_lambda(
    λ::Dict{Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}, Float64},
    n_nodes::Int,
)
    λvals = collect(values(λ))
    λcust = falses(length(λ), n_nodes)
    λmemory = falses(length(λ), n_nodes)
    for (i, k) in enumerate(keys(λ))
        (S, M) = k
        λcust[i, collect(S)] .= true
        λmemory[i, collect(M)] .= true
    end
    return λvals, λcust, λmemory
end


function plot_instance(
    data::EVRPData,
)
    p = plot(
        # xlim = (0, 1), ylim = (0, 1),
        aspect_ratio = :equal, 
        fmt = :png, 
    )
    plot!(
        data.customer_coords[1,:], data.customer_coords[2,:],
        seriestype = :scatter, 
        label = "Customer",
        color = :green
    )
    annotate!.(
        data.customer_coords[1,:] .+ 0.1, data.customer_coords[2,:], 
        text.(
            collect(string(i) for i in 1:data.n_customers), 
            :green, :left, 11
        )
    )
    plot!(
        data.depot_coords[1,:], data.depot_coords[2,:],
        seriestype = :scatter, 
        label = "Depots",
        color = :black
    )
    annotate!.(
        data.depot_coords[1,:] .+ 0.1, data.depot_coords[2,:], 
        text.(
            collect("M" * string(i) for i in 1:data.n_depots), 
            :black, :left, 11
        )
    )
    plot!(
        data.charging_coords[1,:], data.charging_coords[2,:],
        seriestype = :scatter, 
        label = "Charging stations",
        color = :grey
    )
    annotate!.(
        data.charging_coords[1,:] .+ 0.1, data.charging_coords[2,:], 
        text.(
            collect("R" * string(i) for i in 1:data.n_charging), 
            :grey, :left, 11
        )
    )

    plot!(legend = :outerright)
    return p
end

function compute_path_metrics(
    some_paths::Dict{T, Vector{Path}},
) where T
    total_subpath_length = 0.0
    num_subpaths = 0.0
    total_path_length = 0.0
    num_paths = 0.0
    total_ps_length = 0.0
    for path_l in values(some_paths)
        for p in path_l
            if (
                length(p.subpaths) == 1 
                && (
                    p.subpaths[1].artificial # artificial path
                    # || length(p.subpaths[1].arcs) == 1 # path from depot to depot
                )
            )
                continue
            end
            total_subpath_length += sum(sum(s.served) + 1 for s in p.subpaths) 
            num_subpaths += length(p.subpaths)
            total_path_length += sum(p.served) + length(p.subpaths)
            num_paths += 1
            total_ps_length += length(p.subpaths)
        end
    end

    return Dict(
        "mean_subpath_length" => total_subpath_length / num_subpaths,
        "mean_path_length" => total_path_length / num_paths,
        "mean_ps_length" => total_ps_length / num_paths,
    )
end