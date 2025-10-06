include("utils.jl")

Base.@kwdef mutable struct HPathLabel
    arcs::Vector{Vector{Tuple{Int, Int}}}
    times::Vector{Int}
    total_time::Int
    charges::Vector{Vector{Int}}
    total_charge::Int
    load::Int
    served::BitVector
end

function compute_insertion_charging_stations(
    r::HPathLabel,
    cust::Int,
    s_ind::Int, # Subpath index
    a_ind::Int, # Arc index within subpath to insert cust
    graph::EVRPGraph,
    charging_stations_shortest_paths::Dict{Tuple{Int, Int}, Tuple{Int, Vector{Int}}},
)
    subpath = r.arcs[s_ind]
    arc = subpath[a_ind]
    (i, j) = arc
    # Quick time horizon check
    if (
        r.total_time + graph.t[i, cust] + graph.t[cust, j] - graph.t[i, j]
        + (
            r.total_charge 
            + graph.q[i, cust] + graph.q[cust, j] - graph.q[i, j]
            - graph.B
        )
    ) > graph.T
        return nothing
    end
    
    q_left = sum(r.charges[s_ind][1:a_ind-1])
    q_right = sum(r.charges[s_ind][a_ind+1:end])
    if (q_left + graph.q[i, cust] + graph.q[cust, j] + q_right) ≤ graph.B
        return (
            graph.q[i, cust] 
            + graph.q[cust, j], 
            Int[], Int[],
        )
    end

    # Require charging station insertion
    options = Tuple{Int, Vector{Int}, Vector{Int}}[]
    add_left = false
    add_right = false
    if q_left + graph.q[i, cust] > graph.B
        add_left = true
    end
    if graph.q[cust, j] + q_right > graph.B
        add_right = true
    end

    # No insertions left of cust
    added = false
    if !add_left
        for cs_l in intersect(graph.N_charging, outneighbors(graph.G, cust))
            for cs_r in intersect(graph.N_charging, outneighbors(graph.G, j))
                if (
                    (q_left + graph.q[i, cust] + graph.q[cust, cs_l]) ≤ graph.B
                ) && (
                    (graph.q[cs_r, j] + q_right) ≤ graph.B
                )
                    push!(
                        options, 
                        (
                            (
                                graph.q[i, cust] 
                                + graph.q[cust, cs_l] 
                                + charging_stations_shortest_paths[(cs_l, cs_r)][1]
                                + graph.q[cs_r, j]
                            ),
                            Int[], charging_stations_shortest_paths[(cs_l, cs_r)][2],
                        )
                    )
                end
            end
        end
    end

    # No insertions right of cust
    added = false
    if !add_right
        for cs_l in intersect(graph.N_charging, outneighbors(graph.G, i))
            for cs_r in intersect(graph.N_charging, outneighbors(graph.G, cust))
                if (
                    (q_left + graph.q[i, cs_l]) ≤ graph.B
                ) && (
                    (graph.q[cs_r, cust] + graph.q[cust, j] + q_right) ≤ graph.B
                )
                    push!(
                        options, 
                        (
                            (
                                graph.q[i, cs_l]
                                + charging_stations_shortest_paths[(cs_l, cs_r)][1]
                                + graph.q[cs_r, cust] 
                                + graph.q[cust, j]
                            ),
                            charging_stations_shortest_paths[(cs_l, cs_r)][2], Int[],
                        )
                    )
                end
            end
        end
    end

    # Insert both left and right of cust
    for cs_Ll in intersect(graph.N_charging, outneighbors(graph.G, i))
        for cs_Lr in intersect(graph.N_charging, outneighbors(graph.G, cust))
            for cs_Rl in intersect(graph.N_charging, outneighbors(graph.G, cust))
                for cs_Rr in intersect(graph.N_charging, outneighbors(graph.G, j))
                    if (
                        (q_left + graph.q[i, cs_Ll]) ≤ graph.B
                    ) && (
                        (graph.q[cs_Lr, cust] + graph.q[cust, cs_Rl]) ≤ graph.B
                    ) && (
                        (graph.q[cs_Rr, j] + q_right) ≤ graph.B
                    )
                        push!(
                            options, 
                            (
                                (
                                    graph.q[i, cs_Ll]
                                    + charging_stations_shortest_paths[(cs_Ll, cs_Lr)][1]
                                    + graph.q[cs_Lr, cust] 
                                    + graph.q[cust, cs_Rl]
                                    + charging_stations_shortest_paths[(cs_Rl, cs_Rr)][1]
                                    + graph.q[cs_Rr, j]
                                ),
                                charging_stations_shortest_paths[(cs_Ll, cs_Lr)][2],
                                charging_stations_shortest_paths[(cs_Rl, cs_Rr)][2],
                            )
                        )
                    end
                end
            end
        end
    end

    if length(options) == 0
        return nothing
    end
    sort!(options)
    return options[1]
end




function generate_option(
    r::HPathLabel,
    r_ind::Int,
    cust::Int,
    graph::EVRPGraph,
    charging_stations_shortest_paths::Dict{Tuple{Int, Int}, Tuple{Int, Vector{Int}}},
)
    options = Tuple{Int, Int, Int, Int, Int, Vector{Int}, Vector{Int}}[]
    for (s_ind, subpath) in enumerate(r.arcs)
        for (a_ind, (i, j)) in enumerate(subpath)
            option = compute_insertion_charging_stations(r, cust, s_ind, a_ind, graph, charging_stations_shortest_paths)
            if isnothing(option)
                continue
            end
            push!(options, (option[1] - graph.q[i,j], r_ind, cust, s_ind, a_ind, option[2], option[3]))
        end
    end
    sort!(options)
    return options
end




# unserved_customers

function convert_heuristic_path_label_to_path(
    r::HPathLabel,
    graph::EVRPGraph,
)
    current_time, current_charge = (0, graph.B)
    prev_time, prev_charge = current_time, current_charge
    p = Path(
        subpaths = Subpath[],
        charging_arcs = ChargingArc[],
        served = r.served,
        load = r.load,
        arcs = NTuple{2, Int}[],
        customer_arcs = NTuple{2, Int}[],
    )
    for (s_ind, subpath) in enumerate(r.arcs)
        prev_time = current_time
        prev_charge = current_charge
        @assert sum(graph.t[i,j] for (i,j) in subpath) == r.times[s_ind]
        @assert [graph.q[i,j] for (i,j) in subpath] == r.charges[s_ind]
        current_node = subpath[end][end]
        current_time = current_time + r.times[s_ind]
        current_charge = current_charge - sum(r.charges[s_ind])
        served = [count(x -> x[1] == i, subpath) for i in graph.N_customers]
        s = Subpath(
            n_customers = graph.n_customers,
            starting_node = subpath[1][1],
            starting_time = prev_time,
            starting_charge = prev_charge,
            current_node = current_node,
            arcs = subpath,
            current_time = current_time,
            current_charge = current_charge,
            load = sum(graph.d[i] for (i, _) in subpath),
            served = served,
        )
        push!(p.subpaths, s)
        if s_ind == length(r.arcs)
            break
        end
        prev_time = current_time
        prev_charge = current_charge
        (delta, current_time, current_charge) = charge_to_specified_level(
            prev_charge,
            sum(r.charges[s_ind+1]),
            prev_time,
        )

        a = ChargingArc(
            starting_node = current_node, 
            starting_time = prev_time, 
            starting_charge = prev_charge, 
            delta = delta,
            charge_cost_coeff = data.charge_cost_coeffs[current_node],
            current_time = current_time, 
            current_charge = current_charge,
        )
        push!(p.charging_arcs, a)
    end
    p.served = sum(s.served for s in p.subpaths)
    @assert p.served == r.served
    p.load = sum(s.load for s in p.subpaths)
    @assert p.load == r.load
    p.arcs = vcat([s.arcs for s in p.subpaths]...)
    @assert p.arcs == vcat(r.arcs...)
    customers = [a[1] for a in p.arcs if a[1] in graph.N_customers]
    p.customer_arcs = collect(zip(customers[1:end-1], customers[2:end]))
    return p

end

function compute_heuristic_paths_notimewindows(
    data::EVRPData,
    graph::EVRPGraph,
    ;
    use_load::Bool = true,
    verbose::Bool = false,
)

    ### Initialization
    # Initialize start and end depots
    start_depots = vcat([repeat([k], v) for (k, v) in data.v_start]...)
    end_depots = repeat(collect(data.N_depots), outer = Int(ceil(data.n_vehicles / data.n_depots)))[1:data.n_vehicles]
    unserved_customers = copy(graph.N_customers)

    # Precompute shortest paths between charging stations
    G_depots_charging, vmap = induced_subgraph(graph.G, graph.N_depots_charging)
    rmap = Dict(j => i for (i, j) in pairs(vmap))

    dstate = Graphs.floyd_warshall_shortest_paths(G_depots_charging, graph.q[graph.N_depots_charging, graph.N_depots_charging])

    charging_stations_shortest_paths = Dict{Tuple{Int, Int}, Tuple{Int, Vector{Int}}}()
    for cs1 in vertices(G_depots_charging)
        for cs2 in vertices(G_depots_charging)
            j = cs2
            nodeseq = Int[cs2]
            cost = 0
            while j != cs1
                j_prev = dstate.parents[cs1, j]
                cost += graph.q[vmap[j_prev], vmap[j]]
                j = j_prev
                push!(nodeseq, j)
            end
            charging_stations_shortest_paths[(vmap[cs1], vmap[cs2])] = (cost, [vmap[i] for i in nodeseq[end:-1:1]])
        end
    end

    # Initialize partial routes
    partial_routes = HPathLabel[]
    for (start_depot, end_depot) in zip(start_depots, end_depots)
        if start_depot == end_depot
            push!(partial_routes, HPathLabel(
                arcs = [[(start_depot, end_depot)]],
                times = [0],
                total_time = 0,
                charges = [[0]],
                total_charge = 0,
                load = 0,
                served = falses(graph.n_customers),
            ))
        else
            c = charging_stations_shortest_paths[(start_depot, end_depot)]
            arcs = [
                [(i, j)]
                for (i, j) in zip(c[2][1:end-1], c[2][2:end])
            ]
            times = [sum(graph.t[i,j] for (i,j) in subpath) for subpath in arcs]
            charges = [[graph.q[i,j] for (i,j) in subpath] for subpath in arcs]
            push!(partial_routes, HPathLabel(
                arcs = arcs,
                times = times,
                total_time = sum(times),
                charges = charges,
                total_charge = sum(sum(c) for c in charges),
                load = 0,
                served = falses(graph.n_customers),
            ))
        end
    end

    # Repeat until all customers are inserted
    terminate_early = false
    while length(unserved_customers) > 0
        options = Tuple{Int, Int, Int, Int, Int, Vector{Int}, Vector{Int}}[]
        for (r_ind, r) in enumerate(partial_routes)
            for cust in unserved_customers
                if use_load && (r.load + graph.d[cust] > graph.C)
                    continue
                end
                append!(options, generate_option(r, r_ind, cust, graph, charging_stations_shortest_paths))
            end
        end
        sort!(options)

        while true
            if length(options) == 0
                terminate_early = true
                break
            end
            option = popfirst!(options)
            (_, r_ind, cust, s_ind, a_ind, CS_left, CS_right) = option
            r = partial_routes[r_ind]
            (i, j) = r.arcs[s_ind][a_ind]


            if length(CS_left) == 0 && length(CS_right) == 0
                new_subpaths = [
                    vcat(r.arcs[s_ind][1:a_ind-1], [(i, cust), (cust, j)], r.arcs[s_ind][a_ind+1:end])
                ]
            elseif length(CS_left) > 0 && length(CS_right) == 0
                new_subpaths = [
                    vcat(r.arcs[s_ind][1:a_ind-1], [(i, CS_left[1])]),
                    [[(l,r)] for (l,r) in zip(CS_left[1:end-1], CS_left[2:end])]...,
                    vcat([(CS_left[end], cust), (cust, j)], r.arcs[s_ind][a_ind+1:end]),
                ]
            elseif length(CS_left) == 0 && length(CS_right) > 0
                new_subpaths = [
                    vcat(r.arcs[s_ind][1:a_ind-1], [(i, cust), (cust, CS_right[1])]),
                    [[(l,r)] for (l,r) in zip(CS_right[1:end-1], CS_right[2:end])]...,
                    vcat([(CS_right[end], j)], r.arcs[s_ind][a_ind+1:end]),
                ]
            else
                new_subpaths = [
                    vcat(r.arcs[s_ind][1:a_ind-1], [(i, CS_left[1])]),
                    [[(l,r)] for (l,r) in zip(CS_left[1:end-1], CS_left[2:end])]...,
                    vcat([(CS_left[end], cust), (cust, CS_right[1])]),
                    [[(l,r)] for (l,r) in zip(CS_right[1:end-1], CS_right[2:end])]...,
                    vcat([(CS_right[end], j)], r.arcs[s_ind][a_ind+1:end]),
                ]
            end
            new_subpaths_times = [sum(graph.t[i,j] for (i,j) in s) for s in new_subpaths]
            new_subpaths_charges = [[graph.q[i,j] for (i,j) in s] for s in new_subpaths]
            new_total_time = sum(r.times) + sum(new_subpaths_times) - sum(r.times[s_ind])
            new_total_charge = r.total_charge + sum(sum(c) for c in new_subpaths_charges) - sum(r.charges[s_ind])
            if new_total_time + (new_total_charge - graph.B) > graph.T
                continue
            end
            
            verbose && println("Inserting customer $cust into route $r_ind, subpath $s_ind, arc ($i, $j)")
            r.arcs = vcat(r.arcs[1:s_ind-1], new_subpaths, r.arcs[s_ind+1:end])
            r.times = vcat(r.times[1:s_ind-1], new_subpaths_times, r.times[s_ind+1:end])
            r.total_time = new_total_time
            r.charges = vcat(r.charges[1:s_ind-1], new_subpaths_charges, r.charges[s_ind+1:end])
            r.total_charge = new_total_charge

            r.load += graph.d[cust]
            r.served[cust] = true
            unserved_customers = setdiff(unserved_customers, [cust])
            verbose && println("Unserved customers: ", unserved_customers)
            break
        end
        if terminate_early
            return Path[]
        end
    end

    heuristic_paths = Dict{NTuple{2, NTuple{3, Int}}, Vector{Path}}()
    for r in partial_routes
        p = convert_heuristic_path_label_to_path(r, graph)
        key = (
            (p.subpaths[1].starting_node, p.subpaths[1].starting_time, p.subpaths[1].starting_charge),
            (p.subpaths[end].current_node, p.subpaths[end].current_time, p.subpaths[end].current_charge),
        )
        if haskey(heuristic_paths, key)
            push!(heuristic_paths[key], p)
        else
            heuristic_paths[key] = [p]
        end
    end
    return heuristic_paths
end
