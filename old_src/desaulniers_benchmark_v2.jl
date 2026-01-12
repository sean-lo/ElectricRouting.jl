using Pkg
Pkg.activate(".")
include("utils.jl")

abstract type BPathLabel <: Label end

function get_bpath_vkey_type(
    use_load::Load,
    customer_service::CustomerService,
    cuts::Cuts,
)
    if use_load isa YesLoad
        if customer_service isa Elementary
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                Int, # load
                BitVector, # served
            }
        elseif customer_service isa NoService
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                Int, # load
            }
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return Tuple{
                    Float64, # cost
                    Int, Int, Int,
                    Int, # load
                    BitVector, # ng_fset
                }
            elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                return Tuple{
                    Float64, # cost
                    BitVector, # cut_flabels
                    Int, Int, Int,
                    Int, # load
                    BitVector, # ng_fset
                }
            end
        end
    elseif use_load isa NoLoad
        if customer_service isa Elementary
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                BitVector, # served
            }
        elseif customer_service isa NoService
            return Tuple{
                Float64, # cost
                Int, Int, Int,
            }
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return Tuple{
                    Float64, # cost
                    Int, Int, Int,
                    BitVector, # ng_fset
                }
            elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                return Tuple{
                    Float64, # cost
                    BitVector, # cut_flabels
                    Int, Int, Int,
                    BitVector, # ng_fset
                }
            end
        end
    end
end

function get_bpath_vkey_fkey_type(
    use_load::Load,
    customer_service::CustomerService,
    cuts::Cuts,
)
    """
    (
        cost::Float64,
        if (
            customer_service isa NgRoute
            && (cuts isa SR3Cuts || cuts isa LmSR3Cuts)
        )
            cut_flabels::BitVector (n_cuts,)
        end
        min time::Int,
        negative of max charge::Int,
        difference between min time and min charge::Int,
        if use_load isa YesLoad
            load::Int,
        end
        if customer_service isa Elementary
            served::BitVector (n_customers,)
        elseif customer_service isa NgRoute
            ng_fset::BitVector (n_nodes,)
        end
        starting_node::Int,
        current_node::Int,
    )
    """

    if use_load isa YesLoad
        if customer_service isa Elementary
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                Int, # load
                BitVector, # served
                Int, Int, # starting_node, current_node
            }
        elseif customer_service isa NoService
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                Int, # load
                Int, Int, # starting_node, current_node
            }
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return Tuple{
                    Float64, # cost
                    Int, Int, Int,
                    Int, # load
                    BitVector, # ng_fset
                    Int, Int, # starting_node, current_node
                }
            elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                return Tuple{
                    Float64, # cost
                    BitVector, # cut_flabels
                    Int, Int, Int,
                    Int, # load
                    BitVector, # ng_fset
                    Int, Int, # starting_node, current_node
                }
            end
        end
    elseif use_load isa NoLoad
        if customer_service isa Elementary
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                BitVector, # served
                Int, Int, # starting_node, current_node
            }
        elseif customer_service isa NoService
            return Tuple{
                Float64, # cost
                Int, Int, Int,
                Int, Int, # starting_node, current_node
            }
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return Tuple{
                    Float64, # cost
                    Int, Int, Int,
                    BitVector, # ng_fset
                    Int, Int, # starting_node, current_node
                }
            elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                return Tuple{
                    Float64, # cost
                    BitVector, # cut_flabels
                    Int, Int, Int,
                    BitVector, # ng_fset
                    Int, Int, # starting_node, current_node
                }
            end
        end
    end
end

abstract type YesLoadBPathLabel <: BPathLabel end

mutable struct YesLoadElementaryBPathLabel <: YesLoadBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    load::Int
    served::BitVector
end
get_vkey(p::YesLoadElementaryBPathLabel) = (
    p.cost,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.load,
    p.served,
)
mutable struct YesLoadNoServiceBPathLabel <: YesLoadBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    load::Int
end
get_vkey(p::YesLoadNoServiceBPathLabel) = (
    p.cost,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.load,
)

abstract type YesLoadNgRouteBPathLabel <: YesLoadBPathLabel end

mutable struct YesLoadNgRouteNoCutsBPathLabel <: YesLoadNgRouteBPathLabel 
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    load::Int
    ng_fset::BitVector
end
get_vkey(p::YesLoadNgRouteNoCutsBPathLabel) = (
    p.cost,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.load,
    p.ng_fset,
)
mutable struct YesLoadNgRouteSR3CutsBPathLabel <: YesLoadNgRouteBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    load::Int
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::YesLoadNgRouteSR3CutsBPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.load,
    p.ng_fset,
)
mutable struct YesLoadNgRouteLmSR3CutsBPathLabel <: YesLoadNgRouteBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    load::Int
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::YesLoadNgRouteLmSR3CutsBPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.load,
    p.ng_fset,
)

abstract type NoLoadBPathLabel <: BPathLabel end

mutable struct NoLoadElementaryBPathLabel <: NoLoadBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    served::BitVector
end
get_vkey(p::NoLoadElementaryBPathLabel) = (
    p.cost,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.served,
)
mutable struct NoLoadNoServiceBPathLabel <: NoLoadBPathLabel 
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
end
get_vkey(p::NoLoadNoServiceBPathLabel) = (
    p.cost,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
)

abstract type NoLoadNgRouteBPathLabel <: NoLoadBPathLabel end

mutable struct NoLoadNgRouteNoCutsBPathLabel <: NoLoadNgRouteBPathLabel 
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    ng_fset::BitVector
end
get_vkey(p::NoLoadNgRouteNoCutsBPathLabel) = (
    p.cost,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.ng_fset,
)
mutable struct NoLoadNgRouteSR3CutsBPathLabel <: NoLoadNgRouteBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::NoLoadNgRouteSR3CutsBPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.ng_fset,
)
mutable struct NoLoadNgRouteLmSR3CutsBPathLabel <: NoLoadNgRouteBPathLabel
    cost::Float64
    nodes::Vector{Int}
    excesses::Vector{Int}
    slacks::Vector{Int}
    time_mincharge::Int
    time_maxcharge::Int
    charge_mincharge::Int
    charge_maxcharge::Int
    ng_fset::BitVector
    cut_flabels::BitVector
end
get_vkey(p::NoLoadNgRouteLmSR3CutsBPathLabel) = (
    p.cost,
    p.cut_flabels,
    p.time_mincharge,
    - p.charge_maxcharge,
    p.time_mincharge - p.charge_mincharge,
    p.ng_fset,
)

function get_bpath_label_type(
    use_load::Load,
    customer_service::CustomerService,
    cuts::Cuts,
)
    if use_load isa YesLoad
        if customer_service isa Elementary
            return YesLoadElementaryBPathLabel
        elseif customer_service isa NoService
            return YesLoadNoServiceBPathLabel
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return YesLoadNgRouteNoCutsBPathLabel
            elseif cuts isa SR3Cuts
                return YesLoadNgRouteSR3CutsBPathLabel
            elseif cuts isa LmSR3Cuts
                return YesLoadNgRouteLmSR3CutsBPathLabel
            end
        end
    elseif use_load isa NoLoad
        if customer_service isa Elementary
            return NoLoadElementaryBPathLabel
        elseif customer_service isa NoService
            return NoLoadNoServiceBPathLabel
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return NoLoadNgRouteNoCutsBPathLabel
            elseif cuts isa SR3Cuts
                return NoLoadNgRouteSR3CutsBPathLabel
            elseif cuts isa LmSR3Cuts
                return NoLoadNgRouteLmSR3CutsBPathLabel
            end
        end
    end
end

function create_new_bpath_label(
    use_load::Load,
    customer_service::CustomerService,
    cuts::Cuts,
    starting_node::Int,
    data::EVRPData,
    ;
    n_cuts::Int = 0,
)
    if customer_service isa NgRoute
        ng_fset = falses(data.n_nodes)
        ng_fset[starting_node] = true
    end

    if use_load isa YesLoad
        if customer_service isa Elementary
            return YesLoadElementaryBPathLabel(
                0.0, [starting_node], Int[], Int[],
                0, 0, data.B, data.B, 
                0, 
                falses(data.n_customers),
            )
        elseif customer_service isa NoService
            return YesLoadNoServiceBPathLabel(
                0.0, [starting_node], Int[], Int[],
                0, 0, data.B, data.B, 
                0, 
            )
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return YesLoadNgRouteNoCutsBPathLabel(
                    0.0, [starting_node], Int[], Int[],
                    0, 0, data.B, data.B, 
                    0, 
                    ng_fset,
                )
            elseif cuts isa SR3Cuts
                return YesLoadNgRouteSR3CutsBPathLabel(
                    0.0, [starting_node], Int[], Int[],
                    0, 0, data.B, data.B, 
                    0, 
                    ng_fset,
                    falses(n_cuts),
                )
            elseif cuts isa LmSR3Cuts
                return YesLoadNgRouteLmSR3CutsBPathLabel(
                    0.0, [starting_node], Int[], Int[],
                    0, 0, data.B, data.B, 
                    0, 
                    ng_fset,
                    falses(n_cuts),
                )
            end
        end
    elseif use_load isa NoLoad
        if customer_service isa Elementary
            return NoLoadElementaryBPathLabel(
                0.0, [starting_node], Int[], Int[],
                0, 0, data.B, data.B, 
                falses(data.n_customers),
            )
        elseif customer_service isa NoService
            return NoLoadNoServiceBPathLabel(
                0.0, [starting_node], Int[], Int[],
                0, 0, data.B, data.B, 
            )
        elseif customer_service isa NgRoute
            if cuts isa NoCuts
                return NoLoadNgRouteNoCutsBPathLabel(
                    0.0, [starting_node], Int[], Int[],
                    0, 0, data.B, data.B, 
                    ng_fset,
                )
            elseif cuts isa SR3Cuts
                return NoLoadNgRouteSR3CutsBPathLabel(
                    0.0, [starting_node], Int[], Int[],
                    0, 0, data.B, data.B, 
                    ng_fset,
                    falses(n_cuts),
                )
            elseif cuts isa LmSR3Cuts
                return NoLoadNgRouteLmSR3CutsBPathLabel(
                    0.0, [starting_node], Int[], Int[],
                    0, 0, data.B, data.B, 
                    ng_fset,
                    falses(n_cuts),
                )
            end
        end
    end
end


function update_path_ngroute!(
    new_path::BPathLabel,
    next_node::Int,
    neighborhoods::BitMatrix
)
    # update forward ng-set
    ngroute_update_fset!(new_path, next_node, neighborhoods)
    return 
end

function update_path_cut_labels_SR3!(
    new_path::BPathLabel,
    next_node::Int,
    λvals::Vector{Float64},
    λcust::BitMatrix,
)
    ## IMPORTANT: update cost before flabels
    ## since flabels affect cost
    SR3_update_cost!(new_path, next_node, λvals, λcust)
    SR3_update_cut_flabels!(new_path, next_node, λcust)
    return
end

function update_path_cut_labels_LmSR3!(
    new_path::BPathLabel,
    next_node::Int,
    λvals::Vector{Float64},
    λcust::BitMatrix,
    λmemory::BitMatrix,
)
    ## IMPORTANT: update cost before flabels
    ## since flabels affect cost
    lmSR3_update_cost!(new_path, next_node, λvals, λcust, λmemory)
    lmSR3_update_cut_flabels!(new_path, next_node, λcust, λmemory)
    return
end


function compute_new_bpath(
    current_path::BPathLabel,
    data::EVRPData,
    graph::EVRPGraph,
    modified_costs::Matrix{Float64},
    current_node::Int,
    next_node::Int,
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
    ;
    α::Vector{Int} = zeros(graph.n_nodes),
    β::Vector{Int} = fill(graph.T, graph.n_nodes),
    neighborhoods::BitMatrix = falses(graph.n_nodes, graph.n_nodes),
    n_cuts::Int = 0,
    λvals::Vector{Float64} = Float64[],
    λcust::BitMatrix = falses(n_cuts, graph.n_nodes),
    λmemory::BitMatrix = falses(n_cuts, graph.n_nodes),
)
    # customer service
    if customer_service isa Elementary
        if next_node in graph.N_customers && current_path.served[next_node]
            return (false, current_path)
        end
    elseif customer_service isa NgRoute
        if current_path.ng_fset[next_node]
            return (false, current_path)
        end
    end

    # time and charge feasibility
    # (1) battery
    excess = max(
        0, 
        graph.q[current_node,next_node] - current_path.charge_mincharge 
    )
    # (2) time windows
    if use_time_windows isa YesTimeWindows
        if current_path.time_mincharge + excess + graph.t[current_node,next_node] > β[next_node]
            # println("$(current_path.time_mincharge), $excess, $(t[current_node,next_node]), $(β[next_node])")
            # println("not time windows feasible")
            return (false, current_path)
        end
    end
    if current_path.time_mincharge + excess + graph.t[current_node,next_node] + graph.min_t[next_node] > graph.T
        return (false, current_path)
    end
    # (3) charge interval 
    if current_node in graph.N_charging
        if excess > max(graph.B - current_path.charge_mincharge, 0)
            # println("$excess, $(B), $(current_path.charge_mincharge)")
            # println("not charge feasible")
            return (false, current_path)
        end
    else
        if excess > max(current_path.charge_maxcharge - current_path.charge_mincharge, 0)
            # println("$excess, $(current_path.charge_maxcharge), $(current_path.charge_mincharge)")
            # println("not charge feasible")
            return (false, current_path)
        end
    end

    # load feasibility
    if use_load isa YesLoad
        if current_path.load + graph.d[next_node] > graph.C
            return (false, current_path)
        end
    end

    new_path = copy(current_path)
    push!(new_path.nodes, next_node)
    push!(new_path.excesses, excess)
    
    if use_time_windows isa YesTimeWindows
        new_path.time_mincharge = max(
            α[next_node],
            current_path.time_mincharge + graph.t[current_node,next_node] + excess
        )
        if current_node in graph.N_charging
            slack = min(
                new_path.time_mincharge - (current_path.time_mincharge + graph.t[current_node,next_node] + excess),
                graph.B - (current_path.charge_mincharge + excess),
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[next_node],
                max(
                    α[next_node],
                    current_path.time_maxcharge + (graph.B - current_path.charge_maxcharge) + graph.t[current_node,next_node],
                )
            )
        else
            slack = min(
                new_path.time_mincharge - (current_path.time_mincharge + graph.t[current_node,next_node] + excess),
                current_path.charge_maxcharge - (current_path.charge_mincharge + excess),
            )
            push!(new_path.slacks, slack)
            new_path.time_maxcharge = min(
                β[next_node],
                max(
                    α[next_node],
                    current_path.time_maxcharge + graph.t[current_node,next_node],
                )
            )
        end
    else
        new_path.time_mincharge = current_path.time_mincharge + graph.t[current_node,next_node] + excess
        slack = 0

        if current_node in graph.N_charging
            # slack = min(0, graph.B - max(current_path.charge_mincharge, graph.q[current_node,next_node]))
            new_path.time_maxcharge = min(
                graph.T,
                current_path.time_maxcharge + (graph.B - current_path.charge_maxcharge) + graph.t[current_node,next_node],
            )
        else
            # slack = min(
            #     0,
            #     current_path.charge_maxcharge - (current_path.charge_mincharge + excess),
            # )
            new_path.time_maxcharge = min(
                graph.T,
                current_path.time_maxcharge + graph.t[current_node,next_node],
            )
        end
        push!(new_path.slacks, slack)
    end

    new_path.charge_mincharge = (
        current_path.charge_mincharge 
        + excess 
        + slack
        - graph.q[current_node,next_node]
    )
    new_path.charge_maxcharge = (
        new_path.charge_mincharge 
        + new_path.time_maxcharge 
        - new_path.time_mincharge
    )

    new_path.cost += modified_costs[current_node,next_node]
    # Assume only homogenous charging costs for benchmark method
    # Practically - these charging amounts are performed at previously visited CSes
    new_path.cost += data.charge_cost_coeff * (excess + slack)


    # Load
    if use_load isa YesLoad
        new_path.load += graph.d[next_node]
    end

    # Customer service
    if customer_service isa Elementary
        if next_node in graph.N_customers
            new_path.served[next_node] = true
        end
    elseif customer_service isa NgRoute
        update_path_ngroute!(new_path, next_node, neighborhoods)
    end

    # Cuts
    if cuts isa SR3Cuts
        update_path_cut_labels_SR3!(new_path, next_node, λvals, λcust)
    elseif cuts isa LmSR3Cuts
        update_path_cut_labels_LmSR3!(new_path, next_node, λvals, λcust, λmemory)
    end

    return (true, new_path)
end

function generate_path_labels_benchmark(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64}, 
    use_load::Load,
    use_time_windows::TimeWindows,
    customer_service::CustomerService,
    cuts::Cuts,
    ;
    neighborhoods::BitMatrix = falses(graph.n_nodes, graph.n_nodes),
    λ::Dict{<:Tuple, Float64} = Dict(),
    time_limit::Float64 = Inf, 
)
    start_time = time()
    modified_costs = compute_arc_modified_costs(graph, data, ν)
    if cuts isa NoCuts
        λvals = Float64[]
        λcust = falses(length(λ), graph.n_nodes)
        λmemory = falses(length(λ), graph.n_nodes)
    elseif cuts isa SR3Cuts
        λvals, λcust = prepare_lambda(λ, graph.n_nodes)
        λmemory = falses(length(λ), graph.n_nodes)
    elseif cuts isa LmSR3Cuts
        λvals, λcust, λmemory = prepare_lambda(λ, graph.n_nodes)
    end

    path_labels = Dict(
        (starting_node, current_node) => SortedDict{
            get_bpath_vkey_type(use_load, customer_service, cuts),
            get_bpath_label_type(use_load, customer_service, cuts),
            Base.Order.ForwardOrdering,
        }(Base.Order.ForwardOrdering())
        for starting_node in graph.N_depots,
            current_node in graph.N_nodes
    )
    unexplored_states = SortedSet{get_bpath_vkey_fkey_type(use_load, customer_service, cuts)}()

    for depot in graph.N_depots
        p = create_new_bpath_label(
            use_load,
            customer_service,
            cuts,
            depot,
            data,
            ;
            n_cuts = length(λ),
        )
        vkey = get_vkey(p)
        path_labels[(depot, depot)][vkey] = p 
        push!(unexplored_states, (vkey..., depot, depot,))
    end

    while length(unexplored_states) > 0
        if time_limit < time() - start_time
            throw(TimeLimitException())
        end
        # Retrieve most promising unexplored state
        state = pop!(unexplored_states)
        (current_vkey..., starting_node, current_node) = state
        if !(current_vkey in keys(path_labels[(starting_node, current_node)]))
            continue
        end
        current_path = path_labels[(starting_node, current_node)][current_vkey]
        for next_node in setdiff(outneighbors(graph.G, current_node), current_node)
            (feasible, new_path) = compute_new_bpath(
                current_path,
                data,
                graph,
                modified_costs,
                current_node,
                next_node,
                use_load,
                use_time_windows,
                customer_service,
                cuts,
                ;
                α = data.α,
                β = data.β,
                neighborhoods = neighborhoods,
                n_cuts = length(λ),
                λvals = λvals,
                λcust = λcust,
                λmemory = λmemory,
            )
            !feasible && continue

            new_vkey = get_vkey(new_path)
            if cuts isa NoCuts
                added = add_label_to_collection!(
                    path_labels[(starting_node, next_node)],
                    new_vkey, new_path,
                )
            elseif cuts isa SR3Cuts || cuts isa LmSR3Cuts
                added = add_label_to_collection_cuts!(
                    path_labels[(starting_node, next_node)],
                    new_vkey, new_path, λvals,
                )
            end
            if added && !(next_node in graph.N_depots)
                push!(unexplored_states, (new_vkey..., starting_node, next_node))
            end
        end
    end

    for depot in graph.N_depots
        for path in values(path_labels[(depot, depot)])
            if length(path.nodes) == 1
                path.nodes = [depot, depot]
                path.excesses = [0]
                path.slacks = [0]
            end
        end
    end

    for starting_node in graph.N_depots
        for end_node in setdiff(graph.N_nodes, graph.N_depots)
            delete!(path_labels, (starting_node, end_node))
        end
    end
    for starting_node in graph.N_depots
        for end_node in graph.N_depots
            for path in values(path_labels[(starting_node, end_node)])
                path.cost -= (κ[starting_node] + μ[end_node])
            end
        end
    end

    return path_labels
end


unwrap_path_labels(p::Label) = Label[p]

function unwrap_path_labels(d::AbstractDict)
    u = Label[]
    for v in values(d)
        append!(u, unwrap_path_labels(v))
    end
    return u
end

function get_negative_path_labels_from_path_labels(
    path_labels::Dict{
        Tuple{Int, Int}, 
        T,
    },
) where {T <: AbstractDict}
    return Label[
        path_label
        for path_label in unwrap_path_labels(path_labels)
            if path_label.cost < -1e-4
    ]
end

function subproblem_iteration_benchmark(
    data::EVRPData,
    graph::EVRPGraph,
    κ::Dict{Int, Float64},
    μ::Dict{Int, Float64},
    ν::Vector{Float64},
    λ::Dict{<:Tuple, Float64},
    ;
    load::Bool = false,
    time_windows::Bool = false,
    charge_cost_heterogenous::Bool = false,
    neighborhoods::Union{Nothing, BitMatrix} = nothing,
    ngroute::Bool = false,
    elementary::Bool = true,
    time_limit::Float64 = Inf,
)

    if charge_cost_heterogenous
        throw(ErrorException("Benchmark method does not support heterogenous charging costs."))
    end

    start_time = time()
    use_load = load ? YesLoad() : NoLoad()
    use_time_windows = time_windows ? YesTimeWindows() : NoTimeWindows()
    if ngroute
        customer_service = NgRoute()
        if length(λ) == 0
            cuts = NoCuts()
        elseif keytype(λ) == NTuple{3, Int}
            cuts = SR3Cuts()
        elseif keytype(λ) == Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}}
            cuts = LmSR3Cuts()
        else
            error("Unrecognized key type for λ: $(keytype(λ))")
        end
    else
        customer_service = elementary ? Elementary() : NoService()
        cuts = NoCuts()
    end

    path_labels_result = @timed generate_path_labels_benchmark(
        data, 
        graph,
        κ,
        μ,
        ν,
        use_load,
        use_time_windows,
        customer_service,
        cuts,
        ;
        neighborhoods = neighborhoods,
        λ = λ,
        time_limit = time_limit - (time() - start_time),
    )
    path_labels_time = path_labels_result.time
    negative_path_labels = get_negative_path_labels_from_path_labels(path_labels_result.value)
    negative_path_labels_count = length(negative_path_labels)
    return (negative_path_labels, negative_path_labels_count, path_labels_time)

end

function convert_path_label_to_path(
    path_label::BPathLabel,
    data::EVRPData,
    graph::EVRPGraph,
    ;
    use_load::Bool = false,
)
    p = Path(
        subpaths = Subpath[],
        charging_arcs = ChargingArc[],
        served = zeros(Int, graph.n_customers),
        load = 0,
        arcs = NTuple{2, Int}[],
        customer_arcs = NTuple{2, Int}[],
    )
    states = NTuple{3, Int}[]
    current_subpath = Subpath(
        n_customers = graph.n_customers,
        starting_node = path_label.nodes[1],
        starting_time = 0, 
        starting_charge = graph.B,
    )
    i = path_label.nodes[1]
    for (j, e, s) in zip(path_label.nodes[2:end], path_label.excesses, path_label.slacks)
        current_subpath.current_node = j
        push!(current_subpath.arcs, (i, j))
        current_subpath.starting_time += (e + s)
        current_subpath.starting_charge += (e + s) 
        current_subpath.current_time += (graph.t[i,j] + e + s)
        current_subpath.current_charge += (- graph.q[i,j] + e + s)
        if use_load
            current_subpath.load += graph.d[j]
        end
        if j in graph.N_charging
            push!(
                states, 
                (current_subpath.starting_node, current_subpath.starting_time, current_subpath.starting_charge), 
                (current_subpath.current_node, current_subpath.current_time, current_subpath.current_charge), 
            )
            push!(
                p.subpaths,
                current_subpath,
            )
            current_subpath = Subpath(
                n_customers = graph.n_customers,
                starting_node = j,
                starting_time = current_subpath.current_time, 
                starting_charge = current_subpath.current_charge,
            )
        elseif j in graph.N_customers
            current_subpath.served[j] += 1
        end
        i = j
    end
    push!(
        states, 
        (current_subpath.starting_node, current_subpath.starting_time, current_subpath.starting_charge), 
        (current_subpath.current_node, current_subpath.current_time, current_subpath.current_charge), 
    )
    push!(
        p.subpaths,
        current_subpath,
    )
    for i in 1:(length(states)÷2)-1
        push!(
            p.charging_arcs, 
            ChargingArc(
                states[2*i][1],
                states[2*i][2],
                states[2*i+1][2],
                states[2*i+1][2] - states[2*i][2],
                states[2*i][3],
                states[2*i+1][3],
                states[2*i+1][3] - states[2*i][3],
            )
        )
    end
    p.served = sum(s.served for s in p.subpaths)
    p.load = sum(s.load for s in p.subpaths)
    p.arcs = vcat([s.arcs for s in p.subpaths]...)
    customers = [a[1] for a in p.arcs if a[1] in graph.N_customers]
    p.customer_arcs = collect(zip(customers[1:end-1], customers[2:end]))
    return p
end
