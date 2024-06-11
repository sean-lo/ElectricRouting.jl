using JuMP, Gurobi
using Clp
using Random, Distributions
using Combinatorics, Suppressor
using BenchmarkTools

all_costs = collect(permutations(1:4))
B = 100
n = 10
charges = min.(round.(rand(Normal(5, 2), n+1)) .* 10, 90)
costs = round.(rand(Normal(5, 2), n)) .* 10
@btime solve_charge_rebalancing_recursive(B, Int.(charges), Int.(costs))
@btime solve_charge_rebalancing_analytic(B, Int.(charges), Int.(costs))

function solve_charge_rebalancing_recursive(
    B::Int,
    m::Int,
    charges::Vector{Int},
    costs::Vector{Int},
)
    if m ≤ 1
        return Int[]
    end
    if sum(charges) ≤ B
        return zeros(Int, m - 1)
    end

    i = m - argmin(reverse(costs))
    val = min(
        sum(charges) - B, # l_{m-1}
        sum(charges[1:i]), # u_i
        sum(charges[i+1:end]),
        B,
    )
    return vcat(
        solve_charge_rebalancing_recursive(B, i, charges[1:i], costs[1:i-1]),
        [val],
        solve_charge_rebalancing_recursive(B, m-i, charges[i+1:end], costs[i+1:end]),
    )
end

function solve_charge_rebalancing(B::Int, charges::Vector{Int}, costs::Vector{Int})
    # model = Model(Gurobi.Optimizer)
    model = Model(Clp.Optimizer)
    set_silent(model)
    set_string_names_on_creation(model, false)
    n = length(costs)
    @variable(model, x[1:n] ≥ 0)
    @constraint(model, [i=1:n], sum(x[j] for j in 1:i) ≥ sum(charges[j] for j in 1:(i+1)) - B)
    @constraint(model, [i=1:n], sum(x[j] for j in 1:i) ≤ sum(charges[j] for j in 1:i))
    @objective(model, Min, sum(x[i] * costs[i] for i in 1:n))
    optimize!(model)
    x_sol = value.(model[:x])
    return Int.(x_sol)
end


function rebalance_path_label!(p::PathLabel, data::EVRPData)
    charge_costs = [
        data.charge_cost_coeffs[s.nodes[end]]
        for s in p.subpath_labels[1:end-1]
    ]
    charges = [
        s.charge_taken
        for s in p.subpath_labels
    ]
    n = length(charge_costs)
    m = @suppress Model(Clp.Optimizer)
    # m = @suppress Model(() -> Gurobi.Optimizer(env))
    JuMP.set_silent(m)
    JuMP.set_string_names_on_creation(m, false)
    @variable(m, x[1:n] ≥ 0)
    @constraint(m, [i=1:n], sum(x[j] for j in 1:i) ≥ sum(charges[j] for j in 1:(i+1)) - data.B)
    @constraint(m, [i=1:n], sum(x[j] for j in 1:i) ≤ sum(charges[j] for j in 1:i))
    @objective(m, Min, sum(x[i] * charge_costs[i] for i in 1:n))
    optimize!(m)
    p.charging_actions = Int.(value.(m[:x]))
    return
end

include("utils.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")

using CSV, DataFrames
using Test

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 12,
    n_charging = 5,
    n_vehicles = 4,
    depot_pattern = "grid",    
    customer_pattern = "random_box",
    charging_pattern = "grid_clipped",
    customer_spread = 0.1,
    xmin = 0.0, xmax = 2.0,
    ymin = 0.0, ymax = 2.0,
    T = 60000,
    seed = 1,
    B = 20000,
    μ = 5,
    travel_cost_coeff = 70,
    charge_cost_coeff = 0,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
    charge_cost_heterogenous = false,
    # charge_cost_stddev = 10.0,
)
graph = generate_graph_from_data(data)
plot_instance(data)



SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
SR3_list = Tuple{Float64, NTuple{3, Int}}[]
artificial_paths = generate_artificial_paths(data, graph)
some_paths = deepcopy(artificial_paths)
path_costs = compute_path_costs(data, graph, some_paths)
path_service = compute_path_service(graph, some_paths)
model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
)
z




optimize!(model)

CGLP_results = Dict(
    "errored" => false,
    "objective" => objective_value(model),
    "z" => Dict(
        (key, p) => value.(z[(key, p)])
        for (key, p) in keys(z)
    ),
    "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
    "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
    "ν" => dual.(model[:ν]).data,
    "λ" => Dict{keytype(SR3_constraints), Float64}(),
)
[
    (key, p, CGLP_results["z"][(key, p)], some_paths[key][p].artificial)
    for (key, p) in keys(z)
        if CGLP_results["z"][(key, p)] > 1e-3
]
[
    value.(z[(key, p)])
    for key in keys(artificial_paths)
        for p in 1:length(artificial_paths[key])
]

CGLP_results["artificial"] = any(
    value.(z[(key, p)]) > 1e-3
    for key in keys(artificial_paths)
        for p in 1:length(artificial_paths[key])
)

(negative_full_labels, _, base_labels_time, full_labels_time) = subproblem_iteration_ours(
    data, graph, 
    CGLP_results["κ"], 
    CGLP_results["μ"], 
    CGLP_results["ν"], 
    CGLP_results["λ"], 
    ;
    charge_cost_heterogenous = false,
    neighborhoods = nothing,
    ngroute = false,
    elementary = false,
)
generated_paths = get_paths_from_negative_path_labels(
    data, graph, negative_full_labels,
)
mp_constraint_time = add_paths_to_path_model!(
    model,
    z,
    some_paths, 
    path_costs,
    path_service,
    generated_paths,
    SR3_constraints,
    data, graph,
)



for (charge_cost_heterogenous, elementary, ngroute, use_adaptive_ngroute, use_SR3_cuts, use_lmSR3_cuts) in [
    (false, false, false, false, false, false), # FIXME 
    (false, true, false, false, false, false), # FIXME
    (false, false, true, false, false, false),
    (false, false, true, true, false, false),
    (false, false, true, true, true, false),
    (false, false, true, true, true, true),
]
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = @time path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        charge_cost_heterogenous = false,
        method = "benchmark",
        elementary = elementary,
        ngroute = ngroute,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        verbose = true,
        use_adaptive_ngroute = use_adaptive_ngroute,
        use_SR3_cuts = use_SR3_cuts,
        use_lmSR3_cuts = use_lmSR3_cuts,
    );
end


data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 12,
    n_charging = 5,
    n_vehicles = 4,
    depot_pattern = "grid",    
    customer_pattern = "random_box",
    charging_pattern = "grid_clipped",
    customer_spread = 0.1,
    xmin = 0.0, xmax = 2.0,
    ymin = 0.0, ymax = 2.0,
    T = 60000,
    seed = 1,
    B = 20000,
    μ = 5,
    travel_cost_coeff = 70,
    charge_cost_coeff = 30,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
    charge_cost_heterogenous = false,
    # charge_cost_stddev = 10.0,
)
graph = generate_graph_from_data(data)
plot_instance(data)

for (charge_cost_heterogenous, elementary, ngroute, use_adaptive_ngroute, use_SR3_cuts, use_lmSR3_cuts) in [
    (false, false, false, false, false, false),
    (false, true, false, false, false, false),
    (false, false, true, false, false, false),
    (false, false, true, true, false, false),
    (false, false, true, true, true, false),
    (false, false, true, true, true, true),
]
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = @time path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        charge_cost_heterogenous = charge_cost_heterogenous,
        method = "ours",
        elementary = elementary,
        ngroute = ngroute,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        verbose = true,
        use_adaptive_ngroute = use_adaptive_ngroute,
        use_SR3_cuts = use_SR3_cuts,
        use_lmSR3_cuts = use_lmSR3_cuts,
    );
end


data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 12,
    n_charging = 5,
    n_vehicles = 4,
    depot_pattern = "grid",    
    customer_pattern = "random_box",
    charging_pattern = "grid_clipped",
    customer_spread = 0.1,
    xmin = 0.0, xmax = 2.0,
    ymin = 0.0, ymax = 2.0,
    T = 60000,
    seed = 1,
    B = 20000,
    μ = 5,
    travel_cost_coeff = 70,
    charge_cost_coeff = 30,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
    charge_cost_heterogenous = true,
    charge_cost_stddev = 10.0,
)
graph = generate_graph_from_data(data)
plot_instance(data)
for (charge_cost_heterogenous, elementary, ngroute, use_adaptive_ngroute, use_SR3_cuts, use_lmSR3_cuts) in [
    (true, false, false, false, false, false),
    (true, true, false, false, false, false),
    (true, false, true, false, false, false),
    (true, false, true, true, false, false),
    (true, false, true, true, true, false),
    (true, false, true, true, true, true),
]
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = @time path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        charge_cost_heterogenous = charge_cost_heterogenous,
        method = "ours",
        elementary = elementary,
        ngroute = ngroute,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        verbose = true,
        use_adaptive_ngroute = use_adaptive_ngroute,
        use_SR3_cuts = use_SR3_cuts,
        use_lmSR3_cuts = use_lmSR3_cuts,
    );
end




(charge_cost_heterogenous, elementary, ngroute, use_adaptive_ngroute, use_SR3_cuts, use_lmSR3_cuts) = 
(true, false, true, false, false, false)


neighborhoods = compute_ngroute_neighborhoods(
    graph,
    Int(ceil(sqrt(graph.n_customers))),
    ;
    depots_size = "small", 
    charging_size = "small", 
)
some_paths = generate_artificial_paths(data, graph)
path_costs = compute_path_costs(
    data, graph, 
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)
model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
    ; 
    # Env = Env,
)
if use_lmSR3_cuts
    SR3_constraints = Dict{
        Tuple{NTuple{3, Int}, Tuple{Vararg{Int}}},
        ConstraintRef,
    }()
    SR3_list = Tuple{Float64, NTuple{3, Int}, Tuple{Vararg{Int}}}[]
else
    SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
    SR3_list = Tuple{Float64, NTuple{3, Int}}[]
end
CGLP_all_results = []
CGIP_all_results = []
all_params = []
CG_all_params = []
if ngroute
    CG_all_neighborhoods = BitMatrix[]
else
    CG_all_neighborhoods = nothing
end
printlist = String[]

CG_params = Dict{String, Any}()
CG_params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
CG_params["objective"] = Float64[]
CG_params["κ"] = Dict{Int, Float64}[]
CG_params["μ"] = Dict{Int, Float64}[]
CG_params["ν"] = Vector{Float64}[]
CG_params["λ"] = Dict{keytype(SR3_constraints), Float64}[]
CG_params["lp_relaxation_solution_time_taken"] = Float64[]
CG_params["sp_base_time_taken"] = Float64[]
CG_params["sp_full_time_taken"] = Float64[]
CG_params["sp_total_time_taken"] = Float64[]
CG_params["lp_relaxation_constraint_time_taken"] = Float64[]
CG_params["number_of_new_paths"] = Int[]
CG_params["converged"] = false


@suppress optimize!(model)
CGLP_results = Dict(
    "objective" => objective_value(model),
    "z" => Dict(
        (key, p) => value.(z[(key, p)])
        for (key, p) in keys(z)
    ),
    "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
    "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
    "ν" => dual.(model[:ν]).data,
    "λ" => Dict{keytype(SR3_constraints), Float64}(
        S => dual(SR3_constraints[S])
        for S in keys(SR3_constraints)
    ),
)
push!(CG_params["objective"], CGLP_results["objective"])
push!(CG_params["κ"], CGLP_results["κ"])
push!(CG_params["μ"], CGLP_results["μ"])
push!(CG_params["ν"], CGLP_results["ν"])
push!(CG_params["λ"], CGLP_results["λ"])


base_labels_result = @timed generate_base_labels_ngroute(
    data, graph, neighborhoods, 
    CGLP_results["κ"], 
    CGLP_results["μ"], 
    CGLP_results["ν"], 
    ;
    # time_limit = time_limit - (time() - start_time),
)

for (key, value) in pairs(base_labels_result.value)
    println(key)
    println(length(value))
end

full_labels = Dict(
    (starting_node, current_node) => SortedDict{
        Tuple{Float64, Int, Int, Vector{Int}, BitVector}, 
        PathLabel,
        Base.Order.ForwardOrdering,
    }(Base.Order.ForwardOrdering())
    for starting_node in graph.N_depots,
        current_node in graph.N_depots_charging
)
unexplored_states = SortedSet{Tuple{Float64, Int, Int, Vector{Int}, BitVector, Int, Int}}()
for depot in graph.N_depots
    node_labels = falses(graph.n_nodes)
    node_labels[depot] = true
    # label key here has the following fields:
    # 0) reduced cost
    # 1) time
    # 2) - charge
    # 3) - rebalancing slacks
    # 4) whether i-th node is in forward ng-set
    key = (0.0, 0, -graph.B, zeros(Int, data.charge_cost_nlevels - 1), node_labels,)
    full_labels[(depot, depot)][key] = PathLabel(
        0.0,
        BaseSubpathLabel[],
        Int[],
        zeros(Int, data.charge_cost_nlevels),
        zeros(Int, graph.n_customers),
    )
    push!(
        unexplored_states, 
        (
            key..., 
            depot, # starting_node
            depot, # current_node
        )
    )
end

unexplored_states
full_labels[16,16]

# begin
    state = pop!(unexplored_states)
    starting_node = state[end-1]
    current_node = state[end]
    current_set = state[end-2]
    current_key = state[1:end-2]
    if !(current_key in keys(full_labels[(starting_node, current_node)]))
        continue
    end
    (current_path_reduced_cost, current_time, negative_current_charge, negative_rebalancing_slacks, current_set) = current_key
    current_path = full_labels[(starting_node, current_node)][current_key]
    
    
    for next_node in graph.N_depots_charging
        for s in values(base_labels_result.value[(current_node, next_node)])
            # ngroute stitching subpaths check
            (feasible, new_set) = ngroute_extend_partial_path_check(
                neighborhoods, current_set, s.nodes,
            )
            !feasible && continue
            (feasible, new_path, end_time, end_charge, new_rebalancing_slacks) = compute_new_path_heterogenous_charging(
                current_path, 
                s, 
                current_node,
                current_time,
                - negative_current_charge,
                - negative_rebalancing_slacks,
                next_node, 
                data, 
                graph,
            )
            !feasible && continue

            new_key = (
                new_path.cost,
                end_time, 
                - end_charge,
                - new_rebalancing_slacks,
                new_set,
            )
            println(new_key)
            added = add_label_to_collection!(
                full_labels[(starting_node, next_node)],
                new_key, new_path,
                ;
            )
            if added && next_node in graph.N_charging
                new_state = (new_key..., starting_node, next_node)
                push!(unexplored_states, new_state)
            end
        end
    end
# end

base_labels_result.value[19,19]
full_labels
unexplored_state