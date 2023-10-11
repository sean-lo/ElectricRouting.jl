include("decomposition_heuristic.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Cthulhu
using Test

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 20,
    n_charging = 5,
    n_vehicles = 6,
    depot_pattern = "grid",
    customer_pattern = "random_box",
    charging_pattern = "grid_clipped",
    customer_spread = 0.1,
    xmin = 0.0,
    xmax = 2.0,
    ymin = 0.0,
    ymax = 2.0,
    T = 72000,
    seed = 2,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)


(
    CGIP_result, 
    heuristic_results,
) = path_formulation_decomposition_heuristic_new(
    data, graph;
    elementary = false,
    ngroute = true,
    ngroute_alt = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
    use_lmSR3_cuts = false,
    time_heuristic_slack = 0.9,
);

CGIP_result

plot_solution(CGIP_result, data)
plot_solution(heuristic_results, data)


(
    CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
    some_paths, model, z, SR3_constraints
) = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    elementary = false,
    ngroute = true,
    ngroute_alt = true,    
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = true,
);

plot_solution(CGLP_all_results[end], data)

CGIP_all_results[end]["objective"]
heuristic_results["objective"]


# debug

include("decomposition_heuristic.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

(
    CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
    some_paths, model, z, SR3_constraints
) = path_formulation_decomposition_heuristic_new(
    data, graph;
    elementary = false,
    ngroute = true,
    ngroute_alt = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
    use_lmSR3_cuts = false,
    time_heuristic_slack = 0.9,
);

neighborhoods = CG_all_neighborhoods[end]

CGLP_results = CGLP_all_results[end]
max_SR3_cuts = 10
generated_WSR3_list = enumerate_violated_path_WSR3_inequalities(
    CGLP_results["paths"], 
    graph,
)
generated_SR3_list = enumerate_violated_path_SR3_inequalities(
    CGLP_results["paths"],
    graph,
)
implemented_SR3_list = sample(
    generated_SR3_list, 
    Weights([val for (val, _) in generated_SR3_list]),
    max_SR3_cuts, 
    replace = false,
)
add_SR3_constraints_to_path_model!(
    model, z, some_paths, 
    SR3_constraints, implemented_SR3_list, 
)

path_costs = compute_path_costs(
    data, graph, 
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)

time_heuristic_slack = 0.9
T_heuristic = Int(round(time_heuristic_slack * (graph.T + graph.B) * graph.μ / (1 + graph.μ)))

optimize!(model)
CGLP_results = Dict(
    "objective" => objective_value(model),
    "z" => Dict(
        (key, p) => value.(z[(key, p)])
        for (key, p) in keys(z)
    ),
    "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
    "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
    "ν" => dual.(model[:ν]).data,
)
(negative_base_labels, _, base_labels_time) = subproblem_iteration_nocharge(
    data, graph, 
    CGLP_results["κ"], 
    CGLP_results["μ"], 
    CGLP_results["ν"], 
    T_heuristic,
    ;
    neighborhoods = neighborhoods,
    ngroute = true,
    ngroute_alt = true,
    elementary = false,
    # time_limit = (time_limit - (time() - start_time)),
)
generated_paths = get_paths_from_negative_base_labels(
    graph, negative_base_labels,
)

negative_base_labels

[
    compute_path_label_modified_cost(
        p, data, graph,
        CGLP_results["κ"], 
        CGLP_results["μ"], 
        CGLP_results["ν"], 
        ; 
        verbose = true
    )
    for p in negative_base_labels
]


CGLP_results, CG_params = path_formulation_column_generation_nocharge!(
    model, z, 
    data, graph,
    some_paths, path_costs, path_service,
    printlist,
    ;
    time_heuristic_slack = 0.9,
    elementary = false,
    neighborhoods = neighborhoods,
    ngroute = true,
    ngroute_alt = true,
    verbose = true,
    # time_limit = time_limit - (time() - start_time),
    # max_iters = max_iters,
)