include("decomposition_heuristic.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Test

const GRB_ENV = Gurobi.Env()

args_df = CSV.read("$(@__DIR__)/../experiments/heuristic_benchmark/01a/args.csv", DataFrame)
row_index = 8

begin
    n_depots = args_df[row_index, :n_depots]
    n_customers = args_df[row_index, :n_customers]
    n_charging = args_df[row_index, :n_charging]
    depot_pattern = String(args_df[row_index, :depot_pattern])
    customer_pattern = String(args_df[row_index, :customer_pattern])
    charging_pattern = String(args_df[row_index, :charging_pattern])
    customer_spread = args_df[row_index, :customer_spread]
    xmin = args_df[row_index, :xmin]
    xmax = args_df[row_index, :xmax]
    ymin = args_df[row_index, :ymin]
    ymax = args_df[row_index, :ymax]
    n_vehicles = args_df[row_index, :n_vehicles]
    T = args_df[row_index, :T]
    B = args_df[row_index, :B]
    μ = args_df[row_index, :μ]
    seed = args_df[row_index, :seed]
    travel_cost_coeff = args_df[row_index, :travel_cost_coeff]
    charge_cost_coeff = args_df[row_index, :charge_cost_coeff]
    load_scale = args_df[row_index, :load_scale]
    load_shape = args_df[row_index, :load_shape]
    load_tolerance = args_df[row_index, :load_tolerance]
    batch = args_df[row_index, :batch]
    permissiveness = args_df[row_index, :permissiveness]

    use_load = args_df[row_index, :use_load]
    use_time_windows = args_df[row_index, :use_time_windows]

    heuristic_use_adaptive_ngroute = args_df[row_index, :heuristic_use_adaptive_ngroute]
    heuristic_use_SR3_cuts = args_df[row_index, :heuristic_use_SR3_cuts]
    heuristic_use_lmSR3_cuts = args_df[row_index, :heuristic_use_lmSR3_cuts]
    method = String(args_df[row_index, :method])
    ngroute_neighborhood_charging_size = String(args_df[row_index, :ngroute_neighborhood_charging_size])
    use_adaptive_ngroute = args_df[row_index, :use_adaptive_ngroute]
    use_SR3_cuts = args_df[row_index, :use_SR3_cuts]
    use_lmSR3_cuts = args_df[row_index, :use_lmSR3_cuts]    
    max_SR3_cuts = args_df[row_index, :max_SR3_cuts]

    data = generate_instance(
        n_depots = n_depots,
        n_customers = n_customers,
        n_charging = n_charging,
        n_vehicles = n_vehicles,
        depot_pattern = depot_pattern,
        customer_pattern = customer_pattern,
        charging_pattern = charging_pattern,
        customer_spread = customer_spread,
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        T = T,
        seed = seed,
        B = B,
        μ = μ,
        travel_cost_coeff = travel_cost_coeff,
        charge_cost_coeff = charge_cost_coeff,
        load_scale = load_scale,
        load_shape = load_shape,
        load_tolerance = load_tolerance,
        batch = batch,
        permissiveness = permissiveness,
        ;
        data_dir = "../../../data/",
    )
    graph = generate_graph_from_data(data)
end

heuristic_run = @timed path_formulation_decomposition_heuristic(
    data, graph;
    elementary = false,
    ngroute = true,
    use_adaptive_ngroute = heuristic_use_adaptive_ngroute,
    use_SR3_cuts = heuristic_use_SR3_cuts,
    use_lmSR3_cuts = heuristic_use_lmSR3_cuts,
    max_SR3_cuts = max_SR3_cuts,
    time_heuristic_slack = 0.9,
);

(
    hCGIP_result, 
    heuristic_results,
    time_heuristic_slack,
) = heuristic_run.value;


plot_solution(hCGIP_result, data)
plot_solution(heuristic_results, data)

optimal_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    Env = GRB_ENV,
    method = method,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small",
    ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size,
    verbose = true,
    use_adaptive_ngroute = use_adaptive_ngroute,
    use_SR3_cuts = use_SR3_cuts,
    use_lmSR3_cuts = use_lmSR3_cuts,
    max_SR3_cuts = max_SR3_cuts,
    time_limit = 3600.0,
);
(
    CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
    some_paths, model, z, SR3_constraints
) = optimal_run.value;

plot_solution(CGIP_all_results[end], data)

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
    verbose = true,
    # time_limit = time_limit - (time() - start_time),
    # max_iters = max_iters,
)