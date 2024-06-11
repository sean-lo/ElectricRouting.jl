
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Test

const GRB_ENV = Gurobi.Env()

args_df = CSV.read("$(@__DIR__)/../experiments/heuristic_benchmark/02a/args.csv", DataFrame)
row_index = 1

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
    charge_cost_coeff_increment = args_df[row_index, :charge_cost_coeff_increment]
    charge_cost_nlevels = args_df[row_index, :charge_cost_nlevels]
    load_scale = args_df[row_index, :load_scale]
    load_shape = args_df[row_index, :load_shape]
    load_tolerance = args_df[row_index, :load_tolerance]
    batch = args_df[row_index, :batch]
    permissiveness = args_df[row_index, :permissiveness]

    use_load = args_df[row_index, :use_load]
    use_time_windows = args_df[row_index, :use_time_windows]

    method = String(args_df[row_index, :method])
    ngroute_neighborhood_charging_size = String(args_df[row_index, :ngroute_neighborhood_charging_size])
    use_adaptive_ngroute = args_df[row_index, :use_adaptive_ngroute]
    use_SR3_cuts = args_df[row_index, :use_SR3_cuts]
    use_lmSR3_cuts = args_df[row_index, :use_lmSR3_cuts]    
    max_SR3_cuts = args_df[row_index, :max_SR3_cuts]


    data_het = generate_instance(
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
        charge_cost_heterogenous = true,
        charge_cost_nlevels = charge_cost_nlevels,
        charge_cost_coeff_increment = charge_cost_coeff_increment,
    )
    graph_het = generate_graph_from_data(data_het)

    data_hom = generate_instance(
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
    graph_hom = generate_graph_from_data(data_hom)
end

optimal_EVRPhet_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data_het, graph_het,
    ;
    Env = GRB_ENV,
    method = method,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph_het.n_customers))),
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
    CGLP_all_results_het, CGIP_all_results_het, CG_all_params_het, CG_all_neighborhoods_het, all_params_het, printlist_het, 
    some_paths_het, model_het, z_het, SR3_constraints_het
) = optimal_EVRPhet_run.value;

some_paths_metrics_het = compute_path_metrics(some_paths_het)
all_params_het_df = DataFrame(all_params_het)
ind_het = findfirst(x -> x in ["use_SR3_cuts", "use_lmSR3_cuts"], all_params_het_df.method)
if isnothing(ind_het)
    ind_het = length(all_params_het)
end

LP_objective_first_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_het[1]["paths"]
)
LP_objective_last_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_het[ind_het]["paths"]
)
LP_objective_last_SR3_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_het[end]["paths"]
)
IP_objective_first_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_het[1]["paths"]
)
IP_objective_last_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_het[ind_het]["paths"]
)
IP_objective_last_SR3_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_het[end]["paths"]
)


optimal_EVRPhom_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data_hom, graph_hom,
    ;
    Env = GRB_ENV,
    method = method,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph_hom.n_customers))),
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
    CGLP_all_results_hom, CGIP_all_results_hom, CG_all_params_hom, CG_all_neighborhoods_hom, all_params_hom, printlist_hom, 
    some_paths_hom, model_hom, z_hom, SR3_constraints_hom
) = optimal_EVRPhom_run.value;

some_paths_metrics_hom = compute_path_metrics(some_paths_hom)
all_params_hom_df = DataFrame(all_params_hom)
ind_hom = findfirst(x -> x in ["use_SR3_cuts", "use_lmSR3_cuts"], all_params_hom_df.method)
if isnothing(ind_hom)
    ind_hom = length(all_params_hom)
end

LP_objective_first_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_hom[1]["paths"]
)
LP_objective_last_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_hom[ind_hom]["paths"]
)
LP_objective_last_SR3_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_hom[end]["paths"]
)
IP_objective_first_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_hom[1]["paths"]
)
IP_objective_last_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_hom[ind_hom]["paths"]
)
IP_objective_last_SR3_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_hom[end]["paths"]
)

records = [
    (
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
        charge_cost_coeff_increment = charge_cost_coeff_increment,
        charge_cost_nlevels = charge_cost_nlevels,
        load_scale = load_scale,
        load_shape = load_shape,
        load_tolerance = load_tolerance,
        batch = batch,
        permissiveness = permissiveness,
        use_load = use_load,
        use_time_windows = use_time_windows,
        method = method,
        elementary = false,
        ngroute = true,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small",
        ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size,
        
        use_adaptive_ngroute = true,
        use_SR3_cuts = true,
        use_lmSR3_cuts = use_lmSR3_cuts,
        # Time taken
        time_taken_het = optimal_EVRPhet_run.time,
        time_taken_hom = optimal_EVRPhom_run.time,
        # Objective values and gaps
        LP_objective_first_het = LP_objective_first_het,
        LP_objective_last_het = LP_objective_last_het,
        LP_objective_last_SR3_het = LP_objective_last_SR3_het,
        IP_objective_first_het = IP_objective_first_het,
        IP_objective_last_het = IP_objective_last_het,
        IP_objective_last_SR3_het = IP_objective_last_SR3_het,
        LP_objective_first_hom = LP_objective_first_hom,
        LP_objective_last_hom = LP_objective_last_hom,
        LP_objective_last_SR3_hom = LP_objective_last_SR3_hom,
        IP_objective_first_hom = IP_objective_first_hom,
        IP_objective_last_hom = IP_objective_last_hom,
        IP_objective_last_SR3_hom = IP_objective_last_SR3_hom,
    )
]

LP_objective_first_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_het[1]["paths"]
)
LP_objective_last_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_het[ind_het]["paths"]
)
LP_objective_last_SR3_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_het[end]["paths"]
)
IP_objective_first_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_het[1]["paths"]
)
IP_objective_first_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_het[ind_het]["paths"]
)
IP_objective_first_het = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_het[end]["paths"]
)

CGLP_all_results_het[1]["objective"]
CGLP_all_results_het[ind_het]["objective"]
CGLP_all_results_het[end]["objective"]
CGIP_all_results_het[1]["objective"]
CGIP_all_results_het[ind_het]["objective"]
CGIP_all_results_het[end]["objective"]


LP_objective_first_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_hom[1]["paths"]
)
LP_objective_last_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_hom[ind_hom]["paths"]
)
LP_objective_last_SR3_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGLP_all_results_hom[end]["paths"]
)
IP_objective_first_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_hom[1]["paths"]
)
IP_objective_first_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_hom[ind_hom]["paths"]
)
IP_objective_first_hom = sum(
    val * compute_path_cost(data_het, graph_het, p)
    for (val, p) in CGIP_all_results_hom[end]["paths"]
)