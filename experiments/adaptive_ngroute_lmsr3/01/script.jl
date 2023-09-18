include("../../../src/path_formulation.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

using JuMP, Gurobi

const GRB_ENV = Gurobi.Env()

function run_instance(
    args_df, 
    row_index,
    ;
    write_log::Bool = true,
)
    # Get paramters from args_df at row row_index
    n_depots = args_df[row_index, :n_depots]
    n_customers = args_df[row_index, :n_customers]
    n_charging = args_df[row_index, :n_charging]
    n_vehicles = args_df[row_index, :n_vehicles]
    depot_pattern = String(args_df[row_index, :depot_pattern])
    customer_pattern = String(args_df[row_index, :customer_pattern])
    charging_pattern = String(args_df[row_index, :charging_pattern])
    shrinkage_depots = args_df[row_index, :shrinkage_depots]
    shrinkage_charging = args_df[row_index, :shrinkage_charging]
    T = args_df[row_index, :T]
    B = args_df[row_index, :B]
    seed = args_df[row_index, :seed]
    μ = args_df[row_index, :μ]
    travel_cost_coeff = args_df[row_index, :travel_cost_coeff]
    charge_cost_coeff = args_df[row_index, :charge_cost_coeff]
    load_scale = args_df[row_index, :load_scale]
    load_shape = args_df[row_index, :load_shape]
    load_tolerance = args_df[row_index, :load_tolerance]
    batch = args_df[row_index, :batch]
    permissiveness = args_df[row_index, :permissiveness]

    use_load = args_df[row_index, :use_load]
    use_time_windows = args_df[row_index, :use_time_windows]
    formulation = String(args_df[row_index, :formulation])
    method = String(args_df[row_index, :method])
    ngroute = args_df[row_index, :ngroute]
    ngroute_alt = args_df[row_index, :ngroute_alt]
    ngroute_neighborhood_depots_size = String(args_df[row_index, :ngroute_neighborhood_depots_size])
    ngroute_neighborhood_charging_size = String(args_df[row_index, :ngroute_neighborhood_charging_size])
    use_lmSR3_cuts = args_df[row_index, :use_lmSR3_cuts]

    data = generate_instance(
        n_depots = n_depots,
        n_customers = n_customers,
        n_charging = n_charging,
        n_vehicles = n_vehicles,
        depot_pattern = depot_pattern,    
        customer_pattern = customer_pattern,
        charging_pattern = charging_pattern,
        shrinkage_depots = shrinkage_depots,
        shrinkage_charging = shrinkage_charging,
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

    run = @timed @suppress path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        ngroute_alt = ngroute_alt,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = ngroute_neighborhood_depots_size,
        ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size,
        use_adaptive_ngroute = true,
        use_SR3_cuts = true,
        use_lmSR3_cuts = use_lmSR3_cuts,
    );
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = run.value;
    # last non-SR3 iteration index
    all_params_df = DataFrame(all_params)
    ind = findfirst(x -> x == "use_SR3_cuts", all_params_df.method)
    if isnothing(ind)
        ind = length(all_params)
    end

    records = [
        (
            n_depots = n_depots,
            n_customers = n_customers,
            n_charging = n_charging,
            n_vehicles = n_vehicles,
            depot_pattern = depot_pattern,    
            customer_pattern = customer_pattern,
            charging_pattern = charging_pattern,
            shrinkage_depots = shrinkage_depots,
            shrinkage_charging = shrinkage_charging,
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
            use_load = use_load,
            use_time_windows = use_time_windows,
            formulation = formulation,
            method = method,
            ngroute = ngroute,
            ngroute_alt = ngroute_alt,
            ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
            ngroute_neighborhood_depots_size = ngroute_neighborhood_depots_size,
            ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size,
            use_adaptive_ngroute = true,
            use_SR3_cuts = true,
            use_lmSR3_cuts = use_lmSR3_cuts,
            # Time taken
            time_taken_first = CG_all_params[1]["time_taken"],
            time_taken_total = sum(all_params_df.CG_time_taken[1:ind]),
            time_taken_total_SR3 = sum(all_params_df.CG_time_taken),
            # n_iterations
            n_iterations = length(CG_all_params),
            n_CG_iterations = sum(CG_params["counter"] for CG_params in CG_all_params[1:ind]),
            n_CG_iterations_SR3 = sum(CG_params["counter"] for CG_params in CG_all_params),
            # Time taken (subproblem)
            sp_time_taken_mean_first = CG_all_params[1]["sp_time_taken_mean"],
            sp_time_taken_mean_last = CG_all_params[ind]["sp_time_taken_mean"],
            sp_time_taken_mean_last_SR3 = CG_all_params[end]["sp_time_taken_mean"],
            # Objective values and gaps
            LP_objective_first = all_params_df.CGLP_objective[1],
            LP_objective_last = all_params_df.CGLP_objective[ind],
            LP_objective_last_SR3 = all_params_df.CGLP_objective[end],
            IP_objective_first = all_params_df.CGIP_objective[1],
            IP_objective_last = all_params_df.CGIP_objective[ind],
            IP_objective_last_SR3 = all_params_df.CGIP_objective[end],
            LP_IP_gap_first = all_params_df.CG_LP_IP_gap[1],
            LP_IP_gap_last = all_params_df.CG_LP_IP_gap[ind],
            LP_IP_gap_last_SR3 = all_params_df.CG_LP_IP_gap[end],
            # Converged stats
            converged_first = (all_params_df.CG_LP_IP_gap[1] ≈ 0.0),
            converged_last = (all_params_df.CG_LP_IP_gap[ind] ≈ 0.0),
            converged_last_SR3 = (all_params_df.CG_LP_IP_gap[end] ≈ 0.0),
            # Length metrics 
            lp_mean_subpath_length = get(CGLP_all_results[end], "mean_subpath_length", missing),
            lp_weighted_mean_subpath_length = get(CGLP_all_results[end], "weighted_mean_subpath_length", missing),
            lp_mean_path_length = get(CGLP_all_results[end], "mean_path_length", missing),
            lp_weighted_mean_path_length = get(CGLP_all_results[end], "weighted_mean_path_length", missing),
            lp_mean_ps_length = get(CGLP_all_results[end], "mean_ps_length", missing),
            lp_weighted_mean_ps_length = get(CGLP_all_results[end], "weighted_mean_ps_length", missing),
            ip_mean_subpath_length = get(CGIP_all_results[end], "mean_subpath_length", missing),
            ip_weighted_mean_subpath_length = get(CGIP_all_results[end], "weighted_mean_subpath_length", missing),
            ip_mean_path_length = get(CGIP_all_results[end], "mean_path_length", missing),
            ip_weighted_mean_path_length = get(CGIP_all_results[end], "weighted_mean_path_length", missing),
            ip_mean_ps_length = get(CGIP_all_results[end], "mean_ps_length", missing),
            ip_weighted_mean_ps_length = get(CGIP_all_results[end], "weighted_mean_ps_length", missing),
        )
    ]
    if write_log
        CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
    end
    println("Row $row_index processed: $(args_df[row_index, :])")
end

# simple test case to quickly compile 

begin
    test_args_df = DataFrame(CSV.File("$(@__DIR__)/test_args.csv"))
    for i in 1:nrow(test_args_df)
        run_instance(test_args_df, i, write_log = false)
    end
end

println("Compilation complete.")

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

task_index = parse(Int, ARGS[1]) + 1
n_tasks = parse(Int, ARGS[2])

println("Processing rows: $(collect(task_index:n_tasks:size(args_df, 1)))")

for row_index in task_index:n_tasks:size(args_df, 1)
    run_instance(args_df, row_index, write_log = true)
end