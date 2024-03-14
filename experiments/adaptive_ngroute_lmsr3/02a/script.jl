using Pkg
Pkg.activate("../..")
println(Pkg.status())

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
    time_limit,
    ;
    write_log::Bool = true,
    write_results::Bool = true,
)
    # Get paramters from args_df at row row_index
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
    
    method = String(args_df[row_index, :method])
    ngroute_neighborhood_charging_size = String(args_df[row_index, :ngroute_neighborhood_charging_size])
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

    run = @timed @suppress path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
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
        
        use_adaptive_ngroute = true,
        use_SR3_cuts = true,
        use_lmSR3_cuts = use_lmSR3_cuts,
        max_SR3_cuts = max_SR3_cuts,
        time_limit = time_limit,
    );
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = run.value;

    some_paths_metrics = compute_path_metrics(some_paths)
    
    all_params_df = DataFrame(all_params)
    ind = findfirst(x -> x in ["use_SR3_cuts", "use_lmSR3_cuts"], all_params_df.method)
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
            time_taken_first = all_params_df[1, :CG_time_taken],
            time_taken_total = sum(all_params_df[1:ind, :CG_time_taken]),
            time_taken_total_SR3 = sum(all_params_df[1:end, :CG_time_taken]),
            sp_time_taken_mean_first = CG_all_params[1]["sp_time_taken_mean"],
            sp_time_taken_mean_last = CG_all_params[ind]["sp_time_taken_mean"],
            sp_time_taken_mean_last_SR3 = CG_all_params[end]["sp_time_taken_mean"],
            # n_iterations
            converged = all_params_df[end, :converged],
            time_limit = time_limit,
            time_limit_reached = all_params_df[end, :time_limit_reached],
            n_iterations = length(CG_all_params),
            n_CG_iterations = sum(CG_params["counter"] for CG_params in CG_all_params),
            # Objective values and gaps
            LP_objective_first = all_params[1]["CGLP_objective"],
            LP_objective_last = all_params[ind]["CGLP_objective"],
            LP_objective_last_SR3 = all_params[end]["CGLP_objective"],
            IP_objective_first = all_params[1]["CGIP_objective"],
            IP_objective_last = all_params[ind]["CGIP_objective"],
            IP_objective_last_SR3 = all_params[end]["CGIP_objective"],
            LP_IP_gap_first = all_params[1]["CG_LP_IP_gap"],
            LP_IP_gap_last = all_params[ind]["CG_LP_IP_gap"],
            LP_IP_gap_last_SR3 = all_params[end]["CG_LP_IP_gap"],
            # neighborhood size
            neighborhood_size_mean_first = all_params[1]["ngroute_neighborhood_size"],
            neighborhood_size_mean_last = all_params[ind]["ngroute_neighborhood_size"],
            neighborhood_size_mean_last_SR3 = all_params[end]["ngroute_neighborhood_size"],
            # cuts
            implemented_SR3_cuts_count_total = sum(all_params_df.implemented_SR3_cuts_count),
            # Length metrics 
            mean_subpath_length = get(some_paths_metrics, "mean_subpath_length", missing),
            mean_path_length = get(some_paths_metrics, "mean_path_length", missing),
            mean_ps_length = get(some_paths_metrics, "mean_ps_length", missing),
        )
    ]
    if write_log
        CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
    end
    if write_results
        CSV.write("$(@__DIR__)/results/$(row_index).csv", all_params_df)    
    end
    println("Row $row_index processed: $(args_df[row_index, :])")
end



# simple test case to quickly compile 

begin
    test_args_df = DataFrame(CSV.File("$(@__DIR__)/test_args.csv"))
    for i in 1:nrow(test_args_df)
        run_instance(
            test_args_df, i, 100.0,
            ;
            write_log = false,
            write_results = false,
        )
    end
end

println("Compilation complete.")



args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

task_index = parse(Int, ARGS[1]) + 1
n_tasks = parse(Int, ARGS[2])

println("Processing rows: $(collect(task_index:n_tasks:size(args_df, 1)))")

for row_index in task_index:n_tasks:size(args_df, 1)
    sleep(rand() * 30)
    run_instance(
        args_df, row_index, 3600.0,           
        ;
        write_log = true,
        write_results = true,
    )
end