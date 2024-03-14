using Pkg
Pkg.activate("$(@__DIR__)/../../..")
println(Pkg.status())

include("$(@__DIR__)/../../../src/utils.jl")
include("$(@__DIR__)/../../../src/path_formulation.jl")

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
    write_results::Bool = false,
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
    elementary = args_df[row_index, :elementary]
    ngroute = args_df[row_index, :ngroute]

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
        charge_cost_heterogenous = false,
        ;
        data_dir = "../../../data/",
    )
    graph = generate_graph_from_data(data)

    run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        Env = GRB_ENV,
        method = method,
        charge_cost_heterogenous = false,
        elementary = elementary,
        ngroute = ngroute,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small",
        ngroute_neighborhood_charging_size = "small",
        verbose = true,
        use_adaptive_ngroute = false,
        use_SR3_cuts = false,
        use_lmSR3_cuts = false,
        time_limit = time_limit,
    );
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = run.value;

    local records
    if write_log
        if length(all_params) == 0
            return
        end
        if all_params[end]["errored"]
            println("Row index $(row_index) encountered errored model.")
            return
        end
        some_paths_metrics = compute_path_metrics(some_paths)
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
                elementary = elementary,
                ngroute = ngroute,
                ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
                ngroute_neighborhood_depots_size = "small",
                ngroute_neighborhood_charging_size = "small",
                use_adaptive_ngroute = false,
                use_SR3_cuts = false,
                use_lmSR3_cuts = false,
                # Time taken
                time_taken = CG_all_params[end]["time_taken"],
                sp_base_time_taken_total = CG_all_params[end]["sp_base_time_taken_total"],
                sp_full_time_taken_total = CG_all_params[end]["sp_full_time_taken_total"],
                sp_time_taken_total = CG_all_params[end]["sp_time_taken_total"],
                lp_relaxation_time_taken_total = CG_all_params[end]["lp_relaxation_time_taken_total"],
                sp_base_time_taken_mean = CG_all_params[end]["sp_base_time_taken_mean"],
                sp_full_time_taken_mean = CG_all_params[end]["sp_full_time_taken_mean"],
                sp_time_taken_mean = CG_all_params[end]["sp_time_taken_mean"],
                lp_relaxation_time_taken_mean = CG_all_params[end]["lp_relaxation_time_taken_mean"],
                # n_iterations
                counter = CG_all_params[end]["counter"],
                converged = CG_all_params[end]["converged"],
                time_limit_reached = CG_all_params[end]["time_limit_reached"],
                # Objective values and gaps
                LP_artificial = all_params[end]["CGLP_artificial"],
                LP_objective = all_params[end]["CGLP_objective"],
                IP_artificial = all_params[end]["CGIP_artificial"],
                IP_objective = all_params[end]["CGIP_objective"],
                LP_IP_gap = all_params[end]["CG_LP_IP_gap"],
                # Length metrics 
                mean_subpath_length = get(some_paths_metrics, "mean_subpath_length", missing),
                mean_path_length = get(some_paths_metrics, "mean_path_length", missing),
                mean_ps_length = get(some_paths_metrics, "mean_ps_length", missing),
            )
        ]
        CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
        println("$records")
    end
    if write_results
        mkpath("$(@__DIR__)/results/$(row_index)")
        savefig(plot_instance(data), "$(@__DIR__)/results/$(row_index)/plot.png")
        savefig(plot_path_solution(CGLP_all_results[end], data, graph, some_paths), "$(@__DIR__)/results/$(row_index)/plot_path_solution_LP.png")
        savefig(plot_path_solution(CGIP_all_results[end], data, graph, some_paths), "$(@__DIR__)/results/$(row_index)/plot_path_solution_IP.png")
    end
    println("Row $row_index processed. ") 
    println()
end



# simple test case to quickly compile 

begin
    test_args_df = DataFrame(CSV.File("$(@__DIR__)/test_args.csv"))
    # run_instance(test_args_df, 2, 3600.0, write_log = false)
    for i in 1:nrow(test_args_df)
        run_instance(test_args_df, i, 50.0, write_log = false)
    end
end

println("Compilation complete.")

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

task_index = parse(Int, ARGS[1]) + 1
n_tasks = parse(Int, ARGS[2])

println("Processing rows: $(collect(task_index:n_tasks:size(args_df, 1)))")

for row_index in task_index:n_tasks:size(args_df, 1)
    if !isfile("$(@__DIR__)/records/$row_index.csv")
        run_instance(args_df, row_index, 3600.0, write_log = true)
    end
end