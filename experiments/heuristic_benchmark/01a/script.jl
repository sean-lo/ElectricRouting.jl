using Pkg
Pkg.activate("$(@__DIR__)/../../..")
println(Pkg.status())

include("$(@__DIR__)/../../../src/utils.jl")
include("$(@__DIR__)/../../../src/path_formulation.jl")
include("$(@__DIR__)/../../../src/decomposition_heuristic.jl")

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
        time_limit = time_limit,
    );
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
        some_paths, model, z, SR3_constraints
    ) = optimal_run.value;

    some_paths_metrics = compute_path_metrics(some_paths)
    
    local heuristic_run

    try
        heuristic_run = @timed path_formulation_decomposition_heuristic(
            data, graph;
            elementary = false,
            ngroute = true,
            use_adaptive_ngroute = heuristic_use_adaptive_ngroute,
            use_SR3_cuts = heuristic_use_SR3_cuts,
            use_lmSR3_cuts = heuristic_use_lmSR3_cuts,
            max_SR3_cuts = max_SR3_cuts,
            time_limit = time_limit,
            time_heuristic_slack = 0.9,
        );
    catch e
        if isa(e, CGException)
            println("Row $row_index errored due to CGException.:\n $(args_df[row_index, :])")
            return
        else
            throw(e)
        end
    end
    (
        hCGIP_result, 
        heuristic_results,
        time_heuristic_slack,
    ) = heuristic_run.value;

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
            ngroute_alt = true,
            ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
            ngroute_neighborhood_depots_size = "small",
            ngroute_neighborhood_charging_size = ngroute_neighborhood_charging_size,
            use_adaptive_ngroute = use_adaptive_ngroute,
            use_SR3_cuts = use_SR3_cuts,
            max_SR3_cuts = max_SR3_cuts,
            use_lmSR3_cuts = use_lmSR3_cuts,
            heuristic_use_adaptive_ngroute = heuristic_use_adaptive_ngroute,
            heuristic_use_SR3_cuts = heuristic_use_SR3_cuts,
            heuristic_max_SR3_cuts = max_SR3_cuts,
            heuristic_use_lmSR3_cuts = heuristic_use_lmSR3_cuts,
            # Time taken
            time_limit = time_limit,
            h_time = heuristic_run.time,
            opt_time = optimal_run.time,
            # Objective values and gaps
            time_heuristic_slack = time_heuristic_slack,
            h_IP_objective = hCGIP_result["objective"],
            hc_IP_objective = heuristic_results["objective"],
            opt_IP_objective = CGIP_all_results[end]["objective"],
            opt_LP_objective = CGLP_all_results[end]["objective"],
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
        run_instance(
            test_args_df, i, 300.0,
            ;
            write_log = true,
        )
    end
end


println("Compilation complete.")

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

task_index = parse(Int, ARGS[1]) + 1
n_tasks = parse(Int, ARGS[2])

println("Processing rows: $(collect(task_index:n_tasks:size(args_df, 1)))")

for row_index in task_index:n_tasks:size(args_df, 1)
    if !isfile("$(@__DIR__)/records/$row_index.csv") & !(row_index in [427, 407])
        sleep(rand() * 30)
        run_instance(
            args_df, row_index, 3600.0,
            ;
            write_log = true,
        )
    end
end