include("../../../src/arc_formulation.jl")
include("../../../src/utils.jl")

using StatsBase
using Suppressor
using CSV
using DataFrames

# simple test case to quickly compile
_, sample_data = generate_instance_pair(
    n_depots = 2, 
    n_customers = 9,
    n_charging = 2,
    charging_repeats = 1,
    n_vehicles = 3,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.6,
    T = 1350.0,
    seed = 1,
    B = 800.0,
    μ = 5.0,
    batch = 3,
)
sample_data = preprocess_arcs(sample_data, true, false)
sample_G = construct_graph(sample_data)
sample_T_range = 0:50.0:sample_data["T"]
sample_B_range = 0:50.0:sample_data["B"]

(
    sample_arc_ip_results, 
    sample_arc_ip_params,
) = arc_formulation(
    sample_data,
    true,
    ;
    time_limit = 300.0,
)
(
    sample_arc_lp_results, 
    sample_arc_lp_params,
) = arc_formulation(
    sample_data,
    true,
    ;
    integral = false,
    time_limit = 300.0,
)
println("Compilation complete.")

args_df = DataFrame(CSV.File("$(@__DIR__)/args.csv"))

task_index = parse(Int, ARGS[1]) + 1
n_tasks = parse(Int, ARGS[2])

println("Processing rows: $(collect(task_index:n_tasks:size(args_df, 1)))")

for row_index in task_index:n_tasks:size(args_df, 1)

    # Get paramters from args_df at row row_index
    n_depots = args_df[row_index, :n_depots]
    n_customers = args_df[row_index, :n_customers]
    n_charging = args_df[row_index, :n_charging]
    charging_repeats = args_df[row_index, :charging_repeats]
    n_vehicles = args_df[row_index, :n_vehicles]
    shrinkage_depots = args_df[row_index, :shrinkage_depots]
    shrinkage_charging = args_df[row_index, :shrinkage_charging]
    T = args_df[row_index, :T]
    T_step = args_df[row_index, :T_step]
    B = args_df[row_index, :B]
    B_step = args_df[row_index, :B_step]
    seed = args_df[row_index, :seed]
    μ = args_df[row_index, :μ]
    batch = args_df[row_index, :batch]
    time_limit = args_df[row_index, :time_limit]

    _, data = generate_instance_pair(
        n_depots = n_depots,
        n_customers = n_customers,
        n_charging = n_charging,
        charging_repeats = charging_repeats,
        n_vehicles = n_vehicles,
        shrinkage_depots = shrinkage_depots,
        shrinkage_charging = shrinkage_charging,
        T = T,
        seed = seed,
        B = B,
        μ = μ,
        batch = batch,
    )
    data = preprocess_arcs(data, true, false)
    G = construct_graph(data)
    T_range = 0:T_step:data["T"]
    B_range = 0:B_step:data["B"]

    (
        arc_ip_results, 
        arc_ip_params,
    ) = arc_formulation(
        data,
        true,
        ;
        integral = true,
        time_limit = time_limit,
    )

    (
        arc_lp_results, 
        arc_lp_params,
    ) = arc_formulation(
        data,
        true,
        ;
        integral = false,
        time_limit = time_limit,
    )

    records = [
        (
            n_depots = n_depots, 
            n_customers = n_customers, 
            n_charging = n_charging, 
            charging_repeats = charging_repeats, 
            n_vehicles = n_vehicles, 
            shrinkage_depots = shrinkage_depots, 
            shrinkage_charging = shrinkage_charging, 
            T = T, 
            T_step = T_step, 
            B = B, 
            B_step = B_step, 
            seed = seed, 
            μ = μ, 
            batch = batch, 
            time_limit = time_limit, 
            arc_ip_objective = (
                arc_ip_results["status"] in ["optimal", "time_limit_with_values"] ? 
                arc_ip_results["objective"] : missing
            ), 
            arc_ip_time_taken = arc_ip_params["time_taken"],
            arc_ip_solution_time_taken = arc_ip_params["solution_time_taken"],
            arc_ip_constraint_time_taken = arc_ip_params["constraint_time_taken"],
            arc_lp_objective = (
                arc_lp_results["status"] in ["optimal", "time_limit_with_values"] ? 
                arc_lp_results["objective"] : missing
            ), 
            arc_lp_time_taken = arc_lp_params["time_taken"],
            arc_lp_solution_time_taken = arc_lp_params["solution_time_taken"],
            arc_lp_constraint_time_taken = arc_lp_params["constraint_time_taken"],
        )
    ]
    CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
end