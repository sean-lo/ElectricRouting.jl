include("../../../src/subpath_formulation.jl")
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
    sample_cg_results,
    sample_cg_params,
    sample_cg_printlist,
    sample_cg_subpaths,
) = subpath_formulation_column_generation_from_paths(
    sample_G,
    sample_data, 
    sample_T_range,
    sample_B_range,
    ;
    charging_in_subpath = true,
    verbose = true,
)
sample_cg_number_of_subpaths = sum(
    length(v) for v in values(sample_cg_subpaths)
)
sample_cg_subpath_costs = compute_subpath_costs(
    sample_data, 
    sample_cg_subpaths,
)
sample_cg_subpath_service = compute_subpath_service(
    sample_data, 
    sample_cg_subpaths,
)
(
    sample_cg_lpip_results,
    sample_cg_lpip_params,
) = subpath_formulation(
    sample_data,
    sample_cg_subpaths,
    sample_cg_subpath_costs,
    sample_cg_subpath_service,
    sample_T_range,
    sample_B_range,
    ;
    integral = true,
    charging_in_subpath = true,
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
        cg_results,
        cg_params,
        cg_printlist,
        cg_subpaths,
    ) = subpath_formulation_column_generation_from_paths(
        G,
        data, 
        T_range,
        B_range,
        ;
        charging_in_subpath = true,
        verbose = true,
        time_limit = time_limit,
    )
    cg_number_of_subpaths = sum(
        length(v) for v in values(cg_subpaths)
    )
    filename = "cg_lp_$(n_customers)_$(seed)_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))"
    open("$(@__DIR__)/logs/$filename.txt", "w") do io
        for message in cg_printlist
            write(io, message)
        end
    end
    cg_number_of_iterations = length(cg_params["sp_total_time_taken"])
    groupedbar(
        hcat(
            cg_params["lp_relaxation_time_taken"],
            cg_params["sp_total_time_taken"],
        ),
        group = repeat(
            ["LP relaxation solve time", "Subproblem solve time"], 
            inner = cg_number_of_iterations,
        ),
        bar_position = :stack,
        framestyle = :box,
        xlabel = "Iteration",
        xticks = 2:cg_number_of_iterations,
        ylabel = "Time (s)",
        title = """
        Time of LP relaxation and subproblem, with artificial starting subpaths
        ($(data["n_customers"]) customers, $(data["n_vehicles"]) vehicles, $(data["n_depots"]) depots, $(data["n_charging"]) charging stations)
        """,
        size = (800, 600),
    )
    savefig("$(@__DIR__)/plots/$filename.png")
    cg_subpath_costs = compute_subpath_costs(
        data, 
        cg_subpaths,
    )
    cg_subpath_service = compute_subpath_service(
        data, 
        cg_subpaths,
    )
    (
        cg_lpip_results,
        cg_lpip_params,
    ) = subpath_formulation(
        data,
        cg_subpaths,
        cg_subpath_costs,
        cg_subpath_service,
        T_range,
        B_range,
        ;
        integral = true,
        charging_in_subpath = true,
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
            cg_number_of_subpaths = cg_number_of_subpaths,
            cg_number_of_iterations = cg_number_of_iterations,
            cg_objective = cg_results["objective"], 
            cg_time_taken = cg_params["time_taken"],
            cg_mp_total_time_taken = sum(cg_params["lp_relaxation_time_taken"]),
            cg_sp_total_time_taken = sum(cg_params["sp_total_time_taken"]),
            cg_mp_mean_time_taken = mean(cg_params["lp_relaxation_time_taken"]),
            cg_sp_mean_time_taken = mean(cg_params["sp_total_time_taken"]),
            cg_mp_std_time_taken = std(cg_params["lp_relaxation_time_taken"]),
            cg_sp_std_time_taken = std(cg_params["sp_total_time_taken"]),
            cg_lpip_objective = cg_lpip_results["objective"],
            cg_lpip_time_taken = cg_lpip_params["time_taken"],
            cg_lpip_solution_time_taken = cg_lpip_params["solution_time_taken"],
            cg_lpip_constraint_time_taken = cg_lpip_params["constraint_time_taken"],
        )
    ]
    CSV.write("$(@__DIR__)/records/$(row_index).csv", DataFrame(records))
end