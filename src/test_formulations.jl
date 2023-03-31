include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using DataFrames
using Dates
using CSV, JLD2
using Test
using Plots
using CairoMakie

all_data = Dict(
    "xs" => Dict(),
    "s" => Dict(),
    "m" => Dict(),
    "l" => Dict(),
)
for ind in 1:5
    all_data["xs_B$ind"] = Dict()
    all_data["s_B$ind"] = Dict()
end
for ind in 1:4
    all_data["xs_TW$ind"] = Dict()
    all_data["s_TW$ind"] = Dict()
end

all_data["xs2"] = Dict()
all_data["s2"] = Dict()
all_data["m2"] = Dict()
all_data["l2"] = Dict()

params = [
    ("xs", 2, 9, 2, 3, 1500.0, 1050.0, 10, 7, 3, 0.4),
    ("xs_B1", 2, 9, 2, 3, 1500.0, 1100.0, 10, 7, 3, 0.4),
    ("xs_B2", 2, 9, 2, 3, 1500.0, 1150.0, 10, 7, 3, 0.4),
    ("xs_B3", 2, 9, 2, 3, 1500.0, 1200.0, 10, 7, 3, 0.4),
    ("xs_B4", 2, 9, 2, 3, 1500.0, 1250.0, 10, 7, 3, 0.4),
    ("xs_B5", 2, 9, 2, 3, 1500.0, 1300.0, 10, 7, 3, 0.4),
    ("xs_TW1", 2, 9, 2, 3, 1500.0, 1050.0, 10, 7, 3, 0.45),
    ("xs_TW2", 2, 9, 2, 3, 1500.0, 1050.0, 10, 7, 3, 0.5),
    ("xs_TW3", 2, 9, 2, 3, 1500.0, 1050.0, 10, 7, 3, 0.55),
    ("xs_TW4", 2, 9, 2, 3, 1500.0, 1050.0, 10, 7, 3, 0.6),
    ("s", 2, 12, 2, 3, 1900.0, 1550.0, 10, 7, 4, 0.4),
    ("s_B1", 2, 12, 2, 3, 1900.0, 1600.0, 10, 7, 4, 0.4),
    ("s_B2", 2, 12, 2, 3, 1900.0, 1650.0, 10, 7, 4, 0.4),
    ("s_B3", 2, 12, 2, 3, 1900.0, 1700.0, 10, 7, 4, 0.4),
    ("s_B4", 2, 12, 2, 3, 1900.0, 1750.0, 10, 7, 4, 0.4),
    ("s_B5", 2, 12, 2, 3, 1900.0, 1800.0, 10, 7, 4, 0.4),
    ("s_TW1", 2, 12, 2, 3, 1900.0, 1550.0, 10, 7, 4, 0.45),
    ("s_TW2", 2, 12, 2, 3, 1900.0, 1550.0, 10, 7, 4, 0.5),
    ("s_TW3", 2, 12, 2, 3, 1900.0, 1550.0, 10, 7, 4, 0.55),
    ("s_TW4", 2, 12, 2, 3, 1900.0, 1550.0, 10, 7, 4, 0.6),
    ("m", 2, 15, 2, 3, 2250.0, 1850.0, 10, 7, 5, 0.4),
    ("l", 2, 18, 2, 3, 2550.0, 2050.0, 10, 7, 6, 0.4),
    ("xs2", 2, 9, 2, 3, 1500.0, 750.0, 10, 7, 3, 0.4),
    ("s2", 2, 12, 2, 3, 2000.0, 1000.0, 10, 7, 4, 0.4),
    ("m2", 2, 15, 2, 3, 2400.0, 1200.0, 10, 7, 5, 0.4),
    ("l2", 2, 18, 2, 3, 2700.0, 1350.0, 10, 7, 6, 0.4),
]

for (
    size,
    n_depots, n_customers, n_charging, n_vehicles,
    T, B, 
    travel_cost_coeff, charge_cost_coeff, 
    batch, permissiveness,
) in params, seed in 1:10
    run_index = seed
    _, data = generate_instance_pair(
        n_depots = n_depots, 
        n_customers = n_customers,
        n_charging = n_charging,
        charging_repeats = 1,
        n_vehicles = n_vehicles,
        shrinkage_depots = 1.4,
        shrinkage_charging = 0.6,
        T = T,
        seed = seed,
        B = B,
        μ = 5.0,
        travel_cost_coeff = travel_cost_coeff,
        charge_cost_coeff = charge_cost_coeff,
        batch = batch,
        permissiveness = permissiveness,
    )
    data = preprocess_arcs(data, true, false)
    G = construct_graph(data)
    T_range = 0:50.0:data["T"]
    B_range = 0:50.0:data["B"]
    all_data[size][run_index] = Dict(
        "data" => data,
        "G" => G,
        "T_range" => T_range,
        "B_range" => B_range,
    )
end

function compare_formulations!(
    all_data,
    sizes,
    run_indexes,
    ;
    with_charging_cost::Bool = false,
    time_windows::Bool = true,
    arc::Bool = false,
    arc_sizes::Vector = sizes,
    arc_run_indexes::Vector = run_indexes,
    cg_charge_to_full_only::Bool = false,
    cg_with_heuristic::Bool = true,
    subpath_cg::Bool = false,
    subpath_cg_sizes::Vector = sizes,
    subpath_cg_run_indexes::Vector = run_indexes,
    subpath_cgi::Bool = false,
    subpath_cgi_sizes::Vector = sizes,
    subpath_cgi_run_indexes::Vector = run_indexes,
    subpath_cgi_no_time_windows_naive::Bool = true,
    subpath_enum::Bool = false,
    subpath_enum_sizes::Vector = sizes,
    subpath_enum_run_indexes::Vector = run_indexes,
    path_cg::Bool = false,
    path_cg_sizes::Vector = sizes,
    path_cg_run_indexes::Vector = run_indexes,
)
    datestr = Dates.format(Dates.now(), "yyyymmdd_HHMMSS")
    dir = "$(@__DIR__)/../logs/$datestr/"
    mkpath(dir)

    if arc
        for (size, run_index) in Iterators.product(arc_sizes, arc_run_indexes)
            dirname = joinpath(dir, "arc/$size/$run_index")
            mkpath(dirname)
            (
                all_data[size][run_index]["arc_ip_results"], 
                all_data[size][run_index]["arc_ip_params"],
            ) = arc_formulation(
                all_data[size][run_index]["data"],
                true,
                ;
                with_charging_cost = with_charging_cost,
                time_limit = 3600.0,
            )
            all_data[size][run_index]["arc_ip_printlist"] = arc_results_printout(
                all_data[size][run_index]["arc_ip_results"], 
                all_data[size][run_index]["arc_ip_params"],
                all_data[size][run_index]["data"],
                true,
            )
            open(joinpath(dirname, "log.txt"), "w") do io
                for message in all_data[size][run_index]["arc_ip_printlist"]
                    write(io, message)
                end
            end
            (
                all_data[size][run_index]["arc_lp_results"], 
                all_data[size][run_index]["arc_lp_params"],
            ) = arc_formulation(
                all_data[size][run_index]["data"],
                true,
                ;
                with_charging_cost = with_charging_cost,
                integral = false,
                time_limit = 120.0,
            )
        end
    end

    if subpath_cg
        for (size, run_index) in Iterators.product(subpath_cg_sizes, subpath_cg_run_indexes)
            dirname = joinpath(dir, "cg/$size/$run_index")
            mkpath(dirname)
            (
                all_data[size][run_index]["cg_results"],
                all_data[size][run_index]["cg_params"],
                all_data[size][run_index]["cg_printlist"],
                all_data[size][run_index]["cg_subpaths"],
            ) = subpath_formulation_column_generation_from_paths(
                all_data[size][run_index]["G"],
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["T_range"],
                all_data[size][run_index]["B_range"],
                ;
                charging_in_subpath = true,
                charge_to_full_only = cg_charge_to_full_only,
                time_windows = time_windows,
                with_charging_cost = with_charging_cost,
                with_heuristic = cg_with_heuristic,
                verbose = true,
            );
            all_data[size][run_index]["cg_number_of_subpaths"] = sum(
                length(v) for v in values(all_data[size][run_index]["cg_subpaths"])
            )
            open(joinpath(dirname, "log.txt"), "w") do io
                for message in all_data[size][run_index]["cg_printlist"]
                    write(io, message)
                end
            end
            groupedbar(
                hcat(
                    all_data[size][run_index]["cg_params"]["lp_relaxation_time_taken"],
                    all_data[size][run_index]["cg_params"]["sp_total_time_taken"],
                ),
                group = repeat(
                    ["LP relaxation solve time", "Subproblem solve time"], 
                    inner = length(all_data[size][run_index]["cg_params"]["sp_total_time_taken"]),
                ),
                bar_position = :stack,
                framestyle = :box,
                xlabel = "Iteration",
                xticks = 2:length(all_data[size][run_index]["cg_params"]["sp_total_time_taken"]),
                ylabel = "Time (s)",
                title = """
                Time of LP relaxation and subproblem, with artificial starting subpaths
                ($(all_data[size][run_index]["data"]["n_customers"]) customers, $(all_data[size][run_index]["data"]["n_vehicles"]) vehicles, $(all_data[size][run_index]["data"]["n_depots"]) depots, $(all_data[size][run_index]["data"]["n_charging"]) charging stations)
                """,
                size = (800, 600),
            )
            savefig(joinpath(dirname, "barplot.png"))
            all_data[size][run_index]["cg_subpath_costs"] = compute_subpath_costs(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["cg_subpaths"],
                ;
                with_charging_cost = with_charging_cost,
                time_windows = time_windows,
            )
            all_data[size][run_index]["cg_subpath_service"] = compute_subpath_service(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["cg_subpaths"],
            )
            (
                all_data[size][run_index]["cg_lpip_results"],
                all_data[size][run_index]["cg_lpip_params"],
            ) = subpath_formulation(
                all_data[size][run_index]["data"],
                all_data[size][run_index]["cg_subpaths"],
                all_data[size][run_index]["cg_subpath_costs"],
                all_data[size][run_index]["cg_subpath_service"],
                all_data[size][run_index]["T_range"],
                all_data[size][run_index]["B_range"],
                ;
                integral = true,
                charging_in_subpath = true,
            )
            delete!(all_data[size][run_index], "cg_subpaths")
            delete!(all_data[size][run_index], "cg_subpath_costs")
            delete!(all_data[size][run_index], "cg_subpath_service")
            delete!(all_data[size][run_index], "cg_printlist")
            delete!(all_data[size][run_index]["cg_results"], "z")
            delete!(all_data[size][run_index]["cg_lpip_results"], "λ")
            delete!(all_data[size][run_index]["cg_lpip_results"], "z")
        end
    end

    if subpath_cgi
        for (size, run_index) in Iterators.product(subpath_cgi_sizes, subpath_cgi_run_indexes)
            dirname = joinpath(dir, "cgi/$size/$run_index")
            mkpath(dirname)

            data_ntw = deepcopy(all_data[size][run_index]["data"])
            data_ntw["α"] .= 0.0
            data_ntw["β"] .= data_ntw["T"]

            if time_windows
                (
                    all_data[size][run_index]["cgi_results"],
                    all_data[size][run_index]["cgi_params"],
                    all_data[size][run_index]["cgi_printlist"],
                    all_data[size][run_index]["cgi_subpaths"],
                ) = subpath_formulation_column_generation_integrated_from_paths(
                    all_data[size][run_index]["G"],
                    all_data[size][run_index]["data"], 
                    all_data[size][run_index]["T_range"],
                    all_data[size][run_index]["B_range"],
                    ;
                    charge_to_full_only = cg_charge_to_full_only,
                    time_windows = true,
                    with_charging_cost = with_charging_cost,
                    with_heuristic = cg_with_heuristic,
                    verbose = true,
                )
            else
                if subpath_cgi_no_time_windows_naive
                    
                    (
                        all_data[size][run_index]["cgi_results"],
                        all_data[size][run_index]["cgi_params"],
                        all_data[size][run_index]["cgi_printlist"],
                        all_data[size][run_index]["cgi_subpaths"],
                    ) = subpath_formulation_column_generation_integrated_from_paths(
                        all_data[size][run_index]["G"],
                        data_ntw, 
                        all_data[size][run_index]["T_range"],
                        all_data[size][run_index]["B_range"],
                        ;
                        charge_to_full_only = cg_charge_to_full_only,
                        time_windows = true,
                        with_charging_cost = with_charging_cost,
                        verbose = true,
                    )
                else
                    (
                        all_data[size][run_index]["cgi_results"],
                        all_data[size][run_index]["cgi_params"],
                        all_data[size][run_index]["cgi_printlist"],
                        all_data[size][run_index]["cgi_subpaths"],
                    ) = subpath_formulation_column_generation_integrated_from_paths(
                        all_data[size][run_index]["G"],
                        all_data[size][run_index]["data"], 
                        all_data[size][run_index]["T_range"],
                        all_data[size][run_index]["B_range"],
                        ;
                        charge_to_full_only = cg_charge_to_full_only,
                        time_windows = false,
                        with_charging_cost = with_charging_cost,
                        verbose = true,
                    )
                end
            end
            all_data[size][run_index]["cgi_number_of_subpaths"] = sum(
                length(v) for v in values(all_data[size][run_index]["cgi_subpaths"])
            )
            open(joinpath(dirname, "log.txt"), "w") do io
                for message in all_data[size][run_index]["cgi_printlist"]
                    write(io, message)
                end
            end
            groupedbar(
                hcat(
                    all_data[size][run_index]["cgi_params"]["lp_relaxation_time_taken"],
                    all_data[size][run_index]["cgi_params"]["sp_total_time_taken"],
                ),
                group = repeat(
                    ["LP relaxation solve time", "Subproblem solve time"], 
                    inner = length(all_data[size][run_index]["cgi_params"]["sp_total_time_taken"]),
                ),
                bar_position = :stack,
                framestyle = :box,
                xlabel = "Iteration",
                xticks = 2:length(all_data[size][run_index]["cgi_params"]["sp_total_time_taken"]),
                ylabel = "Time (s)",
                title = """
                Time of LP relaxation and subproblem, with artificial starting subpaths
                ($(all_data[size][run_index]["data"]["n_customers"]) customers, $(all_data[size][run_index]["data"]["n_vehicles"]) vehicles, $(all_data[size][run_index]["data"]["n_depots"]) depots, $(all_data[size][run_index]["data"]["n_charging"]) charging stations)
                """,
                size = (800, 600),
            )
            savefig(joinpath(dirname, "barplot.png"))
            if time_windows || !subpath_cgi_no_time_windows_naive
                    all_data[size][run_index]["cgi_subpath_costs"] = compute_subpath_costs(
                    all_data[size][run_index]["data"], 
                    all_data[size][run_index]["cgi_subpaths"],
                    ;
                    with_charging_cost = with_charging_cost,
                    time_windows = time_windows,
                )
                all_data[size][run_index]["cgi_subpath_service"] = compute_subpath_service(
                    all_data[size][run_index]["data"], 
                    all_data[size][run_index]["cgi_subpaths"],
                )
                (
                    all_data[size][run_index]["cgi_lpip_results"],
                    all_data[size][run_index]["cgi_lpip_params"],
                ) = subpath_formulation(
                    all_data[size][run_index]["data"],
                    all_data[size][run_index]["cgi_subpaths"],
                    all_data[size][run_index]["cgi_subpath_costs"],
                    all_data[size][run_index]["cgi_subpath_service"],
                    all_data[size][run_index]["T_range"],
                    all_data[size][run_index]["B_range"],
                    ;
                    integral = true,
                    charging_in_subpath = true,
                )
            else
                all_data[size][run_index]["cgi_subpath_costs"] = compute_subpath_costs(
                    data_ntw, 
                    all_data[size][run_index]["cgi_subpaths"],
                    ;
                    with_charging_cost = with_charging_cost,
                    time_windows = time_windows,
                )
                all_data[size][run_index]["cgi_subpath_service"] = compute_subpath_service(
                    data_ntw, 
                    all_data[size][run_index]["cgi_subpaths"],
                )
                (
                    all_data[size][run_index]["cgi_lpip_results"],
                    all_data[size][run_index]["cgi_lpip_params"],
                ) = subpath_formulation(
                    data_ntw,
                    all_data[size][run_index]["cgi_subpaths"],
                    all_data[size][run_index]["cgi_subpath_costs"],
                    all_data[size][run_index]["cgi_subpath_service"],
                    all_data[size][run_index]["T_range"],
                    all_data[size][run_index]["B_range"],
                    ;
                    integral = true,
                    charging_in_subpath = true,
                )
            end

            delete!(all_data[size][run_index], "cgi_subpaths")
            delete!(all_data[size][run_index], "cgi_subpath_costs")
            delete!(all_data[size][run_index], "cgi_subpath_service")
            delete!(all_data[size][run_index], "cgi_printlist")
            delete!(all_data[size][run_index]["cgi_results"], "z")
            delete!(all_data[size][run_index]["cgi_lpip_results"], "λ")
            delete!(all_data[size][run_index]["cgi_lpip_results"], "z")
        end
    end

    if subpath_enum
        for (size, run_index) in Iterators.product(subpath_enum_sizes, subpath_enum_run_indexes)
            (
                all_data[size][run_index]["all_subpaths"],
                all_data[size][run_index]["enumerate_all_subpaths_time_taken"],    
            ) = enumerate_all_subpaths(
                all_data[size][run_index]["G"], 
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["T_range"], 
                all_data[size][run_index]["B_range"]; 
                charging_in_subpath = true,
                charge_to_full_only = cg_charge_to_full_only,
                time_windows = time_windows,
            )
            all_data[size][run_index]["enum_number_of_subpaths"] = sum(
                length(v) for v in values(all_data[size][run_index]["all_subpaths"])
            )
            all_data[size][run_index]["all_subpath_costs"] = compute_subpath_costs(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["all_subpaths"],
                ;
                with_charging_cost = with_charging_cost,
                time_windows = time_windows,
            )
            all_data[size][run_index]["all_subpath_service"] = compute_subpath_service(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["all_subpaths"],
            )
            (
                all_data[size][run_index]["enum_ip_results"],
                all_data[size][run_index]["enum_ip_params"],
            ) = subpath_formulation(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["all_subpaths"],
                all_data[size][run_index]["all_subpath_costs"], 
                all_data[size][run_index]["all_subpath_service"],
                all_data[size][run_index]["T_range"],
                all_data[size][run_index]["B_range"],
                ;
                charging_in_subpath = true,
                integral = true,
            )
            (
                all_data[size][run_index]["enum_lp_results"],
                all_data[size][run_index]["enum_lp_params"],
            ) = subpath_formulation(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["all_subpaths"],
                all_data[size][run_index]["all_subpath_costs"], 
                all_data[size][run_index]["all_subpath_service"],
                all_data[size][run_index]["T_range"],
                all_data[size][run_index]["B_range"],
                ;
                charging_in_subpath = true,
                integral = false,
            )
            println("Number of starting states:\t$(length(all_data[size][run_index]["all_subpaths"],))")
            println("Number of subpaths:\t\t$(all_data[size][run_index]["enum_number_of_subpaths"])")
            @printf("Time taken:\t\t\t%7.1f s\n", all_data[size][run_index]["enumerate_all_subpaths_time_taken"])
            @printf("Objective: \t\t\t%6.1f\n", all_data[size][run_index]["enum_ip_results"]["objective"])
            @printf("LP Objective:\t\t\t%6.1f\n", all_data[size][run_index]["enum_lp_results"]["objective"])
            delete!(all_data[size][run_index], "all_subpaths")
            delete!(all_data[size][run_index], "all_subpath_costs")
            delete!(all_data[size][run_index], "all_subpath_service")
            delete!(all_data[size][run_index]["enum_ip_results"], "λ")
            delete!(all_data[size][run_index]["enum_ip_results"], "z")
            delete!(all_data[size][run_index]["enum_lp_results"], "λ")
            delete!(all_data[size][run_index]["enum_lp_results"], "z")
        end
    end

    if path_cg
        for (size, run_index) in Iterators.product(path_cg_sizes, path_cg_run_indexes)
            dirname = joinpath(dir, "path_cg/$size/$run_index")
            mkpath(dirname)
            (
                all_data[size][run_index]["path_cg_results"],
                all_data[size][run_index]["path_cg_params"],
                all_data[size][run_index]["path_cg_printlist"],
                all_data[size][run_index]["path_cg_paths"],
            ) = path_formulation_column_generation(
                all_data[size][run_index]["G"],
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["T_range"],
                all_data[size][run_index]["B_range"],
                ;
                charge_to_full_only = cg_charge_to_full_only,
                with_charging_cost = with_charging_cost,
                time_windows = time_windows,
                with_heuristic = cg_with_heuristic,
                verbose = true,
            );
            all_data[size][run_index]["path_cg_number_of_paths"] = sum(
                length(v) for v in values(all_data[size][run_index]["path_cg_paths"])
            )
            open(joinpath(dirname, "log.txt"), "w") do io
                for message in all_data[size][run_index]["path_cg_printlist"]
                    write(io, message)
                end
            end
            groupedbar(
                hcat(
                    all_data[size][run_index]["path_cg_params"]["lp_relaxation_time_taken"],
                    all_data[size][run_index]["path_cg_params"]["sp_total_time_taken"],
                ),
                group = repeat(
                    ["LP relaxation solve time", "Subproblem solve time"], 
                    inner = length(all_data[size][run_index]["path_cg_params"]["sp_total_time_taken"]),
                ),
                bar_position = :stack,
                framestyle = :box,
                xlabel = "Iteration",
                xticks = 2:length(all_data[size][run_index]["path_cg_params"]["sp_total_time_taken"]),
                ylabel = "Time (s)",
                title = """
                Time of LP relaxation and subproblem, with artificial starting paths
                ($(all_data[size][run_index]["data"]["n_customers"]) customers, $(all_data[size][run_index]["data"]["n_vehicles"]) vehicles, $(all_data[size][run_index]["data"]["n_depots"]) depots, $(all_data[size][run_index]["data"]["n_charging"]) charging stations)
                """,
                size = (800, 600),
            )
            savefig(joinpath(dirname, "barplot.png"))
            all_data[size][run_index]["path_cg_path_costs"] = compute_path_costs(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["path_cg_paths"],
                ;
                with_charging_cost = with_charging_cost,
                time_windows = time_windows,
            )
            all_data[size][run_index]["path_cg_path_service"] = compute_path_services(
                all_data[size][run_index]["data"], 
                all_data[size][run_index]["path_cg_paths"],
            )
            (
                all_data[size][run_index]["path_cg_lpip_results"],
                all_data[size][run_index]["path_cg_lpip_params"],
            ) = path_formulation(
                all_data[size][run_index]["data"],
                all_data[size][run_index]["path_cg_paths"],
                all_data[size][run_index]["path_cg_path_costs"],
                all_data[size][run_index]["path_cg_path_service"],
                all_data[size][run_index]["T_range"],
                all_data[size][run_index]["B_range"],
                ;
                integral = true,
            )
            delete!(all_data[size][run_index], "path_cg_paths")
            delete!(all_data[size][run_index], "path_cg_path_costs")
            delete!(all_data[size][run_index], "path_cg_path_service")
            delete!(all_data[size][run_index], "path_cg_printlist")
            delete!(all_data[size][run_index]["path_cg_results"], "y")
            delete!(all_data[size][run_index]["path_cg_lpip_results"], "y")
        end
    end

    metrics = [
        :size,
        :run_index,
        :arc_ip_objective, 
        :arc_ip_time_taken,
        :arc_ip_solution_time_taken,
        :arc_ip_constraint_time_taken,
        :arc_lp_objective,
        :arc_lp_time_taken,
        :arc_lp_solution_time_taken,
        :arc_lp_constraint_time_taken,
        :cg_number_of_subpaths,
        :cg_objective, 
        :cg_number_of_iterations,
        :cg_time_taken,
        :cg_mp_total_time_taken,
        :cg_mp_mean_time_taken,
        :cg_mp_total_constraint_time_taken,
        :cg_mp_mean_constraint_time_taken,
        :cg_mp_total_solution_time_taken,
        :cg_mp_mean_solution_time_taken,
        :cg_sp_total_time_taken,
        :cg_sp_mean_time_taken,   
        :cg_lpip_objective,
        :cg_lpip_time_taken,
        :cg_lpip_solution_time_taken,
        :cg_lpip_constraint_time_taken,
        :cgi_number_of_subpaths,
        :cgi_objective, 
        :cgi_number_of_iterations,
        :cgi_time_taken,
        :cgi_mp_total_time_taken,
        :cgi_mp_mean_time_taken,
        :cgi_mp_total_constraint_time_taken,
        :cgi_mp_mean_constraint_time_taken,
        :cgi_mp_total_solution_time_taken,
        :cgi_mp_mean_solution_time_taken,
        :cgi_sp_total_time_taken,
        :cgi_sp_mean_time_taken,   
        :cgi_lpip_objective,
        :cgi_lpip_time_taken,
        :cgi_lpip_solution_time_taken,
        :cgi_lpip_constraint_time_taken,
        :path_cg_number_of_paths,
        :path_cg_objective, 
        :path_cg_number_of_iterations,
        :path_cg_time_taken,
        :path_cg_mp_total_time_taken,
        :path_cg_mp_mean_time_taken,
        :path_cg_mp_total_constraint_time_taken,
        :path_cg_mp_mean_constraint_time_taken,
        :path_cg_mp_total_solution_time_taken,
        :path_cg_mp_mean_solution_time_taken,
        :path_cg_sp_total_time_taken,
        :path_cg_sp_mean_time_taken,   
        :path_cg_lpip_objective,
        :path_cg_lpip_time_taken,
        :path_cg_lpip_solution_time_taken,
        :path_cg_lpip_constraint_time_taken,
        :enum_number_of_subpaths,
        :enumerate_all_subpaths_time_taken,
        :enum_ip_objective, 
        :enum_ip_time_taken,
        :enum_ip_solution_time_taken,
        :enum_ip_constraint_time_taken,
        :enum_lp_objective,
        :enum_lp_time_taken,
        :enum_lp_solution_time_taken,
        :enum_lp_constraint_time_taken,
    ]
    all_metrics_r = []
    for (size, run_index) in Iterators.product(sizes, run_indexes)
        d = Dict(:size => size, :run_index => run_index)
        if "arc_ip_results" in keys(all_data[size][run_index])
            d[:arc_ip_objective] = all_data[size][run_index]["arc_ip_results"]["objective"]
        else
            d[:arc_ip_objective] = missing
        end
        if "arc_lp_results" in keys(all_data[size][run_index])
            d[:arc_lp_objective] = all_data[size][run_index]["arc_lp_results"]["objective"]
        else
            d[:arc_lp_objective] = missing
        end
        if "cg_results" in keys(all_data[size][run_index])
            d[:cg_objective] = all_data[size][run_index]["cg_results"]["objective"]
        else
            d[:cg_objective] = missing
        end
        if "cg_lpip_results" in keys(all_data[size][run_index])
            d[:cg_lpip_objective] = all_data[size][run_index]["cg_lpip_results"]["objective"]
        else
            d[:cg_lpip_objective] = missing
        end
        if "cgi_results" in keys(all_data[size][run_index])
            d[:cgi_objective] = all_data[size][run_index]["cgi_results"]["objective"]
        else
            d[:cgi_objective] = missing
        end
        if "cgi_lpip_results" in keys(all_data[size][run_index])
            d[:cgi_lpip_objective] = all_data[size][run_index]["cgi_lpip_results"]["objective"]
        else
            d[:cgi_lpip_objective] = missing
        end
        if "path_cg_results" in keys(all_data[size][run_index])
            d[:path_cg_objective] = all_data[size][run_index]["path_cg_results"]["objective"]
        else
            d[:path_cg_objective] = missing
        end
        if "path_cg_lpip_results" in keys(all_data[size][run_index])
            d[:path_cg_lpip_objective] = all_data[size][run_index]["path_cg_lpip_results"]["objective"]
        else
            d[:path_cg_lpip_objective] = missing
        end
        if "enum_ip_results" in keys(all_data[size][run_index])
            d[:enum_ip_objective] = all_data[size][run_index]["enum_ip_results"]["objective"]
        else
            d[:enum_ip_objective] = missing
        end
        if "enum_lp_results" in keys(all_data[size][run_index])
            d[:enum_lp_objective] = all_data[size][run_index]["enum_lp_results"]["objective"]
        else
            d[:enum_lp_objective] = missing
        end
        if "arc_ip_params" in keys(all_data[size][run_index])
            d[:arc_ip_time_taken] = all_data[size][run_index]["arc_ip_params"]["time_taken"]
            d[:arc_ip_solution_time_taken] = all_data[size][run_index]["arc_ip_params"]["solution_time_taken"]
            d[:arc_ip_constraint_time_taken] = all_data[size][run_index]["arc_ip_params"]["constraint_time_taken"]
        else
            d[:arc_ip_time_taken] = missing
            d[:arc_ip_solution_time_taken] = missing
            d[:arc_ip_constraint_time_taken] = missing
        end
        if "arc_lp_params" in keys(all_data[size][run_index])
            d[:arc_lp_time_taken] = all_data[size][run_index]["arc_lp_params"]["time_taken"]
            d[:arc_lp_solution_time_taken] = all_data[size][run_index]["arc_lp_params"]["solution_time_taken"]
            d[:arc_lp_constraint_time_taken] = all_data[size][run_index]["arc_lp_params"]["constraint_time_taken"]
        else
            d[:arc_lp_time_taken] = missing
            d[:arc_lp_solution_time_taken] = missing
            d[:arc_lp_constraint_time_taken] = missing
        end
        if "cg_params" in keys(all_data[size][run_index])
            d[:cg_number_of_iterations] = length(all_data[size][run_index]["cg_params"]["number_of_subpaths"])
            d[:cg_number_of_subpaths] = all_data[size][run_index]["cg_params"]["number_of_subpaths"][end]
            d[:cg_time_taken] = all_data[size][run_index]["cg_params"]["time_taken"]
            d[:cg_mp_total_time_taken] = sum(all_data[size][run_index]["cg_params"]["lp_relaxation_time_taken"])
            d[:cg_mp_total_constraint_time_taken] = sum(all_data[size][run_index]["cg_params"]["lp_relaxation_constraint_time_taken"])
            d[:cg_mp_total_solution_time_taken] = sum(all_data[size][run_index]["cg_params"]["lp_relaxation_solution_time_taken"])
            d[:cg_sp_total_time_taken] = sum(all_data[size][run_index]["cg_params"]["sp_total_time_taken"])
            d[:cg_mp_mean_time_taken] = mean(all_data[size][run_index]["cg_params"]["lp_relaxation_time_taken"])
            d[:cg_mp_mean_constraint_time_taken] = mean(all_data[size][run_index]["cg_params"]["lp_relaxation_constraint_time_taken"])
            d[:cg_mp_mean_solution_time_taken] = mean(all_data[size][run_index]["cg_params"]["lp_relaxation_solution_time_taken"])
            d[:cg_sp_mean_time_taken] = mean(all_data[size][run_index]["cg_params"]["sp_total_time_taken"])
        else
            d[:cg_number_of_iterations] = missing
            d[:cg_number_of_subpaths] = missing
            d[:cg_time_taken] = missing
            d[:cg_mp_total_time_taken] = missing
            d[:cg_mp_total_constraint_time_taken] = missing
            d[:cg_mp_total_solution_time_taken] = missing
            d[:cg_sp_total_time_taken] = missing
            d[:cg_mp_mean_time_taken] = missing
            d[:cg_mp_mean_constraint_time_taken] = missing
            d[:cg_mp_mean_solution_time_taken] = missing
            d[:cg_sp_mean_time_taken] = missing
        end
        if "cg_lpip_params" in keys(all_data[size][run_index])
            d[:cg_lpip_time_taken] = all_data[size][run_index]["cg_lpip_params"]["time_taken"]
            d[:cg_lpip_solution_time_taken] = all_data[size][run_index]["cg_lpip_params"]["solution_time_taken"]
            d[:cg_lpip_constraint_time_taken] = all_data[size][run_index]["cg_lpip_params"]["constraint_time_taken"]
        else
            d[:cg_lpip_time_taken] = missing
            d[:cg_lpip_solution_time_taken] = missing
            d[:cg_lpip_constraint_time_taken] = missing
        end
        if "cgi_params" in keys(all_data[size][run_index])
            d[:cgi_number_of_iterations] = length(all_data[size][run_index]["cgi_params"]["number_of_subpaths"])
            d[:cgi_number_of_subpaths] = all_data[size][run_index]["cgi_params"]["number_of_subpaths"][end]
            d[:cgi_time_taken] = all_data[size][run_index]["cgi_params"]["time_taken"]
            d[:cgi_mp_total_time_taken] = sum(all_data[size][run_index]["cgi_params"]["lp_relaxation_time_taken"])
            d[:cgi_mp_total_constraint_time_taken] = sum(all_data[size][run_index]["cgi_params"]["lp_relaxation_constraint_time_taken"])
            d[:cgi_mp_total_solution_time_taken] = sum(all_data[size][run_index]["cgi_params"]["lp_relaxation_solution_time_taken"])
            d[:cgi_sp_total_time_taken] = sum(all_data[size][run_index]["cgi_params"]["sp_total_time_taken"])
            d[:cgi_mp_mean_time_taken] = mean(all_data[size][run_index]["cgi_params"]["lp_relaxation_time_taken"])
            d[:cgi_mp_mean_constraint_time_taken] = mean(all_data[size][run_index]["cgi_params"]["lp_relaxation_constraint_time_taken"])
            d[:cgi_mp_mean_solution_time_taken] = mean(all_data[size][run_index]["cgi_params"]["lp_relaxation_solution_time_taken"])
            d[:cgi_sp_mean_time_taken] = mean(all_data[size][run_index]["cgi_params"]["sp_total_time_taken"])
        else
            d[:cgi_number_of_iterations] = missing
            d[:cgi_number_of_subpaths] = missing
            d[:cgi_time_taken] = missing
            d[:cgi_mp_total_time_taken] = missing
            d[:cgi_mp_total_constraint_time_taken] = missing
            d[:cgi_mp_total_solution_time_taken] = missing
            d[:cgi_sp_total_time_taken] = missing
            d[:cgi_mp_mean_time_taken] = missing
            d[:cgi_mp_mean_constraint_time_taken] = missing
            d[:cgi_mp_mean_solution_time_taken] = missing
            d[:cgi_sp_mean_time_taken] = missing
        end
        if "cgi_lpip_params" in keys(all_data[size][run_index])
            d[:cgi_lpip_time_taken] = all_data[size][run_index]["cgi_lpip_params"]["time_taken"]
            d[:cgi_lpip_solution_time_taken] = all_data[size][run_index]["cgi_lpip_params"]["solution_time_taken"]
            d[:cgi_lpip_constraint_time_taken] = all_data[size][run_index]["cgi_lpip_params"]["constraint_time_taken"]
        else
            d[:cgi_lpip_time_taken] = missing
            d[:cgi_lpip_solution_time_taken] = missing
            d[:cgi_lpip_constraint_time_taken] = missing
        end
        if "path_cg_params" in keys(all_data[size][run_index])
            d[:path_cg_number_of_iterations] = length(all_data[size][run_index]["path_cg_params"]["number_of_paths"])
            d[:path_cg_number_of_paths] = all_data[size][run_index]["path_cg_params"]["number_of_paths"][end]
            d[:path_cg_time_taken] = all_data[size][run_index]["path_cg_params"]["time_taken"]
            d[:path_cg_mp_total_time_taken] = sum(all_data[size][run_index]["path_cg_params"]["lp_relaxation_time_taken"])
            d[:path_cg_mp_total_constraint_time_taken] = sum(all_data[size][run_index]["path_cg_params"]["lp_relaxation_constraint_time_taken"])
            d[:path_cg_mp_total_solution_time_taken] = sum(all_data[size][run_index]["path_cg_params"]["lp_relaxation_solution_time_taken"])
            d[:path_cg_sp_total_time_taken] = sum(all_data[size][run_index]["path_cg_params"]["sp_total_time_taken"])
            d[:path_cg_mp_mean_time_taken] = mean(all_data[size][run_index]["path_cg_params"]["lp_relaxation_time_taken"])
            d[:path_cg_mp_mean_constraint_time_taken] = mean(all_data[size][run_index]["path_cg_params"]["lp_relaxation_constraint_time_taken"])
            d[:path_cg_mp_mean_solution_time_taken] = mean(all_data[size][run_index]["path_cg_params"]["lp_relaxation_solution_time_taken"])
            d[:path_cg_sp_mean_time_taken] = mean(all_data[size][run_index]["path_cg_params"]["sp_total_time_taken"])
        else
            d[:path_cg_number_of_iterations] = missing
            d[:path_cg_number_of_paths] = missing
            d[:path_cg_time_taken] = missing
            d[:path_cg_mp_total_time_taken] = missing
            d[:path_cg_mp_total_constraint_time_taken] = missing
            d[:path_cg_mp_total_solution_time_taken] = missing
            d[:path_cg_sp_total_time_taken] = missing
            d[:path_cg_mp_mean_time_taken] = missing
            d[:path_cg_mp_mean_constraint_time_taken] = missing
            d[:path_cg_mp_mean_solution_time_taken] = missing
            d[:path_cg_sp_mean_time_taken] = missing
        end
        if "path_cg_lpip_params" in keys(all_data[size][run_index])
            d[:path_cg_lpip_time_taken] = all_data[size][run_index]["path_cg_lpip_params"]["time_taken"]
            d[:path_cg_lpip_solution_time_taken] = all_data[size][run_index]["path_cg_lpip_params"]["solution_time_taken"]
            d[:path_cg_lpip_constraint_time_taken] = all_data[size][run_index]["path_cg_lpip_params"]["constraint_time_taken"]
        else
            d[:path_cg_lpip_time_taken] = missing
            d[:path_cg_lpip_solution_time_taken] = missing
            d[:path_cg_lpip_constraint_time_taken] = missing
        end
        if "all_subpaths" in keys(all_data[size][run_index])
            d[:enum_number_of_subpaths] = sum(length(v) for v in values(all_data[size][run_index]["all_subpaths"]))
        else
            d[:enum_number_of_subpaths] = missing
        end
        d[:enumerate_all_subpaths_time_taken] = get(all_data[size][run_index], "enumerate_all_subpaths_time_taken", missing)
        if "enum_ip_params" in keys(all_data[size][run_index])
            d[:enum_ip_time_taken] = all_data[size][run_index]["enum_ip_params"]["time_taken"]
            d[:enum_ip_solution_time_taken] = all_data[size][run_index]["enum_ip_params"]["solution_time_taken"]
            d[:enum_ip_constraint_time_taken] = all_data[size][run_index]["enum_ip_params"]["constraint_time_taken"]
        else
            d[:enum_ip_time_taken] = missing
            d[:enum_ip_solution_time_taken] = missing
            d[:enum_ip_constraint_time_taken] = missing
        end
        if "enum_lp_params" in keys(all_data[size][run_index])
            d[:enum_lp_time_taken] = all_data[size][run_index]["enum_lp_params"]["time_taken"]
            d[:enum_lp_solution_time_taken] = all_data[size][run_index]["enum_lp_params"]["solution_time_taken"]
            d[:enum_lp_constraint_time_taken] = all_data[size][run_index]["enum_lp_params"]["constraint_time_taken"]
        else
            d[:enum_lp_time_taken] = missing
            d[:enum_lp_solution_time_taken] = missing
            d[:enum_lp_constraint_time_taken] = missing
        end
        push!(all_metrics_r, NamedTuple(d))
    end
    all_metrics_df = DataFrame(all_metrics_r)[!,metrics]
    all_metrics_df[!, :cg_objective_integral] = (
        round.(all_metrics_df[!, :cg_objective], digits = 2)
        .== round.(all_metrics_df[!, :cg_lpip_objective], digits = 2)
    )
    CSV.write(joinpath(dir, "all_metrics.csv"), all_metrics_df)
    return all_metrics_df
end

compare_formulations!(
    all_data, ["xs"], [1], 
    with_charging_cost = true,
    cg_charge_to_full_only = false,
    time_windows = true,
    arc = true, 
    subpath_cg = true, 
    subpath_cgi = true, 
    path_cg = true,
)
compare_formulations!(
    all_data, ["xs"], [1], 
    with_charging_cost = true,
    cg_charge_to_full_only = false,
    time_windows = false,
    subpath_cgi = true, 
    subpath_cgi_no_time_windows_naive = true,
)
compare_formulations!(
    all_data, ["xs"], [1], 
    with_charging_cost = true,
    cg_charge_to_full_only = false,
    time_windows = false,
    subpath_cgi = true, 
    subpath_cgi_no_time_windows_naive = false,
)
compare_formulations!(
    all_data, ["xs"], [1], 
    with_charging_cost = true,
    cg_charge_to_full_only = false,
    cg_with_heuristic = false,
    time_windows = true,
    arc = true, 
    subpath_cg = true, 
    subpath_cgi = true, 
    path_cg = true,
)
compare_formulations!(
    all_data, ["xs"], [1], 
    with_charging_cost = true,
    cg_charge_to_full_only = false,    
    cg_with_heuristic = false,
    time_windows = false,
    subpath_cgi = true, 
    subpath_cgi_no_time_windows_naive = true,
)
compare_formulations!(
    all_data, ["xs"], [1], 
    with_charging_cost = true,
    cg_charge_to_full_only = false,
    cg_with_heuristic = false,
    time_windows = false,
    subpath_cgi = true, 
    subpath_cgi_no_time_windows_naive = false,
)


begin
    Bmax_range = [750.0, 800.0, 850.0, 900.0, 950.0]
    Tmul_range = [2.0, 2.5, 3.0]
    for (iB, B) in enumerate(Bmax_range)
        for (iT_mul, T_mul) in enumerate(Tmul_range)
            T = floor(B * T_mul / 50.0) * 50.0
            for seed in 1:10
                run_index = seed
                size = "xs2_$(iB)_$(iT_mul)"
                _, data = generate_instance_pair(
                    n_depots = 2, 
                    n_customers = 9,
                    n_charging = 2,
                    charging_repeats = 1,
                    n_vehicles = 3,
                    shrinkage_depots = 1.4,
                    shrinkage_charging = 0.6,
                    T = T,
                    seed = seed,
                    B = B,
                    μ = 5.0,
                    travel_cost_coeff = 10,
                    charge_cost_coeff = 7,
                    batch = 3,
                    permissiveness = 0.4,
                )
                data = preprocess_arcs(data, true, false)
                G = construct_graph(data)
                T_range = 0:50.0:data["T"]
                B_range = 0:50.0:data["B"]
                if !(size in keys(all_data))
                    all_data[size] = Dict()
                end
                if !(run_index in keys(all_data[size]))
                    all_data[size][run_index] = Dict(
                        "data" => data,
                        "G" => G,
                        "T_range" => T_range,
                        "B_range" => B_range,
                    )
                end
            end
        end
    end
    compare_formulations!(
        all_data, 
        vec(["xs2_$(i)_$(j)" for i in 1:5, j in 1:3]),
        collect(1:10),
        with_charging_cost = true,
        cg_charge_to_full_only = false,
        time_windows = false,
        subpath_cgi = true, 
        subpath_cgi_no_time_windows_naive = false,
    )
end

begin
    Bmax_range = [1000.0, 1050.0, 1100.0, 1150.0, 1200.0]
    Tmul_range = [2.0, 2.5, 3.0]
    for (iB, B) in enumerate(Bmax_range)
        for (iT_mul, T_mul) in enumerate(Tmul_range)
            T = floor(B * T_mul / 50.0) * 50.0
            for seed in 1:10
                run_index = seed
                size = "s2_$(iB)_$(iT_mul)"
                _, data = generate_instance_pair(
                    n_depots = 2, 
                    n_customers = 12,
                    n_charging = 2,
                    charging_repeats = 1,
                    n_vehicles = 3,
                    shrinkage_depots = 1.4,
                    shrinkage_charging = 0.6,
                    T = T,
                    seed = seed,
                    B = B,
                    μ = 5.0,
                    travel_cost_coeff = 10,
                    charge_cost_coeff = 7,
                    batch = 4,
                    permissiveness = 0.4,
                )
                data = preprocess_arcs(data, true, false)
                G = construct_graph(data)
                T_range = 0:50.0:data["T"]
                B_range = 0:50.0:data["B"]
                if !(size in keys(all_data))
                    all_data[size] = Dict()
                end
                if !(run_index in keys(all_data[size]))
                    all_data[size][run_index] = Dict(
                        "data" => data,
                        "G" => G,
                        "T_range" => T_range,
                        "B_range" => B_range,
                    )
                end
            end
        end
    end
    compare_formulations!(
        all_data, 
        vec(["s2_$(i)_$(j)" for i in 1:5, j in 1:3]),
        collect(1:10),
        with_charging_cost = true,
        cg_charge_to_full_only = false,
        time_windows = false,
        subpath_cgi = true, 
        subpath_cgi_no_time_windows_naive = false,
    )
end

compare_formulations!(
    all_data,
    [
        "xs2",
        "s2",
        "m2",
        "l2",
    ], 
    collect(1:10),
    with_charging_cost = true,
    cg_charge_to_full_only = false,
    cg_with_heuristic = true,
    subpath_cgi = true, 
    time_windows = false,
    subpath_cgi_no_time_windows_naive = false,
)

all_metrics_df = compare_formulations!(
    all_data, 
    ["xs", "s", "m"], collect(1:10),
    cg_charge_to_full_only = false,
    subpath_cgi = true, 
    time_windows = false,
    subpath_cgi_no_time_windows_naive = true,
)

compare_formulations!(
    all_data, 
    ["xs", "s", "m"], collect(1:10),
    arc = true,
    cg_charge_to_full_only = false,
    subpath_cgi = true, 
    subpath_cgi_time_windows = false,
)


compare_formulations!(
    all_data, ["xs", "s", "m"], collect(1:10),
    cg_charge_to_full_only = true,
    subpath_cg = true, 
    subpath_cgi = true, 
    path_cg = true,
)
all_metrics_df = compare_formulations!(
    all_data, ["xs", "s", "m"], [2],
    arc = true,
    # subpath_cg = true, 
    # subpath_cgi = true, 
    # path_cg = true,
)

all_metrics_df = compare_formulations!(
    all_data, 
    # ["xs"], [1],
    ["xs", "s"], collect(1:10),
    subpath_cgi = true, 
    subpath_cgi_time_windows = false,
    subpath_cgi_no_time_windows_naive = true,
    cg_charge_to_full_only = true,
)
all_metrics_df = compare_formulations!(
    all_data, 
    # ["xs"], [1],
    ["xs", "s"], collect(1:10),
    subpath_cgi = true, 
    subpath_cgi_time_windows = false,
    subpath_cgi_no_time_windows_naive = false,
    cg_charge_to_full_only = true,
)

## Testing
@test all_data["xs"][1]["cg_results"]["objective"] ≈ all_data["xs"][1]["cgi_results"]["objective"] ≈ all_data["xs"][1]["path_cg_results"]["objective"] ≈ 2584.0
@test all_data["xs"][2]["cg_results"]["objective"] ≈ all_data["xs"][2]["cgi_results"]["objective"] ≈ all_data["xs"][2]["path_cg_results"]["objective"] ≈ 2831.0
@test all_data["xs"][3]["cg_results"]["objective"] ≈ all_data["xs"][3]["cgi_results"]["objective"] ≈ all_data["xs"][3]["path_cg_results"]["objective"] ≈ 2483.5
@test all_data["xs"][4]["cg_results"]["objective"] ≈ all_data["xs"][4]["cgi_results"]["objective"] ≈ all_data["xs"][4]["path_cg_results"]["objective"] ≈ 3357.0
@test all_data["xs"][5]["cg_results"]["objective"] ≈ all_data["xs"][5]["cgi_results"]["objective"] ≈ all_data["xs"][5]["path_cg_results"]["objective"] ≈ 2766.0
@test all_data["xs"][6]["cg_results"]["objective"] ≈ all_data["xs"][6]["cgi_results"]["objective"] ≈ all_data["xs"][6]["path_cg_results"]["objective"] ≈ 2710.5
@test all_data["xs"][7]["cg_results"]["objective"] ≈ all_data["xs"][7]["cgi_results"]["objective"] ≈ all_data["xs"][7]["path_cg_results"]["objective"] ≈ 3600 + 5 / 12
@test all_data["xs"][8]["cg_results"]["objective"] ≈ all_data["xs"][8]["cgi_results"]["objective"] ≈ all_data["xs"][8]["path_cg_results"]["objective"] ≈ 2124.0
@test all_data["xs"][9]["cg_results"]["objective"] ≈ all_data["xs"][9]["cgi_results"]["objective"] ≈ all_data["xs"][9]["path_cg_results"]["objective"] ≈ 3759.0
@test all_data["xs"][10]["cg_results"]["objective"] ≈ all_data["xs"][10]["cgi_results"]["objective"] ≈ all_data["xs"][10]["path_cg_results"]["objective"] ≈ 4019.0

@test all_data["s"][1]["cg_results"]["objective"] ≈ all_data["s"][1]["cgi_results"]["objective"] ≈ all_data["s"][1]["path_cg_results"]["objective"] ≈ 3878.0
@test all_data["s"][2]["cg_results"]["objective"] ≈ all_data["s"][2]["cgi_results"]["objective"] ≈ all_data["s"][2]["path_cg_results"]["objective"] ≈ 3605.5
@test all_data["s"][3]["cg_results"]["objective"] ≈ all_data["s"][3]["cgi_results"]["objective"] ≈ all_data["s"][3]["path_cg_results"]["objective"] ≈ 2941.0
@test all_data["s"][4]["cg_results"]["objective"] ≈ all_data["s"][4]["cgi_results"]["objective"] ≈ all_data["s"][4]["path_cg_results"]["objective"] ≈ 3835.0
@test all_data["s"][5]["cg_results"]["objective"] ≈ all_data["s"][5]["cgi_results"]["objective"] ≈ all_data["s"][5]["path_cg_results"]["objective"] ≈ 3749.0
@test all_data["s"][6]["cg_results"]["objective"] ≈ all_data["s"][6]["cgi_results"]["objective"] ≈ all_data["s"][6]["path_cg_results"]["objective"] ≈ 3500 + 3 / 17
@test all_data["s"][7]["cg_results"]["objective"] ≈ all_data["s"][7]["cgi_results"]["objective"] ≈ all_data["s"][7]["path_cg_results"]["objective"] ≈ 4342 + 2 / 5
@test all_data["s"][8]["cg_results"]["objective"] ≈ all_data["s"][8]["cgi_results"]["objective"] ≈ all_data["s"][8]["path_cg_results"]["objective"] ≈ 3025 + 2 / 5
@test all_data["s"][9]["cg_results"]["objective"] ≈ all_data["s"][9]["cgi_results"]["objective"] ≈ all_data["s"][9]["path_cg_results"]["objective"] ≈ 4199.5
@test all_data["s"][10]["cg_results"]["objective"] ≈ all_data["s"][10]["cgi_results"]["objective"] ≈ all_data["s"][10]["path_cg_results"]["objective"] ≈ 4580.0

## Comparison


CSV.write("$(@__DIR__)/../data/all_metrics.csv", all_metrics_df)
all_metrics_df = CSV.read("$(@__DIR__)/../data/all_metrics.csv", DataFrame)

all_objective_metrics = all_metrics_df[!, [
    :size, 
    :run_index,
    :arc_ip_objective, 
    :arc_lp_objective,
    :cg_objective,
    :cg_lpip_objective,
    :cgi_objective,
    :cgi_lpip_objective,
    :path_cg_objective,
    :path_cg_lpip_objective,
    :enum_ip_objective,
    :enum_lp_objective,
]]

all_objective_metrics

all_time_metrics = all_metrics_df[!, [
    :size, 
    :run_index,
    # :arc_ip_time_taken, 
    # :arc_lp_time_taken,
    :cg_time_taken,
    :cg_mp_total_time_taken,
    :cg_mp_total_constraint_time_taken,
    :cg_mp_total_solution_time_taken,
    :cg_sp_total_time_taken,
    # :cg_mp_mean_time_taken,
    # :cg_mp_mean_constraint_time_taken,
    # :cg_mp_mean_solution_time_taken,
    # :cg_sp_mean_time_taken,
    :cg_lpip_time_taken,
    :cgi_time_taken,
    :cgi_mp_total_time_taken,
    :cgi_mp_total_constraint_time_taken,
    :cgi_mp_total_solution_time_taken,
    :cgi_sp_total_time_taken,
    # :cgi_mp_mean_time_taken,
    # :cgi_mp_mean_constraint_time_taken,
    # :cgi_mp_mean_solution_time_taken,
    # :cgi_sp_mean_time_taken,
    :cgi_lpip_time_taken,
    :path_cg_time_taken,
    :path_cg_mp_total_time_taken,
    :path_cg_mp_total_constraint_time_taken,
    :path_cg_mp_total_solution_time_taken,
    :path_cg_sp_total_time_taken,
    # :path_cg_mp_mean_time_taken,
    # :path_cg_mp_mean_constraint_time_taken,
    # :path_cg_mp_mean_solution_time_taken,
    # :path_cg_sp_mean_time_taken,
    :path_cg_lpip_time_taken,
    # :enumerate_all_subpaths_time_taken,
    # :enum_ip_time_taken,
    # :enum_lp_time_taken,
]]


all_time_metrics
all_time_metrics |>
    x -> filter(r -> (r.size == "s"), x) |>
    x -> select(x, [
        :cgi_time_taken,     
        :cgi_mp_total_time_taken,
        :cgi_mp_total_constraint_time_taken,
        :cgi_mp_total_solution_time_taken,
        :cgi_sp_total_time_taken,
    ])

all_time_metrics |>
    x -> filter(r -> (r.size == "xs"), x) |>
    x -> select(x, [:cg_time_taken, :cgi_time_taken, :path_cg_time_taken]) |>
    x -> describe(x, :detailed)

all_time_metrics |>
    x -> filter(r -> (r.size == "xs"), x) |>
    x -> select(x, Not([:size, :run_index])) |>
    x -> describe(x, :detailed)

all_time_metrics |>
    x -> filter(r -> (r.size == "s"), x) |>
    x -> select(x, [:cg_time_taken, :cgi_time_taken, :path_cg_time_taken]) |>
    x -> describe(x, :detailed)

all_time_metrics |>
    x -> filter(r -> (r.size == "s"), x) |>
    x -> select(x, Not([:size, :run_index])) |>
    x -> describe(x, :detailed)

# Experiment 1: comparing different formulations
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230216_105009/all_metrics.csv", DataFrame)
time_taken_df = all_metrics_df |>
    x -> select(x, [:size, :run_index, :cg_time_taken, :cgi_time_taken, :path_cg_time_taken]) |>
    x -> stack(x, [:cg_time_taken, :cgi_time_taken, :path_cg_time_taken])

@df filter(r -> (r.size == "xs"), time_taken_df) boxplot(
    string.(:variable),
    :value,
    fill_alpha = 0.5,
    title = "Running time of different algorithms (smallest, n = 9)",
    ylabel = "Time (s)",
    legend = :topleft,
)
@df filter(r -> (r.size == "xs"), time_taken_df) dotplot!(
    string.(:variable),
    :value,
    marker=(:black, stroke(0)),
)

@df filter(r -> (r.size == "s"), time_taken_df) boxplot(
    string.(:variable),
    :value,
    fill_alpha = 0.5,
    title = "Running time of different algorithms (small, n = 12)",
    ylabel = "Time (s)",
    legend = :topleft,
)
@df filter(r -> (r.size == "s"), time_taken_df) dotplot!(
    string.(:variable),
    :value,
    marker=(:black, stroke(0))
)

# Experiment 2: varying battery
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230216_121329/all_metrics.csv", DataFrame)

@df all_metrics_df boxplot(
    string.(:size),
    :cgi_time_taken,
    fill_alpha = 0.5,
    title = "Running time, varying battery (smallest, n = 9)",
    ylabel = "Time (s)",
    xlabel = "Battery capacity"
)
@df all_metrics_df dotplot!(
    string.(:size),
    :cgi_time_taken,
    marker=(:black, stroke(0))
)

# Experiment 3: varying length of time windows
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230216_130511/all_metrics.csv", DataFrame)

@df all_metrics_df boxplot(
    string.(:size),
    :cgi_time_taken,
    fill_alpha = 0.5,
    title = "Running time, varying time windows (smallest, n = 9)",
    ylabel = "Time (s)",
    xlabel = "Time window size",
    legend = :topleft,
)
@df all_metrics_df dotplot!(
    string.(:size),
    :cgi_time_taken,
    marker=(:black, stroke(0))
)

# Experiment 4: adding medium-sized entries
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230221_204918/all_metrics.csv", DataFrame)
all_metrics_df |>
    x -> filter(r -> (r.run_index == 2), x) |>
    x -> select(x, [:cg_time_taken, :cgi_time_taken, :path_cg_time_taken])
all_metrics_df |>
    x -> filter(r -> (r.size == "m"), x) |>
    x -> select(x, [:cg_time_taken, :cgi_time_taken, :path_cg_time_taken]) |>
    x -> describe(x, :detailed)

time_taken_df = all_metrics_df |>
    x -> select(x, [:size, :run_index, :cg_time_taken, :cgi_time_taken, :path_cg_time_taken]) |>
    x -> stack(x, [:cg_time_taken, :cgi_time_taken, :path_cg_time_taken])

@df filter(r -> (r.size == "m"), time_taken_df) boxplot(
    string.(:variable),
    :value,
    fill_alpha = 0.5,
    title = "Running time of different algorithms (medium, n = 15)",
    ylabel = "Time (s)",
    legend = :topleft,
    label = false,
)
@df filter(r -> (r.size == "m"), time_taken_df) dotplot!(
    string.(:variable),
    :value,
    marker = (:black, stroke(0)),
    label = false,
)

@df time_taken_df groupedboxplot(
    string.(:variable),
    :value,
    group = :size,
    title = "Comparing different algorithms for n = 9, 12, 15",
    ylabel = "Time (s)",
)
@df time_taken_df groupeddotplot!(
    string.(:variable),
    :value,
    group = :size,
    marker = (:black, stroke(0)),
    label = false,
)

# Experiment 6
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230227_200940/all_metrics.csv", DataFrame)
all_metrics_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :arc_ip_time_taken => mean,
        :cg_time_taken => mean,
        :cgi_time_taken => mean,
        :path_cg_time_taken => mean,
        :cg_mp_total_time_taken => mean,
        :cgi_mp_total_time_taken => mean,
        :path_cg_mp_total_time_taken => mean,
        :cg_sp_mean_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
        :path_cg_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

all_metrics_full_df = CSV.read("$(@__DIR__)/../logs/20230228_011832/all_metrics.csv", DataFrame)
all_metrics_full_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :arc_ip_time_taken => mean,
        :cg_time_taken => mean,
        :cgi_time_taken => mean,
        :path_cg_time_taken => mean,
        :cg_mp_total_time_taken => mean,
        :cgi_mp_total_time_taken => mean,
        :path_cg_mp_total_time_taken => mean,
        :cg_sp_mean_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
        :path_cg_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

# Experiment 7: two different methods in the no-time-windows situation

# naive
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230228_182555/all_metrics.csv", DataFrame)
all_metrics_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :arc_ip_time_taken => mean,
        :cgi_time_taken => mean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

# charge_bounded
all_metrics_df_1 = CSV.read("$(@__DIR__)/../logs/20230301_114713/all_metrics.csv", DataFrame)
all_metrics_df_1 |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :arc_ip_time_taken => mean,
        :cgi_time_taken => mean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

# charge_bounded, charge_to_full_only
all_metrics_df_2 = CSV.read("$(@__DIR__)/../logs/20230303_143104/all_metrics.csv", DataFrame)
all_metrics_df_2 |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :arc_ip_time_taken => mean,
        :cgi_time_taken => mean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")
[
    all_metrics_df_2[!, :cgi_objective]  all_metrics_df_1[!, :cgi_objective]
]

# Experiment 8: two different domination criteria
# (20230303_192829) many subpaths, less iterations, correct convergence
# (20230304_105525) far fewer subpaths, much faster subproblems, more iterations (tailing-off), investigate convergence
# investigate: T:B ratio
all_metrics_df_1 = CSV.read("$(@__DIR__)/../logs/20230303_192829/all_metrics.csv", DataFrame)
all_metrics_df_1 |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

all_metrics_df_2 = CSV.read("$(@__DIR__)/../logs/20230304_105525/all_metrics.csv", DataFrame)
all_metrics_df_2 |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

# Experimetn 11
all_metrics_df_1 = CSV.read("$(@__DIR__)/../logs/20230307_135734/all_metrics.csv", DataFrame)
all_metrics_df_1 |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")
all_metrics_df_2 = CSV.read("$(@__DIR__)/../logs/20230307_114101/all_metrics.csv", DataFrame)
all_metrics_df_2 |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")
[all_metrics_df_1[!, :cgi_objective] all_metrics_df_2[!, :cgi_objective]]

# Experiment 12: include charge time and customer delay in objective
all_metrics_xs_df = CSV.read("$(@__DIR__)/../logs/20230313_110901/all_metrics.csv", DataFrame)
all_metrics_xs_df |>
    x -> filter(r -> r.size == "xs2_3_3", x) |>
    x -> select(x, :cgi_objective)
gdf_xs = all_metrics_xs_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    )
Bmax_range_xs = [750.0, 800.0, 850.0, 900.0, 950.0]
Tmul_range_xs = [2.0, 2.5, 3.0]
gdf_xs[!, :B] .= repeat(Bmax_range_xs, outer = length(Tmul_range_xs))
gdf_xs[!, :T_mul] .= repeat(Tmul_range_xs, inner = length(Bmax_range_xs))
gdf_xs[!, :T] .= round.(gdf_xs[!, :B] .* gdf_xs[!, :T_mul] ./ 50.0) .* 50.0
select!(gdf_xs, [:size, :B, :T_mul, :T], Not([:size, :B, :T_mul, :T]))


begin
    fig = Figure(resolution = (800, 500), fontsize = 18)
    grid = fig[1,1] = GridLayout()
    ax = Axis(
        grid[1,1],
        xticks = Bmax_range_xs, 
        xlabel = "B",
        yticks = Tmul_range_xs,
        ylabel = "Ratio T : B",
        title = "Running time of subproblem (mean over iterations)\n(9 customers, 2 depots, 2 CS)"
    )
    mat = Matrix(
        unstack(gdf_xs, :B, :T_mul, :cgi_sp_mean_time_taken_mean)[:,2:end]
    )
    (cmin, cmax) = extrema(mat)
    CairoMakie.heatmap!(
        ax, Bmax_range_xs, Tmul_range_xs,
        mat,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_xs), j in 1:length(Tmul_range_xs)
        textcolor = mat[i, j] < 0.3 ? :white : :black
        text!(
            ax, "$(round(mat[i,j], digits = 4))",
            position = (Bmax_range_xs[i], Tmul_range_xs[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    Colorbar(grid[1,2], colorrange = (cmin, cmax))
    display(fig)
end

begin
    fig = Figure(resolution = (800, 500), fontsize = 18)
    grid = fig[1,1] = GridLayout()
    ax = Axis(
        grid[1,1],
        xticks = Bmax_range_xs, 
        xlabel = "B",
        yticks = Tmul_range_xs,
        ylabel = "Ratio T : B",
        title = "Total running time\n(9 customers, 2 depots, 2 CS)"
    )
    mat = Matrix(
        unstack(gdf_xs, :B, :T_mul, :cgi_time_taken_geomean)[:,2:end]
    )
    (cmin, cmax) = extrema(mat)
    CairoMakie.heatmap!(
        ax, Bmax_range_xs, Tmul_range_xs,
        mat,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_xs), j in 1:length(Tmul_range_xs)
        textcolor = mat[i, j] < 2.5 ? :white : :black
        text!(
            ax, "$(round(mat[i,j], digits = 4))",
            position = (Bmax_range_xs[i], Tmul_range_xs[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    Colorbar(grid[1,2], colorrange = (cmin, cmax))
    display(fig)
end

all_metrics_xs_r_df = CSV.read("$(@__DIR__)/../logs/20230313_165734/all_metrics.csv", DataFrame)
all_metrics_xs_r_df |>
    x -> filter(r -> r.size == "xs2_5_3", x) |>
    x -> select(x, :cgi_objective)
gdf_xs_r = all_metrics_xs_r_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    )
Bmax_range_xs = [750.0, 800.0, 850.0, 900.0, 950.0]
Tmul_range_xs = [2.0, 2.5, 3.0]
gdf_xs_r[!, :B] .= repeat(Bmax_range_xs, outer = length(Tmul_range_xs))
gdf_xs_r[!, :T_mul] .= repeat(Tmul_range_xs, inner = length(Bmax_range_xs))
gdf_xs_r[!, :T] .= round.(gdf_xs_r[!, :B] .* gdf_xs_r[!, :T_mul] ./ 50.0) .* 50.0
select!(gdf_xs_r, [:size, :B, :T_mul, :T], Not([:size, :B, :T_mul, :T]))

begin
    mat_xs = Matrix(
        unstack(gdf_xs, :B, :T_mul, :cgi_sp_mean_time_taken_mean)[:,2:end]
    )
    mat_xs_r = Matrix(
        unstack(gdf_xs_r, :B, :T_mul, :cgi_sp_mean_time_taken_mean)[:,2:end]
    )
    fig = Figure(resolution = (1200, 500), fontsize = 18)
    grid = fig[1,1] = GridLayout()
    cmin = min(minimum(mat_xs), minimum(mat_xs_r))
    cmax = max(maximum(mat_xs), maximum(mat_xs_r))
    ax_xs = Axis(
        grid[1,1],
        xticks = Bmax_range_xs, 
        xlabel = "B",
        yticks = Tmul_range_xs,
        ylabel = "Ratio T : B",
        title = "Running time of subproblem (mean over iterations)\n(9 customers, 2 depots, 2 CS)"
    )
    CairoMakie.heatmap!(
        ax_xs, Bmax_range_xs, Tmul_range_xs,
        mat_xs,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_xs), j in 1:length(Tmul_range_xs)
        textcolor = mat_xs[i, j] < 0.3 ? :white : :black
        text!(
            ax_xs, "$(round(mat_xs[i,j], digits = 4))",
            position = (Bmax_range_xs[i], Tmul_range_xs[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    ax_xs_r = Axis(
        grid[1,2],
        xticks = Bmax_range_xs, 
        xlabel = "B",
        yticks = Tmul_range_xs,
        ylabel = "Ratio T : B",
        title = "Running time of subproblem (mean over iterations)\n(9 customers, 2 depots, 2 CS)\n(more aggressive domination)"
    )
    CairoMakie.heatmap!(
        ax_xs_r, Bmax_range_xs, Tmul_range_xs,
        mat_xs_r,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_xs), j in 1:length(Tmul_range_xs)
        textcolor = mat_xs_r[i, j] < 0.3 ? :white : :black
        text!(
            ax_xs_r, "$(round(mat_xs_r[i,j], digits = 4))",
            position = (Bmax_range_xs[i], Tmul_range_xs[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    Colorbar(grid[1,3], colorrange = (cmin, cmax))
    display(fig)
end

begin
    mat_xs = Matrix(
        unstack(gdf_xs, :B, :T_mul, :cgi_time_taken_geomean)[:,2:end]
    )
    mat_xs_r = Matrix(
        unstack(gdf_xs_r, :B, :T_mul, :cgi_time_taken_geomean)[:,2:end]
    )
    fig = Figure(resolution = (1200, 500), fontsize = 18)
    grid = fig[1,1] = GridLayout()
    cmin = min(minimum(mat_xs), minimum(mat_xs_r))
    cmax = max(maximum(mat_xs), maximum(mat_xs_r))
    ax_xs = Axis(
        grid[1,1],
        xticks = Bmax_range_xs, 
        xlabel = "B",
        yticks = Tmul_range_xs,
        ylabel = "Ratio T : B",
        title = "Total running time\n(9 customers, 2 depots, 2 CS)"
    )
    CairoMakie.heatmap!(
        ax_xs, Bmax_range_xs, Tmul_range_xs,
        mat_xs,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_xs), j in 1:length(Tmul_range_xs)
        textcolor = mat_xs[i, j] < 2.5 ? :white : :black
        text!(
            ax_xs, "$(round(mat_xs[i,j], digits = 4))",
            position = (Bmax_range_xs[i], Tmul_range_xs[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    ax_xs_r = Axis(
        grid[1,2],
        xticks = Bmax_range_xs, 
        xlabel = "B",
        yticks = Tmul_range_xs,
        ylabel = "Ratio T : B",
        title = "Total running time\n(9 customers, 2 depots, 2 CS)\n(more aggressive domination)"
    )
    CairoMakie.heatmap!(
        ax_xs_r, Bmax_range_xs, Tmul_range_xs,
        mat_xs_r,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_xs), j in 1:length(Tmul_range_xs)
        textcolor = mat_xs_r[i, j] < 2.5 ? :white : :black
        text!(
            ax_xs_r, "$(round(mat_xs_r[i,j], digits = 4))",
            position = (Bmax_range_xs[i], Tmul_range_xs[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    Colorbar(grid[1,3], colorrange = (cmin, cmax))
    display(fig)
end


all_metrics_s_df = CSV.read("$(@__DIR__)/../logs/20230313_120650/all_metrics.csv", DataFrame)
all_metrics_s_df |>
    x -> filter(r -> r.size == "s2_1_1", x) |>
    x -> select(x, :cgi_objective)
gdf_s = all_metrics_s_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    )
Bmax_range_s = [1000.0, 1050.0, 1100.0, 1150.0, 1200.0]
Tmul_range_s = [2.0, 2.5, 3.0]
gdf_s[!, :B] .= repeat(Bmax_range_s, outer = length(Tmul_range_s))
gdf_s[!, :T_mul] .= repeat(Tmul_range_s, inner = length(Bmax_range_s))
gdf_s[!, :T] .= round.(gdf_s[!, :B] .* gdf_s[!, :T_mul] ./ 50.0) .* 50.0
select!(gdf_s, [:size, :B, :T_mul, :T], Not([:size, :B, :T_mul, :T]))

begin
    fig = Figure(resolution = (800, 500), fontsize = 18)
    grid = fig[1,1] = GridLayout()
    ax = Axis(
        grid[1,1],
        xticks = Bmax_range_s, 
        xlabel = "B",
        yticks = Tmul_range_s,
        ylabel = "Ratio T : B",
        title = "Running time of subproblem (mean over iterations)\n(12 customers, 2 depots, 2 CS)"
    )
    mat = Matrix(
        unstack(gdf_s, :B, :T_mul, :cgi_sp_mean_time_taken_mean)[:,2:end]
    )
    (cmin, cmax) = extrema(mat)
    CairoMakie.heatmap!(
        ax, Bmax_range_s, Tmul_range_s,
        mat,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_s), j in 1:length(Tmul_range_s)
        textcolor = mat[i, j] < 2.4 ? :white : :black
        text!(
            ax, "$(round(mat[i,j], digits = 4))",
            position = (Bmax_range_s[i], Tmul_range_s[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    Colorbar(grid[1,2], colorrange = (cmin, cmax))
    display(fig)
end

begin
    fig = Figure(resolution = (800, 500), fontsize = 18)
    grid = fig[1,1] = GridLayout()
    ax = Axis(
        grid[1,1],
        xticks = Bmax_range_s, 
        xlabel = "B",
        yticks = Tmul_range_s,
        ylabel = "Ratio T : B",
        title = "Total running time\n(12 customers, 2 depots, 2 CS)"
    )
    mat = Matrix(
        unstack(gdf_s, :B, :T_mul, :cgi_time_taken_geomean)[:,2:end]
    )
    (cmin, cmax) = extrema(mat)
    CairoMakie.heatmap!(
        ax, Bmax_range_s, Tmul_range_s,
        mat,
        colorrange = (cmin, cmax),
    )
    for i in 1:length(Bmax_range_s), j in 1:length(Tmul_range_s)
        textcolor = mat[i, j] < 25 ? :white : :black
        text!(
            ax, "$(round(mat[i,j], digits = 4))",
            position = (Bmax_range_s[i], Tmul_range_s[j]),
            color = textcolor, 
            align = (:center, :center),
            fontsize = 10,
        )
    end
    Colorbar(grid[1,2], colorrange = (cmin, cmax))
    display(fig)
end

# Experiment 13: with rougher dominance
# naive, not utilizing notimewindows structure
all_metrics_1_df = CSV.read("$(@__DIR__)/../logs/20230313_203258/all_metrics.csv", DataFrame)
gdf1 = all_metrics_1_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

# previous method utilizing notimewindows structure
all_metrics_2_df = CSV.read("$(@__DIR__)/../logs/20230314_100637/all_metrics.csv", DataFrame)
gdf2 = all_metrics_2_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")

# method utilizing notimewindows structure, more aggressive dominance
all_metrics_3_df = CSV.read("$(@__DIR__)/../logs/20230313_200031/all_metrics.csv", DataFrame)
gdf3 = all_metrics_3_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken => geomean,
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken => mean,
    ) |>
    x -> permutedims(x, "size")


## Experiment 14: heuristic 
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230313_200031/all_metrics.csv", DataFrame)
gdf = all_metrics_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths => mean,
        :cgi_time_taken .=> [minimum, maximum, geomean],
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken .=> [minimum, maximum, mean],
    ) |>
    x -> permutedims(x, "size")

all_metrics_df_h = CSV.read("$(@__DIR__)/../logs/20230315_171528/all_metrics.csv", DataFrame)
gdf_h = all_metrics_df_h |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths .=> mean,
        :cgi_time_taken .=> [minimum, maximum, geomean],
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken .=> [minimum, maximum, mean],
    ) |>
    x -> permutedims(x, "size")

## Experiment 15: removing discretization!
all_metrics_df = CSV.read("$(@__DIR__)/../logs/20230315_171528/all_metrics.csv", DataFrame)
gdf = all_metrics_df |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths .=> mean,
        :cgi_time_taken .=> [minimum, maximum, geomean],
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken .=> [minimum, maximum, mean],
    ) |>
    x -> permutedims(x, "size")

all_metrics_df_new = CSV.read("$(@__DIR__)/../logs/20230320_212038/all_metrics.csv", DataFrame)
gdf_new = all_metrics_df_new |>
    x -> groupby(x, :size) |>
    x -> combine(x, 
        :cgi_number_of_iterations => mean,
        :cgi_number_of_subpaths .=> mean,
        :cgi_time_taken .=> [minimum, maximum, geomean],
        :cgi_mp_total_time_taken => mean,
        :cgi_sp_mean_time_taken .=> [minimum, maximum, mean],
    ) |>
    x -> permutedims(x, "size")


### Scratch work