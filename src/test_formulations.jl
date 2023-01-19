include("arc_formulation.jl")
include("subpath_formulation.jl")
include("utils.jl")

using DataFrames
using Dates
using CSV, JLD2

JLD2.save_object("/data/all_data.jld", all_data)
all_data = JLD2.load_object("/data/all_data.jld")

all_data = Dict(
    "xs" => Dict(),
    "s" => Dict(),
    "m" => Dict(),
)

_, data = generate_instance_pair(
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
data = preprocess_arcs(data, true, false)
G = construct_graph(data)
T_range = 0:50.0:data["T"]
B_range = 0:50.0:data["B"]
all_data["xs"][1] = Dict(
    "data" => data,
    "G" => G,
    "T_range" => T_range,
    "B_range" => B_range,
)

_, data = generate_instance_pair(
    n_depots = 2, 
    n_customers = 9,
    n_charging = 2,
    charging_repeats = 1,
    n_vehicles = 3,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.6,
    T = 1350.0,
    seed = 2,
    B = 800.0,
    μ = 5.0,
    batch = 3,
)
data = preprocess_arcs(data, true, false)
G = construct_graph(data)
T_range = 0:50.0:data["T"]
B_range = 0:50.0:data["B"]
all_data["xs"][2] = Dict(
    "data" => data,
    "G" => G,
    "T_range" => T_range,
    "B_range" => B_range,
)



_, data = generate_instance_pair(
    n_depots = 2, 
    n_customers = 12,
    n_charging = 2,
    charging_repeats = 1,
    n_vehicles = 3,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.6,
    T = 1800.0,
    seed = 1,
    B = 1400.0,
    μ = 5.0,
    batch = 4,
)
data = preprocess_arcs(data, true, false)
G = construct_graph(data)
T_range = 0:50.0:data["T"]
B_range = 0:50.0:data["B"]
all_data["s"][1] = Dict(
    "data" => data,
    "G" => G,
    "T_range" => T_range,
    "B_range" => B_range,
)

_, data = generate_instance_pair(
    n_depots = 2, 
    n_customers = 15,
    n_charging = 2,
    charging_repeats = 1,
    n_vehicles = 3,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.6,
    T = 2250.0,
    seed = 1,
    B = 1550.0,
    μ = 5.0,
    batch = 5,
)
data = preprocess_arcs(data, true, false)
G = construct_graph(data)
T_range = 0:50.0:data["T"]
B_range = 0:50.0:data["B"]
all_data["m"][1] = Dict(
    "data" => data,
    "G" => G,
    "T_range" => T_range,
    "B_range" => B_range,
)

### Arc formulation

(
    all_data["xs"][1]["arc_ip_results"], 
    all_data["xs"][1]["arc_ip_params"],
) = arc_formulation(
    all_data["xs"][1]["data"],
    true,
    ;
    time_limit = 300.0,
)
arc_results_printout(
    all_data["xs"][1]["arc_ip_results"], 
    all_data["xs"][1]["arc_ip_params"],
    all_data["xs"][1]["data"],
    true,
)

(
    all_data["xs"][1]["arc_lp_results"], 
    all_data["xs"][1]["arc_lp_params"],
) = arc_formulation(
    all_data["xs"][1]["data"],
    true,
    ;
    integral = false,
    time_limit = 300.0,
)

(
    all_data["s"][1]["arc_ip_results"], 
    all_data["s"][1]["arc_ip_params"],
) = arc_formulation(
    all_data["s"][1]["data"],
    true,
    ;
    time_limit = 300.0,
)
arc_results_printout(
    all_data["s"][1]["arc_ip_results"], 
    all_data["s"][1]["arc_ip_params"],
    all_data["s"][1]["data"],
    true,
)

(
    all_data["s"][1]["arc_lp_results"], 
    all_data["s"][1]["arc_lp_params"],
) = arc_formulation(
    all_data["s"][1]["data"],
    true,
    ;
    integral = false,
    time_limit = 300.0,
)

(
    all_data["m"][1]["arc_ip_results"], 
    all_data["m"][1]["arc_ip_params"],
) = arc_formulation(
    all_data["m"][1]["data"],
    true,
    ;
    time_limit = 300.0,
)
arc_results_printout(
    all_data["m"][1]["arc_ip_results"], 
    all_data["m"][1]["arc_ip_params"],
    all_data["m"][1]["data"],
    true,
)

(
    all_data["m"][1]["arc_lp_results"], 
    all_data["m"][1]["arc_lp_params"],
) = arc_formulation(
    all_data["m"][1]["data"],
    true,
    ;
    integral = false,
    time_limit = 300.0,
)

# Subpath formulations

## Subpath generation

(
    all_data["xs"][2]["cg_results"],
    all_data["xs"][2]["cg_params"],
    all_data["xs"][2]["cg_printlist"],
    all_data["xs"][2]["cg_subpaths"],
) = subpath_formulation_column_generation(
    all_data["xs"][2]["G"],
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["T_range"],
    all_data["xs"][2]["B_range"],
    ;
    charging_in_subpath = true,
    verbose = true,
);
print.(all_data["xs"][2]["cg_printlist"]);
filename = "cg_lp_xs_2_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))"
open("$(@__DIR__)/logs/$filename.txt", "w") do io
    for message in all_data["xs"][2]["cg_printlist"]
        write(io, message)
    end
end
groupedbar(
    hcat(
        all_data["xs"][2]["cg_params"]["lp_relaxation_time_taken"],
        all_data["xs"][2]["cg_params"]["sp_total_time_taken"],
    ),
    group = repeat(
        ["LP relaxation solve time", "Subproblem solve time"], 
        inner = length(all_data["xs"][2]["cg_params"]["sp_total_time_taken"]),
    ),
    bar_position = :stack,
    framestyle = :box,
    xlabel = "Iteration",
    xticks = 2:length(all_data["xs"][2]["cg_params"]["sp_total_time_taken"]),
    ylabel = "Time (s)",
    title = """
    Time of LP relaxation and subproblem, with artificial starting subpaths
    ($(all_data["xs"][2]["data"]["n_customers"]) customers, $(all_data["xs"][2]["data"]["n_vehicles"]) vehicles, $(all_data["xs"][2]["data"]["n_depots"]) depots, $(all_data["xs"][2]["data"]["n_charging"]) charging stations)
    """,
    size = (800, 600),
)
savefig("$(@__DIR__)/plots/$filename.png")


all_data["xs"][2]["cg_subpath_costs"] = compute_subpath_costs(
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["cg_subpaths"],
)
all_data["xs"][2]["cg_subpath_service"] = compute_subpath_service(
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["cg_subpaths"],
)
(
    all_data["xs"][2]["cg_lpip_results"],
    all_data["xs"][2]["cg_lpip_params"],
) = subpath_formulation(
    all_data["xs"][2]["data"],
    all_data["xs"][2]["cg_subpaths"],
    all_data["xs"][2]["cg_subpath_costs"],
    all_data["xs"][2]["cg_subpath_service"],
    all_data["xs"][2]["T_range"],
    all_data["xs"][2]["B_range"],
    ;
    integral = true,
    charging_in_subpath = true,
)



(
    all_data["s"][1]["cg_results"],
    all_data["s"][1]["cg_params"],
    all_data["s"][1]["cg_printlist"],
    all_data["s"][1]["cg_subpaths"],
) = subpath_formulation_column_generation(
    all_data["s"][1]["G"],
    all_data["s"][1]["data"], 
    all_data["s"][1]["T_range"],
    all_data["s"][1]["B_range"],
    ;
    charging_in_subpath = true,
    verbose = true,
);
print.(all_data["s"][1]["cg_printlist"]);
filename = "cg_lp_s_1_$(Dates.format(Dates.now(), "yyyymmdd_HHMMSS"))"
open("$(@__DIR__)/logs/$filename.txt", "w") do io
    for message in all_data["s"][1]["cg_printlist"]
        write(io, message)
    end
end
groupedbar(
    hcat(
        all_data["s"][1]["cg_params"]["lp_relaxation_time_taken"],
        all_data["s"][1]["cg_params"]["sp_total_time_taken"],
    ),
    group = repeat(
        ["LP relaxation solve time", "Subproblem solve time"], 
        inner = length(all_data["s"][1]["cg_params"]["sp_total_time_taken"]),
    ),
    bar_position = :stack,
    framestyle = :box,
    xlabel = "Iteration",
    xticks = 1:length(all_data["s"][1]["cg_params"]["sp_total_time_taken"]),
    ylabel = "Time (s)",
    title = """
    Time of LP relaxation and subproblem, with artificial starting subpaths
    ($(all_data["s"][1]["data"]["n_customers"]) customers, $(all_data["s"][1]["data"]["n_vehicles"]) vehicles, $(all_data["s"][1]["data"]["n_depots"]) depots, $(all_data["s"][1]["data"]["n_charging"]) charging stations)
    """,
    size = (800, 600),
)
savefig("$(@__DIR__)/plots/$filename.png")


all_data["s"][1]["cg_subpath_costs"] = compute_subpath_costs(
    all_data["s"][1]["data"], 
    all_data["s"][1]["cg_subpaths"],
)
all_data["s"][1]["cg_subpath_service"] = compute_subpath_service(
    all_data["s"][1]["data"], 
    all_data["s"][1]["cg_subpaths"],
)
(
    all_data["s"][1]["cg_lpip_results"],
    all_data["s"][1]["cg_lpip_params"],
) = subpath_formulation(
    all_data["s"][1]["data"],
    all_data["s"][1]["cg_subpaths"],
    all_data["s"][1]["cg_subpath_costs"],
    all_data["s"][1]["cg_subpath_service"],
    all_data["s"][1]["T_range"],
    all_data["s"][1]["B_range"],
    ;
    integral = true,
    charging_in_subpath = true,
)



## Subpath enumeration

### Subpath formulation (charging in subpath)

(
    all_data["xs"][2]["all_subpaths"],
    all_data["xs"][2]["enumerate_all_subpaths_time_taken"],    
) = enumerate_all_subpaths(
    all_data["xs"][2]["G"], 
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["T_range"], 
    all_data["xs"][2]["B_range"]; 
    charging_in_subpath = true,
)
all_data["xs"][2]["enum_number_of_subpaths"] = sum(
    length(v) for v in values(all_data["xs"][2]["all_subpaths"])
)
all_data["xs"][2]["all_subpath_costs"] = compute_subpath_costs(
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["all_subpaths"],
)
all_data["xs"][2]["all_subpath_service"] = compute_subpath_service(
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["all_subpaths"],
)
(
    all_data["xs"][2]["enum_ip_results"],
    all_data["xs"][2]["enum_ip_params"],
) = subpath_formulation(
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["all_subpaths"],
    all_data["xs"][2]["all_subpath_costs"], 
    all_data["xs"][2]["all_subpath_service"],
    all_data["xs"][2]["T_range"],
    all_data["xs"][2]["B_range"],
    ;
    charging_in_subpath = true,
    integral = true,
)
(
    all_data["xs"][2]["enum_lp_results"],
    all_data["xs"][2]["enum_lp_params"],
) = subpath_formulation(
    all_data["xs"][2]["data"], 
    all_data["xs"][2]["all_subpaths"],
    all_data["xs"][2]["all_subpath_costs"], 
    all_data["xs"][2]["all_subpath_service"],
    all_data["xs"][2]["T_range"],
    all_data["xs"][2]["B_range"],
    ;
    charging_in_subpath = true,
    integral = false,
)
println("Number of starting states:\t$(length(all_data["xs"][2]["all_subpaths"],))")
println("Number of subpaths:\t\t$(sum(length(x) for x in values(all_data["xs"][2]["all_subpaths"])))")
@printf("Time taken:\t\t\t%7.1f s\n", all_data["xs"][2]["enumerate_all_subpaths_time_taken"])
@printf("Objective: \t\t\t%6.1f\n", all_data["xs"][2]["enum_ip_results"]["objective"])
@printf("LP Objective:\t\t\t%6.1f\n", all_data["xs"][2]["enum_lp_results"]["objective"])

(
    all_data["s"][1]["all_subpaths"],
    all_data["s"][1]["enumerate_all_subpaths_time_taken"],    
) = enumerate_all_subpaths(
    all_data["s"][1]["G"], 
    all_data["s"][1]["data"], 
    all_data["s"][1]["T_range"], 
    all_data["s"][1]["B_range"]; 
    charging_in_subpath = true,
)
all_data["s"][1]["enum_number_of_subpaths"] = sum(
    length(v) for v in values(all_data["s"][1]["all_subpaths"])
)
all_data["s"][1]["all_subpath_costs"] = compute_subpath_costs(
    all_data["s"][1]["data"], 
    all_data["s"][1]["all_subpaths"],
)
all_data["s"][1]["all_subpath_service"] = compute_subpath_service(
    all_data["s"][1]["data"], 
    all_data["s"][1]["all_subpaths"],
)
(
    all_data["s"][1]["enum_ip_results"],
    all_data["s"][1]["enum_ip_params"],
) = subpath_formulation(
    all_data["s"][1]["data"], 
    all_data["s"][1]["all_subpaths"],
    all_data["s"][1]["all_subpath_costs"], 
    all_data["s"][1]["all_subpath_service"],
    all_data["s"][1]["T_range"],
    all_data["s"][1]["B_range"],
    ;
    charging_in_subpath = true,
    integral = true,
)
(
    all_data["s"][1]["enum_lp_results"],
    all_data["s"][1]["enum_lp_params"],
) = subpath_formulation(
    all_data["s"][1]["data"], 
    all_data["s"][1]["all_subpaths"],
    all_data["s"][1]["all_subpath_costs"], 
    all_data["s"][1]["all_subpath_service"],
    all_data["s"][1]["T_range"],
    all_data["s"][1]["B_range"],
    ;
    charging_in_subpath = true,
    integral = false,
)
println("Number of starting states:\t$(length(all_data["s"][1]["all_subpaths"],))")
println("Number of subpaths:\t\t$(sum(length(x) for x in all_data["s"][1]["all_subpaths"],))")
@printf("Time taken:\t\t\t%7.1f s\n", all_data["s"][1]["enumerate_all_subpaths_time_taken"])
@printf("Objective: \t\t\t%6.1f\n", all_data["s"][1]["enum_ip_results"]["objective"])
@printf("LP Objective:\t\t\t%6.1f\n", all_data["s"][1]["enum_lp_results"]["objective"])

## Comparison

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
    :cg_time_taken,
    :cg_lp_relaxation_time_taken,
    :cg_sp_total_time_taken,
    :cg_lpip_objective,
    :cg_lpip_time_taken,
    :cg_lpip_solution_time_taken,
    :cg_lpip_constraint_time_taken,
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
for (size, run_index) in [
    ("xs", 1), 
    ("xs", 2),
    ("s", 1),
]
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
        d[:cg_number_of_subpaths] = all_data[size][run_index]["cg_params"]["number_of_subpaths"][end]
        d[:cg_time_taken] = all_data[size][run_index]["cg_params"]["time_taken"]
        d[:cg_lp_relaxation_time_taken] = sum(all_data[size][run_index]["cg_params"]["lp_relaxation_time_taken"])
        d[:cg_sp_total_time_taken] = sum(all_data[size][run_index]["cg_params"]["sp_total_time_taken"])
    else
        d[:cg_number_of_subpaths] = missing
        d[:cg_time_taken] = missing
        d[:cg_lp_relaxation_time_taken] = missing
        d[:cg_sp_total_time_taken] = missing
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
    if "all_subpaths" in keys(all_data[size][run_index])
        d[:enum_number_of_subpaths] = sum(length(v) for v in values(all_data[size][run_index]["all_subpaths"]))
    else
        d[:enum_number_of_subpaths] = missing
    end
    d[:enum_number_of_subpaths] = get(all_data[size][run_index], "enum_number_of_subpaths", missing)
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
CSV.write("$(@__DIR__)/data/all_metrics.csv", all_metrics_df)


all_metrics_df[!, [
    :size, 
    :run_index,
    :arc_ip_objective, 
    :arc_lp_objective,
    :cg_objective,
    :cg_lpip_objective,
    :enum_ip_objective,
    :enum_lp_objective,
]] |>
    x -> permutedims(x, 1, makeunique = true)

all_metrics_df[!, [
    :size, 
    :run_index,
    :arc_ip_time_taken, 
    :arc_lp_time_taken,
    :cg_time_taken,
    :cg_lpip_time_taken,
    :enumerate_all_subpaths_time_taken,
    :enum_ip_time_taken,
    :enum_lp_time_taken,
]] |>
    x -> permutedims(x, 1, makeunique = true)



all_metrics_df[!, [
    :size, 
    :run_index,
    :cg_number_of_subpaths,
    :enum_number_of_subpaths, 
]] |>
    x -> permutedims(x, 1, makeunique = true)




### Scratch work

all_data["xs"][2]["enum_ip_subpaths"] = [
    all_data["xs"][2]["all_subpaths"][k[1]][k[2]]
    for (k, v) in pairs(all_data["xs"][2]["enum_ip_results"]["z"].data)
    if v > 0.5 
]
for s in all_data["xs"][2]["enum_ip_subpaths"]
    v = compute_subpath_reduced_cost(
        s,
        all_data["xs"][2]["data"],
        all_data["xs"][2]["T_range"],
        all_data["xs"][2]["B_range"],
        all_data["xs"][2]["cg_params"]["κ"][end],
        all_data["xs"][2]["cg_params"]["λ"][end],
        all_data["xs"][2]["cg_params"]["μ"][end],
        all_data["xs"][2]["cg_params"]["ν"][end],       
    )
    println(v)
end

s = all_data["xs"][2]["enum_ip_subpaths"][4]

labels = find_smallest_reduced_cost_subpaths(
    20,
    0.0,
    800.0,
    all_data["xs"][2]["G"],
    all_data["xs"][2]["data"],
    all_data["xs"][2]["T_range"],
    all_data["xs"][2]["B_range"],
    all_data["xs"][2]["cg_params"]["κ"][end],
    all_data["xs"][2]["cg_params"]["λ"][end],
    all_data["xs"][2]["cg_params"]["μ"][end],
    all_data["xs"][2]["cg_params"]["ν"][end],
)
Dict(
    k => v.cost
    for (k, v) in pairs(labels)
)

starting_node = 20
starting_time = 0.0
starting_charge = 800.0
G = all_data["xs"][2]["G"]
data = all_data["xs"][2]["data"]
T_range = all_data["xs"][2]["T_range"]
B_range = all_data["xs"][2]["B_range"]
κ = all_data["xs"][2]["cg_params"]["κ"][end]
λ = all_data["xs"][2]["cg_params"]["λ"][end]
μ = all_data["xs"][2]["cg_params"]["μ"][end]
ν = all_data["xs"][2]["cg_params"]["ν"][end]

modified_costs = Float64.(copy(data["c"]))
for i in 1:data["n_customers"]
    j = data["n_customers"] + i
    modified_costs[i,j] -= ν[i]
    modified_costs[j,i] -= ν[i]
end

# initialize set of labels
initial_cost = 0.0
if starting_node in data["N_depots"]
    if starting_time == 0.0 && starting_charge == data["B"]
        initial_cost = initial_cost - κ[starting_node]
    end
elseif starting_node in data["N_charging"]
    initial_cost = initial_cost - λ[(starting_node, starting_time, starting_charge)]
end
labels = Dict()
labels[starting_node] = SubpathWithCost(
    cost = initial_cost,
    n_customers = data["n_customers"],
    starting_node = starting_node,
    starting_time = starting_time,
    starting_charge = starting_charge,
)
# initialize queue of unexplored nodes
Q = [starting_node]

labels[20]


begin 
    i = popfirst!(Q)
    # iterate over all out-neighbors of i
    for j in setdiff(outneighbors(G, i), i)
        # feasibility check
        if j in data["N_pickups"]
            feasible_service = !(labels[i].served[j])
        else
            feasible_service = true
        end
        if !feasible_service 
            continue
        end
        current_time = max(
            data["α"][j], 
            labels[i].time + data["t"][i,j],
        )
        current_charge = labels[i].charge - data["q"][i,j]
        feasible_timewindow = (current_time ≤ data["β"][j])
        feasible_charge = (current_charge ≥ 0.0)
        feasible = (
            feasible_timewindow
            && feasible_charge
        )
        if !feasible
            continue
        end

        current_cost = labels[i].cost + modified_costs[i,j]
        if j in data["N_charging"]
            charging_options = generate_charging_options(
                current_time, current_charge, data, T_range, B_range,
            )
            if length(charging_options) == 0
                continue
            end
            (val, ind) = findmin(
                λ[(j, round_time, round_charge)]
                for (_, _, _, _, round_time, round_charge) in charging_options
            )
            (delta_time, delta_charge, end_time, end_charge, round_time, round_charge) = charging_options[ind]
            current_cost = current_cost + val
        else
            if j in data["N_depots"]
                current_cost = current_cost - μ[j]
            end
            delta_time = 0
            delta_charge = 0
            end_time = current_time
            end_charge = current_charge
            round_time = dceil(end_time, T_range)
            round_charge = dfloor(end_charge, B_range)
        end

        update = true
        if j in keys(labels)
            # dont update label if current label is 
            # better than the current cost (if the subpath ends here)
            if labels[j].cost ≤ current_cost
                update = false
            end
        end

        if update
            # update labels
            s_j = copy(labels[i])
            s_j.cost = current_cost
            s_j.current_node = j
            push!(s_j.arcs, (i,j))
            s_j.time = current_time
            s_j.charge = current_charge
            if j in data["N_dropoffs"]
                s_j.served[j-data["n_customers"]] = true
            end
            s_j.delta_time = delta_time
            s_j.delta_charge = delta_charge
            s_j.end_time = end_time
            s_j.end_charge = end_charge
            s_j.round_time = round_time
            s_j.round_charge = round_charge 
            labels[j] = s_j
        end

        # add node j to the list of unexplored nodes
        if !(j in Q) && !(j in union(data["N_charging"], data["N_depots"]))
            push!(Q, j) 
        end
    end
    println(Q)
end

Q

labels[10]
labels[9]
labels[18]
labels[3]
labels[12]
labels[22]
labels[20]

"""
PROBLEM:
shortest path from I -> j' (with charging) that passes through k 
does not imply that I -> k is a shortest path!

Possible idea: 
when keeping track of shortest path I -> k (so far), 
"keep track" of all possible charging options? k -> j -> j'?

Possible idea: transform the base network?

Possible idea: change the order over which we iterate over the neighbors j of current node i?
"""

T_fine = 0.0:1.0:T_range[end]
B_fine = 0.0:1.0:B_range[end]

λ_reach_heatmap = fill(Inf, length(T_fine), length(B_fine))
for (it, t) in enumerate(T_fine), (ib, b) in enumerate(B_fine)
    for (_, _, _, _, rt, rc) in generate_charging_options(
        t, b, data, T_range, B_range;
        require_charge = true,
    )
        if λ_reach_heatmap[it, ib] > λ[(21, rt, rc)]
            λ_reach_heatmap[it, ib] = λ[(21, rt, rc)]
        end
    end
end

heatmap(λ_reach_heatmap, size = (600, 800))

λ_reach = Dict(
    (n, t, b) => Inf
    for n in data["N_charging"], t in T_fine, b in B_fine
)
for t in T_fine, b in B_fine
    charging_options = generate_charging_options(
        t, b, data, T_range, B_range,
    )
    if length(charging_options) == 0
        continue
    end
    for n in data["N_charging"]
        λ_reach[(n, t, b)] = minimum(
            λ[(n, rt, rb)]
            for (_, _, _, _, rt, rb) in charging_options
        )
    end
end

starting_node = 20
starting_time = 0.0
starting_charge = 800.0
G = all_data["xs"][2]["G"]
data = all_data["xs"][2]["data"]
T_range = all_data["xs"][2]["T_range"]
B_range = all_data["xs"][2]["B_range"]
κ = all_data["xs"][2]["cg_params"]["κ"][end]
λ = all_data["xs"][2]["cg_params"]["λ"][end]
μ = all_data["xs"][2]["cg_params"]["μ"][end]
ν = all_data["xs"][2]["cg_params"]["ν"][end]

new_labels = find_smallest_reduced_cost_subpaths(
    starting_node, starting_time, starting_charge,
    G, data, T_range, B_range, 
    κ, λ, μ, ν,
)

generated_subpaths_withcharge, smallest_reduced_costs, sp_max_time_taken = generate_subpaths_withcharge(
    G, data, T_range, B_range, 
    κ, λ, μ, ν,
)
[
    compute_subpath_reduced_cost(s, data, T_range, B_range, κ, λ, μ, ν)
    for v in values(generated_subpaths_withcharge) 
        for s in v
]
generated_subpaths_withcharge



modified_costs = Float64.(copy(data["c"]))
for i in 1:data["n_customers"]
    j = data["n_customers"] + i
    modified_costs[i,j] -= ν[i]
    modified_costs[j,i] -= ν[i]
end

initial_cost = 0.0
if starting_node in data["N_depots"]
    if starting_time == 0.0 && starting_charge == data["B"]
        initial_cost = initial_cost - κ[starting_node]
    end
elseif starting_node in data["N_charging"]
    initial_cost = initial_cost - λ[(starting_node, starting_time, starting_charge)]
end


starting_node = 21
length(T_range)
cost_array = fill(Inf, (length(T_range), length(B_range), data["n_nodes"], length(T_range), length(B_range)));
size(cost_array)
for (it, t) in enumerate(T_range), (ib, b) in enumerate(B_range)
    cost_array[it,ib,starting_node,it,ib] = -λ[(starting_node, t, b)]
end




labels = Dict()
labels[starting_node] = Dict(
    (starting_time, starting_charge) => SubpathWithCost(
        cost = initial_cost,
        n_customers = data["n_customers"],
        starting_node = starting_node,
        starting_time = starting_time,
        starting_charge = starting_charge,
    )
)

labels[20][(0.0, 800.0)]

Q = [starting_node]

begin
    i = popfirst!(Q)
    for j in setdiff(outneighbors(G, i), i)
        for ((now_time, now_charge), s) in pairs(labels[i])
            if j in data["N_pickups"]
                feasible_service = !(s.served[j])
            else
                feasible_service = true
            end
            if !feasible_service 
                continue
            end
            current_time = max(
                data["α"][j], 
                now_time + data["t"][i,j],
            )
            current_charge = now_charge - data["q"][i,j]
            feasible_timewindow = (current_time ≤ data["β"][j])
            feasible_charge = (current_charge ≥ 0.0)
            feasible = (
                feasible_timewindow
                && feasible_charge
            )
            if !feasible
                continue
            end
            
            current_cost = s.cost + modified_costs[i,j]
            if j in data["N_charging"]
                charging_options = generate_charging_options(
                    current_time, current_charge, 
                    data, T_range, B_range,
                )
                if length(charging_options) == 0
                    continue
                end
                (val, ind) = findmin(
                    λ[(j, rt, rc)]
                    for (_, _, _, _, rt, rc) in charging_options
                )
                (delta_time, delta_charge, end_time, end_charge, round_time, round_charge) = charging_options[ind]
                current_cost = current_cost + val
            else
                if j in data["N_depots"]
                    current_cost = current_cost - μ[j]
                end
                delta_time = 0
                delta_charge = 0
                end_time = current_time
                end_charge = current_charge
                round_time = dceil(end_time, T_range)
                round_charge = dceil(end_charge, B_range)
            end

            if j in keys(labels)
                if (current_time, current_charge) in keys(labels[j])
                    if labels[j][(current_time, current_charge)].cost ≤ current_cost
                        println("Current label at $j, $current_time, $current_charge superior.")
                        continue
                    end
                end
            end

            # update labels
            s_j = copy(s)
            s_j.cost = current_cost
            s_j.current_node = j
            push!(s_j.arcs, (i, j))
            s_j.time = current_time
            s_j.charge = current_charge
            if j in data["N_dropoffs"]
                s_j.served[j-data["n_customers"]] = true
            end
            s_j.delta_time = delta_time
            s_j.delta_charge = delta_charge
            s_j.end_time = end_time
            s_j.end_charge = end_charge
            s_j.round_time = round_time
            s_j.round_charge = round_charge
            if !(j in keys(labels))
                labels[j] = Dict(
                    (current_time, current_charge) => s_j
                )
            else
                labels[j][(current_time, current_charge)] = s_j
            end
            println("$j: $current_time, $current_charge: from ($i, $now_time, $now_charge)")
            if !(j in Q) && !(j in union(data["N_charging"], data["N_depots"]))
                push!(Q, j)
            end
        end
    end
end

Q

labels[1]
labels[2]
labels[3]
labels[4]
labels[5]
labels[6]
labels[7]
labels[8]
labels[9]
labels[10]
labels[11]
labels[12]
labels[13]
labels[14]
labels[15]
labels[16]
labels[17]
labels[18]
labels[19]
labels[20]
labels[21]
labels[22]

labels[19]


Base.max(v::Vector{Float64}, c::Float64) = [
    max(x, c) for x in v
]
Base.max(c::Float64, v::Vector{Float64}) = [
    max(c, x) for x in v
]