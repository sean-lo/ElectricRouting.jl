using Pkg
Pkg.activate("$(@__DIR__)/..")

using CSV, DataFrames
using StatsBase: mean

include("utils.jl")
include("read_evrptw_data.jl")


function read_solution(
    sol_filepath::String,
    data::EVRPData,
    graph::EVRPGraph,
)

    # Read text solution 
    lines = readlines(sol_filepath)

    # Retrieve lines that start with "arc"
    arcs_inds_lines = [
        (i, l) for (i, l) in enumerate(lines) if 
        startswith(l, "arcs")
    ]
    arcs_inds = [i[1] for i in arcs_inds_lines]
    counts_inds = arcs_inds .- 2
    counts = [
        match(r"(\d+.\d+), Path\(Subpath\[Subpath:", lines[i])
        for i in counts_inds
    ]
    path_inds = findall(x -> !isnothing(x), counts)
    counts = [
        parse(Float64, x[1]) |> round |> Int
        for x in counts[path_inds]
    ]
    append!(path_inds, length(arcs_inds_lines) + 1)

    # Get all lists of arcs
    arcs_lines = [i[2] for i in arcs_inds_lines]
    arcs_lines = [
        match(r"arcs:\s+(\[.+\])", l)[1]
        for l in arcs_lines
    ]
    arcs = [
        eval(Meta.parse(l)) 
        for l in arcs_lines
    ]

    paths = [
        [arcs[i:i_next-1]]
        for (i, i_next) in zip(path_inds[1:end-1], path_inds[2:end])
    ]
    
    # Aggregate paths via counts 
    all_paths = vcat([repeat(path, count) for (path, count) in zip(paths, counts)]...)
    all_paths = [x for x in all_paths if !(length(x) == 1 && length(x[1]) == 1)]

    # Collect summary statistics 
    n_recharges = sum([length(x) for x in all_paths] .- 1, init=0.0)
    n_vehicles = length(all_paths)
    if length(all_paths) > 0
        path_lengths = [sum(length(y) for y in x) for x in all_paths]
        mean_path_length = mean(path_lengths)
        subpath_lengths = [length(y) for x in all_paths for y in x]
        mean_subpath_length = mean(subpath_lengths)
        ps_lengths = [length(x) for x in all_paths]
        mean_ps_length = mean(ps_lengths)
    else
        mean_path_length = mean_subpath_length = mean_ps_length = 0.0
    end

    # Path costs
    if length(all_paths) > 0
        path_costs = [
            sum(
                sum(data.c[a[1], a[2]] for a in arcs)
                for arcs in path
            )
            for path in all_paths
        ]
        total_cost = sum(path_costs)
        mean_path_cost = mean(path_costs)
    else
        path_costs = Int[]
        total_cost = mean_path_cost = 0.0
    end

    return Dict(
        "n_recharges" => n_recharges,
        "n_vehicles" => n_vehicles,
        "mean_path_length" => mean_path_length,
        "mean_subpath_length" => mean_subpath_length,
        "mean_ps_length" => mean_ps_length,
        "total_cost" => total_cost,
        "mean_path_cost" => mean_path_cost,
        "path_costs" => path_costs,
        "all_paths" => all_paths,
    )
end

# row_index = 1
# experiment_dir = "$(@__DIR__)/../experiments/ours_benchmark/12/"

# sol_filepath = joinpath(experiment_dir, "solutions/$row_index.txt")
# args_df = CSV.read(joinpath(experiment_dir, "args.csv"), DataFrame)
# r = args_df[row_index, :]

# data = read_evrptw_instance(
#     r.instance_filename |> String,
#     r.n_vehicles,
#     1,
#     ;
#     n_charging = r.n_charging,
#     n_customers = r.n_customers,
#     data_dir = "$(@__DIR__)/../data/",
# )
# graph = generate_graph_from_data(data)


# result = read_solution(sol_filepath, data, graph)
