include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools
using Cthulhu
using Test

### Experiment

all_results = DataFrame(
    seed = Int[],
    n_customers = Int[],
    ngroute_neighborhood_size = Int[],
    ngroute_neighborhood_depots_size = String[],
    ngroute_neighborhood_charging_size = String[],
    method = String[],
    ngroute_alt = Bool[],
    time_taken_first = Float64[],
    time_taken_total = Float64[],
    time_taken_total_SR3 = Float64[],
    n_iterations = Int[],
    n_CG_iterations = Int[],
    sp_time_taken_mean_first = Float64[],
    sp_time_taken_mean_last = Float64[],
    sp_time_taken_mean_last_SR3 = Float64[],
    LP_objective_first = Float64[],
    LP_objective_last = Float64[],
    LP_objective_last_SR3 = Float64[],
    IP_objective_first = Float64[],
    IP_objective_last = Float64[],
    IP_objective_last_SR3 = Float64[],
    LP_IP_gap_first = Float64[],
    LP_IP_gap_last = Float64[],
    LP_IP_gap_last_SR3 = Float64[],
)

# compilation
for n_customers in [16, 20, 24],
    (method, ngroute_alt) in [
        ("benchmark", false),
        ("benchmark", true),
        ("ours", false),
        ("ours", true),
    ],
    seed in [1]
    data = generate_instance(
        ;
        n_depots = 4,
        n_customers = n_customers,
        n_charging = 7,
        n_vehicles = 6,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "circular_packing",
        shrinkage_depots = 1.0,
        shrinkage_charging = 1.0,
        T = 40000,
        seed = seed,
        B = 15000,
        μ = 5,
        travel_cost_coeff = 7,
        charge_cost_coeff = 3,
        load_scale = 5.0,
        load_shape = 20.0,
        load_tolerance = 1.3,
        batch = 1,
        permissiveness = 0.2,
    )
    graph = generate_graph_from_data(data)
    path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        ngroute_alt = ngroute_alt,
        verbose = true,
        use_adaptive_ngroute = true,
        use_SR3_cuts = true,
    );
end

(n_customers, seed, ngroute_alt, ngroute_neighborhood_charging_size, method) = (20, 3, true, "small", "ours")
data = generate_instance(
    ;
    n_depots = 4,
    n_customers = n_customers,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = seed,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)
r = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = method,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    ngroute_alt = ngroute_alt,
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    max_SR3_cuts = 10,
);
DataFrame(r[5])

for (
    ngroute_alt,
    ngroute_neighborhood_depots_size, 
    ngroute_neighborhood_charging_size,
    seed, 
    n_customers, 
    method, 
) in Iterators.product(
    [true, false],
    ["small",],
    # ["small",],
    ["small", "medium"],
    # [1],
    1:20,
    # [16],
    [16, 20, 24],
    ["ours", "benchmark",],
)
    data = generate_instance(
        ;
        n_depots = 4,
        n_customers = n_customers,
        n_charging = 7,
        n_vehicles = 6,
        depot_pattern = "circular",    
        customer_pattern = "random_box",
        charging_pattern = "circular_packing",
        shrinkage_depots = 1.0,
        shrinkage_charging = 1.0,
        T = 40000,
        seed = seed,
        B = 15000,
        μ = 5,
        travel_cost_coeff = 7,
        charge_cost_coeff = 3,
        load_scale = 5.0,
        load_shape = 20.0,
        load_tolerance = 1.3,
        batch = 1,
        permissiveness = 0.2,
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
        max_SR3_cuts = 10,
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

    push!(all_results, 
        (
            seed,
            n_customers,
            Int(ceil(sqrt(graph.n_customers))),
            ngroute_neighborhood_depots_size,
            ngroute_neighborhood_charging_size,
            method,
            ngroute_alt,
            CG_all_params[1]["time_taken"],
            sum(all_params_df.CG_time_taken[1:ind]),
            sum(all_params_df.CG_time_taken),
            length(CG_all_params),
            sum(CG_params["counter"] for CG_params in CG_all_params),
            CG_all_params[1]["sp_time_taken_mean"],
            CG_all_params[ind]["sp_time_taken_mean"],
            CG_all_params[end]["sp_time_taken_mean"],
            all_params_df.CGLP_objective[1],
            all_params_df.CGLP_objective[ind],
            all_params_df.CGLP_objective[end],
            all_params_df.CGIP_objective[1],
            all_params_df.CGIP_objective[ind],
            all_params_df.CGIP_objective[end],
            all_params_df.CG_LP_IP_gap[1],
            all_params_df.CG_LP_IP_gap[ind],
            all_params_df.CG_LP_IP_gap[end],
        )
    )
    # println("""
    # $method\t$ngroute_alt\t$(n_customers)
    # $(ngroute_neighborhood_depots_size)\t$(ngroute_neighborhood_charging_size)\t$(seed)
    # Completed in $(run.time) s.
    # $(CGLP_all_results[1]["objective"]),\t$(CGLP_all_results[end]["objective"])
    # $(CGIP_all_results[1]["objective"]),\t$(CGIP_all_results[end]["objective"])
    # """)
    @printf("""
    %s    \t%s\t%2d\t%s\t%s\t%2d : \t%6.1f s.

    %9.2f\t%9.2f\t%9.2f
    %9.2f\t%9.2f\t%9.2f

    """,
        method, ngroute_alt, n_customers, 
        ngroute_neighborhood_depots_size, ngroute_neighborhood_charging_size, seed,
        run.time,
        all_params_df.CGLP_objective[1], all_params_df.CGLP_objective[ind], all_params_df.CGLP_objective[end],
        all_params_df.CGIP_objective[1], all_params_df.CGIP_objective[ind], all_params_df.CGIP_objective[end],
    )
end

all_results |>
    x -> sort!(x, [
        order(:n_customers), 
        order(:ngroute_neighborhood_depots_size, rev = true),
        order(:ngroute_neighborhood_charging_size, rev = true),
    ])
all_results.LP_IP_gap_first .*= 100
all_results.LP_IP_gap_last .*= 100
CSV.write("adaptive_ngroute_SR3.csv", all_results)
all_results_SR3 = CSV.read("adaptive_ngroute_SR3.csv", DataFrame)

all_results_SR3[!, "converged_first"] .= (
    all_results_SR3[!, "LP_IP_gap_first"]
    .≈ 0.0
)
all_results_SR3[!, "converged_last"] .= (
    all_results_SR3[!, "LP_IP_gap_last"]
    .≈ 0.0
)
all_results_SR3[!, "converged_last_SR3"] .= (
    all_results_SR3[!, "LP_IP_gap_last_SR3"]
    .≈ 0.0
)
summary_SR3 = (
    all_results_SR3
    |> x -> groupby(x, [
        :method, 
        :n_customers,
        :ngroute_alt,
        :ngroute_neighborhood_charging_size,
    ])
    |> x -> combine(
        x, 
        :time_taken_first => geomean,
        :time_taken_total => geomean,
        :time_taken_total_SR3 => geomean,
        :sp_time_taken_mean_first => geomean,
        :sp_time_taken_mean_last => geomean,
        :sp_time_taken_mean_last_SR3 => geomean,
        :LP_IP_gap_first => mean,
        :LP_IP_gap_last => mean,
        :LP_IP_gap_last_SR3 => mean,
        :converged_first => mean,
        :converged_last => mean,
        :converged_last_SR3 => mean,
        :n_CG_iterations => mean,
        :n_iterations => mean,
    )
)
(
    summary_SR3
    |> x -> unstack(
        x, 
        [:method, :ngroute_alt, :ngroute_neighborhood_charging_size,],
        :n_customers,
        :time_taken_total_geomean,
    )
)
(
    summary_SR3
    |> x -> unstack(
        x, 
        [:method, :ngroute_alt, :ngroute_neighborhood_charging_size,],
        :n_customers,
        :time_taken_total_SR3_geomean,
    )
)

# How the % converged increases as you add adaptive NG-route, and SR3 cuts
(
    summary_SR3
    |> x -> unstack(
        x, 
        [:method, :ngroute_alt, :ngroute_neighborhood_charging_size,],
        :n_customers,
        :converged_first_mean,
    )
)
(
    summary_SR3
    |> x -> unstack(
        x, 
        [:method, :ngroute_alt, :ngroute_neighborhood_charging_size,],
        :n_customers,
        :converged_last_mean,
    )
)
(
    summary_SR3
    |> x -> unstack(
        x, 
        [:method, :ngroute_alt, :ngroute_neighborhood_charging_size,],
        :n_customers,
        :converged_last_SR3_mean,
    )
)

(
    all_results 
    |> x -> unstack(
        x, 
        [:n_customers, :seed], 
        :ngroute_neighborhood_depots_size, 
        :ngroute_neighborhood_charging_size, 
        :LP_objective_last
    )
    |> x -> filter(r -> r.n_customers == 24, x)
)


### Testing

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 16,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 2,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)

for (method, ngroute_alt) in [
    ("benchmark", false),
    ("benchmark", true),
    ("ours", false),
    ("ours", true),
]
    path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        ngroute_alt = ngroute_alt,
        verbose = true,
        use_adaptive_ngroute = true,
        use_SR3_cuts = true,
    );
end

data = generate_instance(
    ;
    n_depots = 4,
    n_customers = 16,
    n_charging = 7,
    n_vehicles = 6,
    depot_pattern = "circular",    
    customer_pattern = "random_box",
    charging_pattern = "circular_packing",
    shrinkage_depots = 1.0,
    shrinkage_charging = 1.0,
    T = 40000,
    seed = 1,
    B = 15000,
    μ = 5,
    travel_cost_coeff = 7,
    charge_cost_coeff = 3,
    load_scale = 5.0,
    load_shape = 20.0,
    load_tolerance = 1.3,
    batch = 1,
    permissiveness = 0.2,
)
graph = generate_graph_from_data(data)

ours_result = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_alt = false,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
);
(
    ours_CGLP_all_results, ours_CGIP_all_results, ours_CG_all_params, ours_CG_all_neighborhoods, ours_all_params, 
    ours_printlist, ours_some_paths, ours_model, ours_z, ours_SR3_constraints
) = ours_result.value;

findfirst(x -> x == "use_SR3_cuts", DataFrame(ours_all_params).method)

plot_path_solution(
    ours_CGLP_all_results[end], data, graph, ours_some_paths
)


ours_alt_result = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_alt = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
);
(
    ours_alt_CGLP_all_results, ours_alt_CGIP_all_results, ours_alt_CG_all_params, ours_alt_printlist, 
    ours_alt_some_paths, ours_alt_model, ours_alt_z, ours_alt_CG_all_neighborhoods, ours_alt_WSR3_constraints
) = ours_alt_result.value;

plot_path_solution(
    ours_CGLP_all_results[end], data, graph, ours_some_paths
)

enumerate_violated_path_SRnk_inequalities(
    ours_CGLP_all_results[end]["paths"],
    graph,
    4, 3
)

for (method, ngroute_alt) in [
    ("benchmark", false),
    ("benchmark", true),
    ("ours", false),
    ("ours", true),
]
    @time @suppress path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
        data, graph,
        ;
        method = method,
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_depots_size = "small", 
        ngroute_neighborhood_charging_size = "small", 
        ngroute_alt = ngroute_alt,
        verbose = true,
        use_adaptive_ngroute = true,
        use_SR3_cuts = false,
    );
end
# todo: debug difference between ours true and everything else

ours_alt_results = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    ngroute_alt = true,
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
);

(
    ours_alt_CGLP_all_results, ours_alt_CGIP_all_results, ours_alt_CG_all_params, ours_alt_printlist, 
    ours_alt_some_paths, ours_alt_model, ours_alt_z, ours_alt_CG_all_neighborhoods, ours_alt_WSR3_constraints
) = ours_alt_results.value;

ours_results = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    ngroute_alt = false,
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
);
(
    ours_CGLP_all_results, ours_CGIP_all_results, ours_CG_all_params, ours_printlist, 
    ours_some_paths, ours_model, ours_z, ours_CG_all_neighborhoods, ours_WSR3_constraints
) = ours_results.value;

ours_alt_results.time
ours_results.time


[
    x["objective"] for x in 
    ours_CGLP_all_results
]
[
    x["objective"] for x in 
    ours_CGIP_all_results
]

[
    x["objective"] for x in 
    ours_alt_CGLP_all_results
]
[
    x["objective"] for x in 
    ours_alt_CGIP_all_results
]

ours_CG_all_params[1]["number_of_paths"]
ours_alt_CG_all_params[1]["number_of_paths"]

[
    compute_path_modified_cost(
        data, graph, p, 
        ours_CG_all_params[2]["κ"][end],
        ours_CG_all_params[2]["μ"][end],
        ours_CG_all_params[2]["ν"][end],
    )
    for (val, p) in ours_alt_CGLP_all_results[2]["paths"]
]
p = ours_alt_CGLP_all_results[2]["paths"][1][2]

compute_path_modified_cost(
    data, graph, p, 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
    ;
    verbose = true,
)

base_labels = @time generate_base_labels_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)
full_labels = @time find_nondominated_paths_notimewindows_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    base_labels, 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ;
)


using Cthulhu
using BenchmarkTools
include("subpath_stitching.jl")
include("utils.jl")

base_labels_ngroute = @btime generate_base_labels_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)

@descend generate_base_labels_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)

base_labels_ngroute_alt = @btime generate_base_labels_ngroute_alt(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)

@descend generate_base_labels_ngroute_alt(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)

full_labels_ngroute = @btime find_nondominated_paths_notimewindows_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    base_labels_ngroute,
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
)

@descend find_nondominated_paths_notimewindows_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    base_labels_ngroute,
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
)

full_labels_ngroute_alt = @btime find_nondominated_paths_notimewindows_ngroute_alt(
    data, graph, ours_CG_all_neighborhoods[2], 
    base_labels_ngroute_alt,
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
)

@descend find_nondominated_paths_notimewindows_ngroute_alt(
    data, graph, ours_CG_all_neighborhoods[2], 
    base_labels_ngroute_alt,
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
)


include("desaulniers_benchmark.jl")

α = zeros(Int, graph.n_nodes)
β = fill(graph.T, graph.n_nodes)
@btime find_nondominated_paths_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
    α, β,
)

@descend find_nondominated_paths_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
    α, β,
)

@btime find_nondominated_paths_ngroute_alt(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
    α, β,
)

@descend find_nondominated_paths_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
    α, β,
)

include("utils.jl")
include("subpath_stitching.jl")

neighborhoods = ours_CG_all_neighborhoods[2]
nodes = [27, 11, 3, 10, 27]
set = falses(graph.n_nodes)
set[27] = 1
@btime ngroute_check_create_fset(
    neighborhoods,
    set,
    nodes[4],
)

set = falses(graph.n_nodes)
set[27] = 1
@btime ngroute_create_bset(
    $neighborhoods,
    $nodes[1:4],
    $set,
)

findall(set)

begin
    set = falses(graph.n_nodes)
    set[27] = 1
    for i in 2:length(nodes)
        (f, set) = ngroute_check_create_fset(
            neighborhoods, 
            set, 
            nodes[i],
        )
    end
    subpath_fset = set
    println(findall(subpath_fset))
end
begin
    set = falses(graph.n_nodes)
    set[27] = 1
    for i in 2:length(nodes)
        set = ngroute_create_bset(
            neighborhoods,
            nodes[1:i],
            set,
        )
    end
    subpath_bset = set
    println(findall(subpath_bset))
end
path_fset = collect(keys(base_labels_ngroute[(17,27)]))[3]
findall(path_fset)


(f1, s1) = @btime ngroute_extend_partial_path_check(
    $neighborhoods,
    $path_fset,
    $nodes,
)


s1 == s2
all(subpath_fset .≤ s1)
all(s1 .≤ subpath_fset .| path_fset)

findall(subpath_fset)
findall(s1)
findall(subpath_fset .| path_fset)

findall(s2)

neighborhoods

compute_subpath_modified_cost(
    data, graph, p.subpaths[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)
p.subpaths[2].current_time - p.subpaths[2].starting_time

base_labels_ngroute[(27,27)]

base_labels_ngroute[(17,27)]
base_labels_ngroute[(27,26)]
base_labels_ngroute[(26,18)]

base_labels[17][27][(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)]
base_labels[27][26]
base_labels[27][26][(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)]
base_labels[26][18][(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

base_labels[27][27]
base_labels[27][26]
keys(base_labels[27][27])
neighborhoods = ours_CG_all_neighborhoods[2]


findall(neighborhoods[:,27])
# 3, 10, 16, 27
findall(neighborhoods[:,9])
# 3, 9, 10, 16
findall(neighborhoods[:,11])
# 11, 16
findall(neighborhoods[:,6])
# 6, 11
findall(neighborhoods[:,26])
# 26


findall(neighborhoods[:,16])

findall(neighborhoods[:,10])

