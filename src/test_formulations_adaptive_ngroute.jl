include("arc_formulation.jl")
include("subpath_formulation.jl")
include("path_formulation.jl")
include("utils.jl")

using CSV, DataFrames
using BenchmarkTools

using Test

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
    n_iterations = Int[],
    n_CG_iterations = Int[],
    sp_time_taken_mean_first = Float64[],
    sp_time_taken_mean_last = Float64[],
    LP_objective_first = Float64[],
    LP_objective_last = Float64[],
    IP_objective_first = Float64[],
    IP_objective_last = Float64[],
    LP_IP_gap_first = Float64[],
    LP_IP_gap_last = Float64[],
)

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
    ["small", "medium", "large"],
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
        use_SR3_cuts = false,
    );
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, printlist, 
        some_paths, model, z, CG_all_neighborhoods, WSR3_constraints
    ) = run.value;

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
            run.time,
            length(CG_all_params),
            sum(CG_params["counter"] for CG_params in CG_all_params),
            CG_all_params[1]["sp_time_taken_mean"],
            CG_all_params[end]["sp_time_taken_mean"],
            CGLP_all_results[1]["objective"],
            CGLP_all_results[end]["objective"],
            CGIP_all_results[1]["objective"],
            CGIP_all_results[end]["objective"],
            CG_all_params[1]["LP_IP_gap"],
            CG_all_params[end]["LP_IP_gap"],
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

    %9.2f\t%9.2f
    %9.2f\t%9.2f

    """,
        method, ngroute_alt, n_customers, 
        ngroute_neighborhood_depots_size, ngroute_neighborhood_charging_size, seed,
        run.time,
        CGLP_all_results[1]["objective"], CGLP_all_results[end]["objective"],
        CGIP_all_results[1]["objective"], CGIP_all_results[end]["objective"],
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
CSV.write("adaptive_ngroute.csv", all_results)

all_results = CSV.read("adaptive_ngroute.csv", DataFrame)

# Sensemaking

for group in (
    all_results
    |> x -> stack(
        x, 
        [:LP_objective_first, :LP_objective_last, :IP_objective_first, :IP_objective_last,],
        [:seed, :n_customers, :method, :ngroute_alt, :ngroute_neighborhood_charging_size,]
    )
    |> x -> sort(
        x, 
        [:seed, :n_customers, :method, :ngroute_alt, :ngroute_neighborhood_charging_size,]
    )
    |> x -> groupby(
        x, 
        [:seed, :n_customers,]
    )
)
    if (group.n_customers[1], group.seed[1]) in [(16, 7), (24, 7), (16, 17)]
        continue
    end
    @test all(
        filter(r -> r.variable == "LP_objective_first", group).value
        .≤ filter(r -> r.variable == "LP_objective_last", group).value .+ 1e-6
    )
    @test all(
        filter(r -> r.variable == "IP_objective_first", group).value
        .≥ filter(r -> r.variable == "IP_objective_last", group).value .- 1e-6
    )
    @test all(
        filter(r -> r.variable == "LP_objective_first", group).value
        .≤ filter(r -> r.variable == "IP_objective_first", group).value .+ 1e-6
    )
    @test all(
        filter(r -> r.variable == "LP_objective_last", group).value
        .≤ filter(r -> r.variable == "IP_objective_last", group).value .+ 1e-6
    )
    group_lp_first = filter(r -> r.variable == "LP_objective_first", group)
    @test all(
        filter(r -> r.ngroute_neighborhood_charging_size == "small", group_lp_first).value
        .≤ filter(r -> r.ngroute_neighborhood_charging_size == "medium", group_lp_first).value .+ 1e-6
        .≤ filter(r -> r.ngroute_neighborhood_charging_size == "large", group_lp_first).value .+ 2e-6
    )
    @test all(
        filter(r -> r.method == "ours", group_lp_first).value
        .≈ filter(r -> r.method == "benchmark", group_lp_first).value
    )
    @test all(
        filter(r -> r.ngroute_alt == true, group_lp_first).value
        .≈ filter(r -> r.ngroute_alt == false, group_lp_first).value
    )
    group_lp_last = filter(r -> r.variable == "LP_objective_last", group)
    @test all(
        filter(r -> r.ngroute_neighborhood_charging_size == "small", group_lp_last).value
        .≈ filter(r -> r.ngroute_neighborhood_charging_size == "medium", group_lp_last).value
        .≈ filter(r -> r.ngroute_neighborhood_charging_size == "large", group_lp_last).value
    )
    @test all(
        filter(r -> r.method == "ours", group_lp_last).value
        .≈ filter(r -> r.method == "benchmark", group_lp_last).value
    )
    @test all(
        filter(r -> r.ngroute_alt == true, group_lp_last).value
        .≈ filter(r -> r.ngroute_alt == false, group_lp_last).value
    )
    println("$(group.n_customers[1]), $(group.seed[1])")
end

(
    all_results
    |> x -> groupby(x, [
        :method, 
        :n_customers,
        :ngroute_alt,
        :ngroute_neighborhood_charging_size,
    ])
    |> x -> combine(
        x, 
        # :time_taken_first => geomean,
        # :time_taken_total => geomean,
        # :sp_time_taken_mean_first => geomean,
        # :sp_time_taken_mean_last => geomean,
        :LP_IP_gap_first => mean,
        :LP_IP_gap_last => mean,
        :n_CG_iterations => mean,
        :n_iterations => mean,
    )
)

(
    all_results 
    |> x -> unstack(
        x, 
        [:n_customers, :seed], 
        :ngroute_neighborhood_depots_size, 
        :ngroute_neighborhood_charging_size, 
        :LP_IP_gap_last
    )
    |> x -> filter(r -> r.n_customers == 24, x)
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
[
    (16, 1),
    (16, 12),
    (20, 16),
    (20, 18),
    (24, 12),
]

for row in all_results |>
    x -> unstack(
        x, 
        [:n_customers, :seed],
        :ngroute_neighborhood_depots_size, 
        :ngroute_neighborhood_charging_size, 
        :LP_objective_last,
    ) |> eachrow
    if row.small ≈ row.medium ≈ row.large
        continue
    end
    println(
        "$(row.n_customers), $(row.seed), $(row.small), $(row.medium), $(row.large)"
    )
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
    seed = 12,
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

s_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
);
(
    s_CGLP_all_results, s_CGIP_all_results, s_CG_all_params, s_printlist, 
    s_some_paths, s_model, s_z, s_CG_all_neighborhoods, s_WSR3_constraints
) = s_run.value;

s_CG_all_params[1]["time_taken"]

m_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "medium", 
    ngroute_neighborhood_charging_size = "medium", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
);
(
    m_CGLP_all_results, m_CGIP_all_results, m_CG_all_params, m_printlist, 
    m_some_paths, m_model, m_z, m_CG_all_neighborhoods, m_WSR3_constraints
) = m_run.value;


l_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "large", 
    ngroute_neighborhood_charging_size = "large", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
);
(
    l_CGLP_all_results, l_CGIP_all_results, l_CG_all_params, l_printlist, 
    l_some_paths, l_model, l_z, l_CG_all_neighborhoods, l_WSR3_constraints
) = l_run.value;

ll_run = @timed path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    ngroute_neighborhood_size = graph.n_customers,
    ngroute_neighborhood_depots_size = "large", 
    ngroute_neighborhood_charging_size = "large", 
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = false,
);
(
    ll_CGLP_all_results, ll_CGIP_all_results, ll_CG_all_params, ll_printlist, 
    ll_some_paths, ll_model, ll_z, ll_CG_all_neighborhoods, ll_WSR3_constraints
) = ll_run.value;

ll_CGLP_all_results[end]["paths"] = collect_path_solution_support(
    ll_CGLP_all_results[end], ll_some_paths, data, graph
)

[
    x["objective"] for x in 
    s_CGLP_all_results
]
[
    x["objective"] for x in 
    m_CGLP_all_results
]
[
    x["objective"] for x in 
    l_CGLP_all_results
]
[
    x["objective"] for x in 
    ll_CGLP_all_results
]

[
    x["objective"] for x in 
    s_CGIP_all_results
]
[
    x["objective"] for x in 
    m_CGIP_all_results
]
[
    x["objective"] for x in 
    l_CGIP_all_results
]
[
    x["objective"] for x in 
    ll_CGIP_all_results
]

for (ind, (val, p)) in enumerate(ll_CGLP_all_results[end]["paths"])
    cost = compute_path_modified_cost(
        data, graph, p,
        s_CGLP_all_results[end]["κ"],
        s_CGLP_all_results[end]["μ"],
        s_CGLP_all_results[end]["ν"],
        ;
    )
    if cost < 0
        println(ind)
        compute_path_modified_cost(
            data, graph, p,
            s_CGLP_all_results[end]["κ"],
            s_CGLP_all_results[end]["μ"],
            s_CGLP_all_results[end]["ν"],
            ;
            verbose = true,
        )
    end   
end

plot_path_solution(
    ll_CGLP_all_results[end], data, graph, ll_some_paths
)

plot_path_solution(
    l_CGLP_all_results[end], data, graph, l_some_paths
)

enumerate_violated_path_SR3_inequalities(
    l_CGLP_all_results[end]["paths"], graph,
)

p = ll_CGLP_all_results[end]["paths"][1][2]

s_CGLP_all_results[end]["paths"][1][2]

compute_path_modified_cost(
    data, graph, p,
    s_CGLP_all_results[end]["κ"],
    s_CGLP_all_results[end]["μ"],
    s_CGLP_all_results[end]["ν"],
    verbose = true,
)

base_labels = generate_base_labels_ngroute(
    data, graph,
    s_CG_all_neighborhoods[end], 
    s_CGLP_all_results[end]["κ"],
    s_CGLP_all_results[end]["μ"],
    s_CGLP_all_results[end]["ν"],
)


base_labels[17][27][(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)]
base_labels[27][26]
base_labels[26][18]





base_labels[17][27][(0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)]
base_labels[27][26][(0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)]
base_labels[26][18][(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

data.N_depots
findall(s_CG_all_neighborhoods[end][17,:])
findall(s_CG_all_neighborhoods[end][18,:])
findall(s_CG_all_neighborhoods[end][19,:])
findall(s_CG_all_neighborhoods[end][20,:])

findall(m_CG_all_neighborhoods[end][17,:])
findall(m_CG_all_neighborhoods[end][18,:])
findall(m_CG_all_neighborhoods[end][19,:])
findall(m_CG_all_neighborhoods[end][20,:])

findall(l_CG_all_neighborhoods[end][17,:])
findall(l_CG_all_neighborhoods[end][18,:])
findall(l_CG_all_neighborhoods[end][19,:])
findall(l_CG_all_neighborhoods[end][20,:])

findall(s_CG_all_neighborhoods[end][17,:])
findall(s_CG_all_neighborhoods[end][18,:])
findall(s_CG_all_neighborhoods[end][19,:])
findall(s_CG_all_neighborhoods[end][20,:])


for (method, ngroute_alt) in [
    # ("benchmark", false),
    # ("benchmark", true),
    # ("ours", false),
    ("ours", true),
]
    @time path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
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

base_labels = generate_base_labels_ngroute(
    data, graph, ours_CG_all_neighborhoods[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)

p
compute_subpath_modified_cost(
    data, graph, p.subpaths[2], 
    ours_CG_all_params[2]["κ"][end],
    ours_CG_all_params[2]["μ"][end],
    ours_CG_all_params[2]["ν"][end],
)
p.subpaths[2].current_time - p.subpaths[2].starting_time


base_labels[17][27]
base_labels[27][26]
base_labels[26][18]

base_labels[17][27][(0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1)]
base_labels[27][26]
base_labels[27][26][(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0)]
base_labels[26][18][(0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0)]

base_labels[27][27]
base_labels[27][26]
keys(base_labels[27][27])
neighborhoods = ours_CG_all_neighborhoods[2]


findall(neighborhoods[27,:])
# 3, 10, 16, 27
findall(neighborhoods[9,:])
# 3, 9, 10, 16
findall(neighborhoods[11,:])
# 11, 16
findall(neighborhoods[6,:])
# 6, 11
findall(neighborhoods[26,:])
# 26


findall(neighborhoods[16,:])

findall(neighborhoods[10,:])

