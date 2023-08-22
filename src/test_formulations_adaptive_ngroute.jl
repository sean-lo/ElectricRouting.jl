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
    ngroute_neighborhood_charging_depots_size = String[],
    time_taken = Float64[],
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

for seed in 1:20, (n_customers,) in [
    (16,)
    (20,)
    (24,)
], ngroute_neighborhood_charging_depots_size in ["small", "medium", "large"]
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
        method = "ours",
        ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
        ngroute_neighborhood_charging_depots_size = ngroute_neighborhood_charging_depots_size,
        verbose = true,
    );
    (
        CGLP_all_results, CGIP_all_results, CG_all_params, printlist, 
        some_paths, model, z, neighborhoods, WSR3_constraints
    ) = run.value;

    push!(all_results, 
        (
            seed,
            n_customers,
            Int(ceil(sqrt(graph.n_customers))),
            ngroute_neighborhood_charging_depots_size,
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
    println("""
    $(n_customers)\t$(ngroute_neighborhood_charging_depots_size)\t$(seed)
    Completed in $(run.time) s.
    Initial gap $(CG_all_params[1]["LP_IP_gap"]), final gap $(CG_all_params[end]["LP_IP_gap"]).
    """)
end

all_results |>
    x -> sort!(x, [
        order(:n_customers), 
        order(:ngroute_neighborhood_charging_depots_size, rev = true),
    ])
all_results.LP_IP_gap_first .*= 100
all_results.LP_IP_gap_last .*= 100
CSV.write("adaptive_ngroute.csv", all_results)

(
    all_results
    |> x -> groupby(x, [:n_customers, :ngroute_neighborhood_charging_depots_size])
    |> x -> combine(
        x, 
        :time_taken => geomean,
        :sp_time_taken_mean_first => geomean,
        :sp_time_taken_mean_last => geomean,
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
        :ngroute_neighborhood_charging_depots_size, 
        :LP_IP_gap_last
    )
    |> x -> filter(r -> r.n_customers == 24, x)
)