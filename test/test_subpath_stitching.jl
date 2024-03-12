using Test
using Suppressor
include("../src/path_formulation.jl")

@testset "Subpath stitching" begin

    (xmin, xmax, ymin, ymax) = (0.0, 2.0, 0.0, 2.0)
    (n_depots, n_customers, n_vehicles) = (4, 12, 6)
    n_charging = Int((xmax - xmin + 1)*(ymax - ymin + 1) - 4)
    depot_pattern = "grid"
    customer_pattern = "random_box"
    charging_pattern = "grid_clipped"
    customer_spread = 0.1
    B = 15000
    μ = 5
    k = 3.0
    T = Int(B * k * (μ + 1) / μ)

    travel_cost_coeff = 7
    charge_cost_coeff = 3
    load_scale = 5.0
    load_shape = 20.0
    load_tolerance = 1.3
    batch = 1
    permissiveness = 0.2

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
        seed = 1,
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
    @test isa(data, EVRPData)

    graph = generate_graph_from_data(data)
    @test isa(graph, EVRPGraph)

    method_params = [
        ("benchmark", false, false)
        ("benchmark", false,  true)
        ("benchmark",  true, false)
        ("ours", false, false)
        ("ours", false,  true)
        ("ours",  true, false)
    ]

    @testset for seed in 1:10
        all_objectives = Dict{Tuple, Float64}()
        all_times = Dict{Tuple, Float64}()
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
        for method_param in method_params
            (method, elementary, ngroute) = method_param
            (
                CGLP_all_results, CGIP_all_results, CG_all_params, CG_all_neighborhoods, all_params, printlist, 
                some_paths, model, z, SR3_constraints
            ) = @suppress path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
                data, graph,
                ;
                Env = Gurobi.Env(),
                method = method,
                elementary = elementary,
                ngroute = ngroute,
                ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
                ngroute_neighborhood_depots_size = "small",
                ngroute_neighborhood_charging_size = "small",
                verbose = true,
                use_adaptive_ngroute = false,
                use_SR3_cuts = false,
                use_lmSR3_cuts = false,
                time_limit = 300.0,
            )
            all_objectives[(method_param..., "LP")] = all_params[1]["CGLP_objective"]
            all_objectives[(method_param..., "IP")] = all_params[1]["CGIP_objective"]
            all_times[method_param] = CG_all_params[1]["time_taken"]
        end
        @printf(
            """
            %.4e %.4e    %.4e %.4e     %6.3f s %6.3f s
            %.4e %.4e    %.4e %.4e     %6.3f s %6.3f s
            %.4e %.4e    %.4e %.4e     %6.3f s %6.3f s

            """,
            all_objectives[("benchmark", false, false, "LP")],
            all_objectives[("ours", false, false, "LP")],
            all_objectives[("benchmark", false, false, "IP")],
            all_objectives[("ours", false, false, "IP")],
            all_times[("benchmark", false, false)],
            all_times[("ours", false, false)],

            all_objectives[("benchmark", false, true, "LP")],
            all_objectives[("ours", false, true, "LP")],
            all_objectives[("benchmark", false, true, "IP")],
            all_objectives[("ours", false, true, "IP")],
            all_times[("benchmark", false, true)],
            all_times[("ours", false, true)],

            all_objectives[("benchmark", true, false, "LP")],
            all_objectives[("ours", true, false, "LP")],
            all_objectives[("benchmark", true, false, "IP")],
            all_objectives[("ours", true, false, "IP")],
            all_times[("benchmark", true, false)],
            all_times[("ours", true, false)],
        )
        for method_param in method_params
            @test all_objectives[(method_param..., "LP")] ≤ all_objectives[(method_param..., "IP")]
        end
        @test all_objectives[("benchmark", false, false, "LP")] ≤ all_objectives[("benchmark", false, true, "LP")]
        @test all_objectives[("benchmark", false, true, "LP")] ≤ all_objectives[("benchmark", true, false, "LP")]
        @test all_objectives[("ours", false, false, "LP")] ≤ all_objectives[("ours", false, true, "LP")]
        @test all_objectives[("ours", false, true, "LP")] ≤ all_objectives[("ours", true, false, "LP")]

        @test all_objectives[("benchmark", false, false, "LP")] ≈ all_objectives[("ours", false, false, "LP")]
        @test all_objectives[("benchmark", false, true, "LP")] ≈ all_objectives[("ours", false, true, "LP")]
        @test all_objectives[("benchmark", true, false, "LP")] ≈ all_objectives[("ours", true, false, "LP")]

        # @test all_objectives[("benchmark", false, false, "IP")] ≈ all_objectives[("ours", false, false, "IP")]
        # @test all_objectives[("benchmark", false, true, "IP")] ≈ all_objectives[("ours", false, true, "IP")]
        # @test all_objectives[("benchmark", true, false, "IP")] ≈ all_objectives[("ours", true, false, "IP")]
    end

end;