include("electric_trucking.jl")
include("utils.jl")

nocharge_data, charge_data = generate_instance_pair(
    n_depots = 2, 
    n_customers = 5,
    n_charging = 2,
    charging_repeats = 2,
    n_vehicles = 3,
    shrinkage_depots = 1.4,
    shrinkage_charging = 0.6,
    T = 900.0,
    seed = 0,
    B = 700.0,
    μ = 5.0,
)
plot_instance(nocharge_data, false)
plot_instance(charge_data, true)

# Results: without preprocessing arcs
arc_nocharge_results, arc_nocharge_params = arc_formulation(
    nocharge_data, 
    false,
);
arc_charge_results, arc_charge_params = arc_formulation(
    charge_data, 
    true, 
    time_limit = 3600
);

paths = construct_paths(arc_nocharge_results, nocharge_data)
arc_chargesep_results, arc_chargesep_params = arc_formulation(
    charge_data, 
    true, 
    true; 
    paths = paths, 
    time_limit = 3600
);

# Results: with preprocessing arcs
nocharge_pp_data = preprocess_arcs(nocharge_data, false)
arc_nocharge_pp_results, arc_nocharge_pp_params = arc_formulation(
    nocharge_pp_data, 
    false,
);
paths_pp = construct_paths(arc_nocharge_pp_results, nocharge_pp_data)

arc_nocharge_results["objective"] ≈ arc_nocharge_pp_results["objective"]
issetequal(values(paths), values(paths_pp))

charge_pp_data = preprocess_arcs(charge_data, true, true)
arc_charge_pp_results, arc_charge_pp_params = arc_formulation(
    charge_pp_data, 
    true, 
    time_limit = 3600
);
arc_chargesep_pp_results, arc_chargesep_pp_params = arc_formulation(
    charge_pp_data, 
    true, 
    true; 
    paths = paths_pp, 
    time_limit = 3600
);

# Stricter version of arc preprocessing
charge_pp_e_data = preprocess_arcs(charge_data, true, false)
arc_charge_pp_e_results, arc_charge_pp_e_params = arc_formulation(
    charge_pp_e_data, 
    true, 
    time_limit = 3600
);
arc_chargesep_pp_e_results, arc_chargesep_pp_e_params = arc_formulation(
    charge_pp_e_data, 
    true, 
    true; 
    paths = paths_pp, 
    time_limit = 3600
);

# Printouts
results_printout(
    arc_nocharge_results, arc_nocharge_params, 
    nocharge_data, 
    false,
)
results_printout(
    arc_nocharge_pp_results, 
    arc_nocharge_pp_params, 
    nocharge_pp_data, 
    false,
)

results_printout(
    arc_charge_results, 
    arc_charge_params, 
    charge_data, 
    true,
)
results_printout(
    arc_charge_pp_results, 
    arc_charge_pp_params, 
    charge_pp_data, 
    true,
)
results_printout(
    arc_charge_pp_e_results, 
    arc_charge_pp_e_params, 
    charge_pp_e_data, 
    true,
)

results_printout(
    arc_chargesep_results, 
    arc_chargesep_params, 
    charge_data, 
    true,
)
results_printout(
    arc_chargesep_pp_results, 
    arc_chargesep_pp_params, 
    charge_pp_data, 
    true,
)
results_printout(
    arc_chargesep_pp_e_results, 
    arc_chargesep_pp_e_params, 
    charge_pp_e_data, 
    true,
)