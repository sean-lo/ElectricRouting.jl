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


arc_nocharge_results, arc_nocharge_params = arc_formulation(nocharge_data, false);
results_printout(arc_nocharge_results, arc_nocharge_params, nocharge_data, false)

arc_charge_results, arc_charge_params = arc_formulation(charge_data, true, time_limit = 900);
results_printout(arc_charge_results, arc_charge_params, charge_data, true)
