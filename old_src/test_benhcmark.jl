using Pkg
Pkg.activate("$(@__DIR__)/..")

include("utils.jl")

using CSV, DataFrames

args_df = CSV.read("$(@__DIR__)/../experiments/ours_benchmark/01a/args.csv", DataFrame)
row_index = 5000

n_depots = args_df[row_index, :n_depots]
n_customers = args_df[row_index, :n_customers]
# n_customers = 20
n_charging = args_df[row_index, :n_charging]
depot_pattern = String(args_df[row_index, :depot_pattern])
customer_pattern = String(args_df[row_index, :customer_pattern])
charging_pattern = String(args_df[row_index, :charging_pattern])
customer_spread = args_df[row_index, :customer_spread]
xmin = args_df[row_index, :xmin]
xmax = args_df[row_index, :xmax]
ymin = args_df[row_index, :ymin]
ymax = args_df[row_index, :ymax]
n_vehicles = args_df[row_index, :n_vehicles]
T = args_df[row_index, :T]
B = args_df[row_index, :B]
inverse_refueling_rate = args_df[row_index, :μ] |> Float64
seed = args_df[row_index, :seed]
travel_cost_coeff = args_df[row_index, :travel_cost_coeff]
charge_cost_coeff = args_df[row_index, :charge_cost_coeff]
charge_cost_coeff = 3
load_scale = args_df[row_index, :load_scale]
load_shape = args_df[row_index, :load_shape]
load_tolerance = args_df[row_index, :load_tolerance]
batch = args_df[row_index, :batch]
permissiveness = args_df[row_index, :permissiveness]

use_load = args_df[row_index, :use_load]
use_time_windows = args_df[row_index, :use_time_windows]

# use_lmSR3_cuts = args_df[row_index, :use_lmSR3_cuts]
# max_SR3_cuts = args_df[row_index, :max_SR3_cuts]
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
    inverse_refueling_rate = inverse_refueling_rate,
    travel_cost_coeff = travel_cost_coeff,
    charge_cost_coeff = charge_cost_coeff,
    load_scale = load_scale,
    load_shape = load_shape,
    load_tolerance = load_tolerance,
    batch = batch,
    permissiveness = permissiveness * 2.0,
    ;
    data_dir = "../../../data/",
)
graph = generate_graph_from_data(data)

plot_instance(data, plot_edges = true, graph = graph)



include("path_formulation.jl")
include("subpath_stitching_v2.jl")
include("desaulniers_benchmark.jl")

@timev ours_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "ours",
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = true,
    max_SR3_cuts = 5,
    time_limit = 120.0,
);



@timev benchmark_results = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "benchmark",
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
    time_limit = 120.0,
);



@timev benchmark_results_tw = path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts(
    data, graph,
    ;
    method = "benchmark",
    time_windows = true,
    elementary = false,
    ngroute = true,
    ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers))),
    ngroute_neighborhood_depots_size = "small", 
    ngroute_neighborhood_charging_size = "small",
    verbose = true,
    use_adaptive_ngroute = true,
    use_SR3_cuts = true,
    use_lmSR3_cuts = false,
    max_SR3_cuts = 5,
    time_limit = 120.0,
);


include("path_formulation.jl")

### Start debug"
Env = Gurobi.Env()
method = "ours"
charge_cost_heterogenous = false
time_windows = use_time_windows
elementary = false
ngroute = true
verbose = true
ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers)))
ngroute_neighborhood_depots_size = "small"
ngroute_neighborhood_charging_size = "small"
max_SR3_cuts = 5

neighborhoods = compute_ngroute_neighborhoods(
    graph,
    ngroute_neighborhood_size; 
    depots_size = ngroute_neighborhood_depots_size,
    charging_size = ngroute_neighborhood_charging_size,
)

artificial_paths = generate_artificial_paths(data, graph)
some_paths = deepcopy(artificial_paths)
path_costs = compute_path_costs(
    data, graph, 
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)

model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
    ; 
    Env = Env,
)
SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
SR3_list = Tuple{Float64, NTuple{3, Int}}[]
printlist = String[]

# include("path_formulation.jl")
# CGLP_results, CG_params = path_formulation_column_generation!(
#     model, z, SR3_constraints,
#     data, graph,
#     artificial_paths,
#     some_paths, path_costs, path_service,
#     printlist,
#     ;
#     method = method,
#     charge_cost_heterogenous = charge_cost_heterogenous,
#     time_windows = time_windows,
#     elementary = elementary,
#     neighborhoods = neighborhoods,
#     ngroute = ngroute,
#     verbose = verbose,
#     time_limit = 120.0,
#     max_iters = Inf,
# )


CGLP_all_results = Dict[]
CGIP_all_results = Dict[]
all_params = Dict[]
CG_all_params = Dict[]
if ngroute
    CG_all_neighborhoods = BitMatrix[]
else
    CG_all_neighborhoods = nothing
end



# start of loop
CGLP_results, CG_params = path_formulation_column_generation!(
    model,
    z,
    SR3_constraints,
    data,
    graph,
    artificial_paths,
    some_paths,
    path_costs,
    path_service,
    printlist,
    ;
    method = method,
    charge_cost_heterogenous = charge_cost_heterogenous,
    time_windows = time_windows,
    elementary = elementary,
    neighborhoods = neighborhoods,
    ngroute = ngroute,
    verbose = verbose,
    time_limit = 120.0,
    max_iters = Inf,
)
CGIP_results = path_formulation_solve_integer_model!(
    model,
    z,
    artificial_paths,
)
push!(CGLP_all_results, CGLP_results)
push!(CGIP_all_results, CGIP_results)
push!(CG_all_params, CG_params)
if ngroute
    push!(CG_all_neighborhoods, copy(neighborhoods))
end
CG_params["CGLP_objective"] = CGLP_results["objective"]
CG_params["CGIP_objective"] = CGIP_results["objective"]
CG_params["LP_IP_gap"] = 1.0 - CGLP_results["objective"] / CGIP_results["objective"]
CGLP_results["paths"] = collect_path_solution_support(
    CGLP_results, some_paths, data, graph
)
CGIP_results["paths"] = collect_path_solution_support(
    CGIP_results, some_paths, data, graph
)

println(CG_params["converged"])
# Expand ng-route neighborhoods
cycles_lookup = detect_cycles_in_path_solution([p for (val, p) in CGLP_results["paths"]], graph)
if length(cycles_lookup) > 0 
    delete_paths_with_found_cycles_from_model!(model, z, some_paths, path_costs, path_service, cycles_lookup, graph)
    modify_neighborhoods_with_found_cycles!(neighborhoods, cycles_lookup)
    continue_flag = true
    # iteration_params["method"] = "use_adaptive_ngroute"
    # iteration_params["cycles_lookup_length"] = length(cycles_lookup)
end
add_message!(printlist, "Expanded ng-route neighborhoods by $(sum(length.(values(cycles_lookup))))\n", verbose)
add_message!(printlist, "\n", verbose)


generated_SR3_list = enumerate_violated_path_SR3_inequalities(
    CGLP_results["paths"],
    graph,
)
generated_SR3_list = select_representative_violated_path_SR3_inequalities(
    generated_SR3_list, data,
)
implemented_SR3_list = generated_SR3_list[1:max_SR3_cuts]

append!(SR3_list, implemented_SR3_list)
# Add violated inequalities to master problem
add_SR3_constraints_to_path_model!(
    model, z, some_paths, 
    SR3_constraints, implemented_SR3_list, 
)
# iteration_params["method"] = "use_SR3_cuts"
# iteration_params["implemented_SR3_cuts_count"] = length(implemented_SR3_list)
# continue_flag = true
add_message!(printlist, "Imposed SR3 cuts:\t\t$(length(implemented_SR3_list))\n", verbose)
for (val, S) in implemented_SR3_list
    add_message!(printlist, "$S: $val\n", verbose)
end
add_message!(printlist, "\n", verbose)












Env = Gurobi.Env()
method = "ours"
charge_cost_heterogenous = false
time_windows = use_time_windows
elementary = false
ngroute = true
verbose = true
ngroute_neighborhood_size = Int(ceil(sqrt(graph.n_customers)))
ngroute_neighborhood_depots_size = "small"
ngroute_neighborhood_charging_size = "small"
max_SR3_cuts = 5

neighborhoods = compute_ngroute_neighborhoods(
    graph,
    ngroute_neighborhood_size; 
    depots_size = ngroute_neighborhood_depots_size,
    charging_size = ngroute_neighborhood_charging_size,
)

artificial_paths = generate_artificial_paths(data, graph)
some_paths = deepcopy(artificial_paths)
path_costs = compute_path_costs(
    data, graph, 
    some_paths,
)
path_service = compute_path_service(
    graph,
    some_paths,
)

model, z = path_formulation_build_model(
    data, graph, some_paths, path_costs, path_service,
    ; 
    Env = Env,
)
SR3_constraints = Dict{NTuple{3, Int}, ConstraintRef}()
SR3_list = Tuple{Float64, NTuple{3, Int}}[]
printlist = String[]


###
CG_params = Dict{String, Any}()
CG_params["number_of_paths"] = [sum(length(v) for v in values(some_paths))]
CG_params["objective"] = Float64[]
CG_params["κ"] = Dict{Int, Float64}[]
CG_params["μ"] = Dict{Int, Float64}[]
CG_params["ν"] = Vector{Float64}[]
CG_params["λ"] = Dict{keytype(SR3_constraints), Float64}[]
CG_params["lp_relaxation_solution_time_taken"] = Float64[]
CG_params["sp_base_time_taken"] = Float64[]
CG_params["sp_full_time_taken"] = Float64[]
CG_params["sp_total_time_taken"] = Float64[]
CG_params["lp_relaxation_constraint_time_taken"] = Float64[]
CG_params["number_of_new_paths"] = Int[]
CG_params["converged"] = false
CG_params["cycled"] = false





## Start loop

optimize!(model)
CGLP_results = Dict(
    "errored" => false,
    "objective" => objective_value(model),
    "z" => Dict(
        (key, p) => value.(z[(key, p)])
        for (key, p) in keys(z)
    ),
    "κ" => Dict(zip(graph.N_depots, dual.(model[:κ]).data)),
    "μ" => Dict(zip(graph.N_depots, dual.(model[:μ]).data)),
    "ν" => dual.(model[:ν]).data,
    "λ" => Dict{keytype(SR3_constraints), Float64}(
        S => dual(SR3_constraints[S])
        for S in keys(SR3_constraints)
    ),
)
CGLP_results["artificial"] = any(
    value.(z[(key, p)]) > 1e-3
    for key in keys(artificial_paths)
        for p in 1:length(artificial_paths[key])
)
push!(CG_params["objective"], CGLP_results["objective"])
push!(CG_params["κ"], CGLP_results["κ"])
push!(CG_params["μ"], CGLP_results["μ"])
push!(CG_params["ν"], CGLP_results["ν"])
push!(CG_params["λ"], CGLP_results["λ"])



κ = CGLP_results["κ"]
μ = CGLP_results["μ"]
ν = CGLP_results["ν"]
λ = CGLP_results["λ"]


include("subpath_stitching.jl")
(negative_full_labels1, negative_full_labels_count1, _, _) = subproblem_iteration_ours(
    data,
    graph,
    κ,
    μ,
    ν,
    λ,
    ;
    charge_cost_heterogenous = false,
    neighborhoods = neighborhoods,
    ngroute = true,
    elementary = false,
    # time_limit = 60.0,
);
generated_paths1 = get_paths_from_negative_path_labels(
    data, graph, negative_full_labels1,
);
num_new_paths1 = sum(length(v) for v in values(generated_paths1))


include("subpath_stitching_v2.jl")
(negative_full_labels2, negative_full_labels_count2, _, _) = subproblem_iteration_ours(
    data,
    graph,
    κ,
    μ,
    ν,
    λ,
    ;
    charge_cost_heterogenous = false,
    neighborhoods = neighborhoods,
    ngroute = true,
    elementary = false,
    # time_limit = 60.0,
);
generated_paths2 = get_paths_from_negative_path_labels(
    data, graph, negative_full_labels2,
);
num_new_paths2 = sum(length(v) for v in values(generated_paths2))

@assert num_new_paths1 == num_new_paths2

setdiff(
    keys(generated_paths1),
    keys(generated_paths2),
)

setdiff(
    keys(generated_paths2),
    keys(generated_paths1),
)

for (n1, n2) in keys(generated_paths1)
    if length(generated_paths1[(n1, n2)]) != length(generated_paths2[(n1, n2)])
        println("Difference in number of path labels for $n1 -> $n2")
        println("  Subpath stitching v1: $(length(generated_paths1[(n1, n2)]))")
        println("  Subpath stitching v2: $(length(generated_paths2[(n1, n2)]))")
    end
    for v1 in keys(generated_paths1[(n1, n2)])
        if all(!isequal(v1, v2) for v2 in keys(generated_paths2[(n1, n2)]))
            println("  - Paths in v1 but not v2: $v1")
        end
    end
    for v2 in keys(generated_paths2[(n1, n2)])
        if all(!isequal(v1, v2) for v1 in keys(generated_paths1[(n1, n2)]))
            println("  - Paths in v2 but not v1: $v2")
        end
    end
end



push!(
    CG_params["number_of_new_paths"],
    num_new_paths1
)
mp_constraint_time = add_paths_to_path_model!(
    model,
    z,
    some_paths, 
    path_costs,
    path_service,
    generated_paths1,
    SR3_constraints,
    data, graph,
)
push!(
    CG_params["number_of_paths"], 
    sum(length(v) for v in values(some_paths))
)
push!(
    CG_params["lp_relaxation_constraint_time_taken"],
    mp_constraint_time,
)








include("subpath_stitching.jl")
κ = CGLP_results["κ"]
μ = CGLP_results["μ"]
ν = CGLP_results["ν"]
λ = CGLP_results["λ"]
base_labels_result1 = @timed generate_base_labels_ngroute_lambda(
    data, graph, neighborhoods, 
    κ, μ, ν, λ,
    ;
    time_limit = 60.0,
);
subpath_labels1 = base_labels_result1.value;

full_labels_result = @timed find_nondominated_paths_notimewindows_ngroute_lambda(
    data, graph, neighborhoods, 
    subpath_labels1, κ, μ, λ,
    ;
    time_limit = 60.0,
);
path_labels1 = full_labels_result.value;

include("subpath_stitching_v2.jl")
subpath_labels2 = generate_subpath_labels(
    data,
    graph,
    κ, μ, ν,
    NoLoad(),
    NgRoute(),
    SR3Cuts(),
    ;
    neighborhoods = neighborhoods,
    λ = λ,
);


for (n1, n2) in keys(subpath_labels1)
    # if length(subpath_labels1[(n1, n2)]) != length(subpath_labels2[(n1, n2)])
    #     println("Difference in number of subpath labels for $n1 -> $n2")
    #     println("  Subpath stitching v1: $(length(subpath_labels1[(n1, n2)]))")
    #     println("  Subpath stitching v2: $(length(subpath_labels2[(n1, n2)]))")
    # end
    for k1 in keys(subpath_labels1[(n1, n2)])
        if all(k1[1:3] != k2[1:3] for k2 in keys(subpath_labels2[(n1, n2)]))
            println("  - $n1 -> $n2: Subpaths in v1 but not v2: $k1")
        end
    end
    for k2 in keys(subpath_labels2[(n1, n2)])
        if all(k1[1:3] != k2[1:3] for k1 in keys(subpath_labels1[(n1, n2)]))
            println("  - $n1 -> $n2: Subpaths in v2 but not v1: $k2")
        end
    end
end

path_labels2 = generate_path_labels_notimewindows(
    data,
    graph,
    NoLoad(),
    HomCharge(),
    NgRoute(),
    NoCuts(),
    subpath_labels2,
    ;
    λ = λ,
    # time_limit = time_limit - (time() - start_time),
);

for (n1, n2) in keys(path_labels1)
    if length(path_labels1[(n1, n2)]) != length(path_labels2[(n1, n2)])
        println("Difference in number of path labels for $n1 -> $n2")
        println("  Subpath stitching v1: $(length(path_labels1[(n1, n2)]))")
        println("  Subpath stitching v2: $(length(path_labels2[(n1, n2)]))")
    end
    for k1 in keys(path_labels1[(n1, n2)])
        if all(k1[1:3] != k2[1:3] for k2 in keys(path_labels2[(n1, n2)]))
            println("  - Paths in v1 but not v2: $k1")
        end
    end
    for k2 in keys(path_labels2[(n1, n2)])
        if all(k1[1:3] != k2[1:3] for k1 in keys(path_labels1[(n1, n2)]))
            println("  - Paths in v2 but not v1: $k2")
        end
    end
end

include("subpath_stitching.jl")
negative_full_labels1 = get_negative_path_labels_from_path_labels(path_labels1)

generated_paths1 = get_paths_from_negative_path_labels(
    data, graph, negative_full_labels1,
)
num_new_paths1 = sum(length(v) for v in values(generated_paths1))

include("subpath_stitching_v2.jl")
negative_full_labels2 = get_negative_path_labels_from_path_labels(path_labels2)
generated_paths2 = get_paths_from_negative_path_labels(
    data, graph, negative_full_labels2,
)
num_new_paths2 = sum(length(v) for v in values(generated_paths2))

@assert num_new_paths1 == num_new_paths2

setdiff(
    keys(generated_paths1),
    keys(generated_paths2),
)

setdiff(
    keys(generated_paths1),
    keys(generated_paths2),
)

for (n1, n2) in keys(generated_paths1)
    if length(generated_paths1[(n1, n2)]) != length(generated_paths2[(n1, n2)])
        println("Difference in number of path labels for $n1 -> $n2")
        println("  Subpath stitching v1: $(length(generated_paths1[(n1, n2)]))")
        println("  Subpath stitching v2: $(length(generated_paths2[(n1, n2)]))")
    end
    for v1 in keys(generated_paths1[(n1, n2)])
        if all(!isequal(v1, v2) for v2 in keys(generated_paths2[(n1, n2)]))
            println("  - Paths in v1 but not v2: $v1")
        end
    end
    for v2 in keys(generated_paths2[(n1, n2)])
        if all(!isequal(v1, v2) for v1 in keys(generated_paths1[(n1, n2)]))
            println("  - Paths in v2 but not v1: $v2")
        end
    end
end


push!(
    CG_params["number_of_new_paths"],
    num_new_paths2
)
mp_constraint_time = add_paths_to_path_model!(
    model,
    z,
    some_paths, 
    path_costs,
    path_service,
    generated_paths2,
    SR3_constraints,
    data, graph,
)
push!(
    CG_params["number_of_paths"], 
    sum(length(v) for v in values(some_paths))
)
push!(
    CG_params["lp_relaxation_constraint_time_taken"],
    mp_constraint_time,
)

# ---

κ = CGLP_results["κ"]
μ = CGLP_results["μ"] 

s1 = subpath_labels1[35, 40] |> values |> first
s2 = subpath_labels2[35, 40] |> values |> first

s1.cost
s2.cost

(s1.cost - s2.cost) / 2

λ

subpath_labels1[35, 40] |> values |> collect |> last
subpath_labels2[35, 40] |> values |> collect |> x -> x[end-1]

subpath_labels2[35, 40]

modified_costs = compute_arc_modified_costs(graph, data, ν)
λvals, λcust = prepare_lambda(λ, graph.n_nodes)
use_load = NoLoad()
charge_costs = HomCharge()
customer_service = NgRoute()
cuts = SR3Cuts()

nodeseq = [35, 13, 5, 9, 40]

starting_node = nodeseq[1]
current_subpath = create_new_subpath_label(
    use_load,
    customer_service,
    cuts,
    starting_node,
    data,
    ;
    n_cuts = length(λ), 
    neighborhoods = neighborhoods,
    λmemory = (
        cuts isa LmSR3Cuts 
        ? λmemory 
        : falses(length(λ), graph.n_nodes)
    ), 
)

i = 1
current_node, next_node = nodeseq[i], nodeseq[i+1]
(feasible, new_subpath) = compute_new_subpath(
    current_subpath,
    graph,
    modified_costs,
    current_node,
    next_node, 
    NoLoad(),
    NgRoute(),
    SR3Cuts(),
    ;
    neighborhoods = neighborhoods,
    n_cuts = length(λ),
    λvals = λvals,
    λcust = λcust,
    λmemory = falses(length(λ), graph.n_nodes)
)
println(new_subpath.cut_labels)
@show current_subpath.cost + modified_costs[current_node, next_node]
@show new_subpath.cost
println(new_subpath.cost)
println(current_subpath.cost)

current_subpath = new_subpath
i += 1
