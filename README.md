# Subpath-based Column Generation for the Electric Routing-Scheduling Problem (ERSP)

This repository implements the algorithm and experiment suite for the paper: "Subpath-based Column Generation for the Electric Routing-Scheduling Problem" by Alexandre Jacquillat and Sean Lo. The paper can be found [here](https://arxiv.org/abs/2407.02640).

Please cite this paper as:
```
@misc{jacquillat2024subpath,
      title={Subpath-Based Column Generation for the Electric Routing-Scheduling Problem}, 
      author={Alexandre Jacquillat and Sean Lo},
      year={2024},
      eprint={2407.02640},
      archivePrefix={arXiv},
      primaryClass={math.OC},
      url={https://arxiv.org/abs/2407.02640}, 
}
```
<!-- TODO: add links to paper -->

## Usage

1) Create an instance `data::EVRPData` (found in `src/utils.jl`), such as via `generate_instance()`.
2) Create a graph `graph::EVRPGraph` from the problem instance using `generate_graph_from_data()`.
3) Implement column generation on the set-partitioning formulation, augmented with *ng*-routes and *lm*-SRI cuts, via the function `path_formulation_column_generation_with_adaptve_ngroute_SR3_cuts()` (found in `src/path_formulation.jl`)

Important arguments:
- `method::String = "ours"`: can be either `"ours"` or `"benchmark"`, where the benchmark implements the method by [Desaulniers et al (2016)](https://pubsonline.informs.org/doi/10.1287/opre.2016.1535)
- `charge_cost_heterogenous::Bool = false`: whether heterogenous charging costs at charging stations are modelled
- `verbose::Bool = true`: whether to print output as the algorithm progresses
- `time_limit::Float64 = Inf`: time limit in seconds for the overall algorithm.
- Path elementarity:
    - `elementary::Bool = false`: whether to impose binary labels for each customer / task, according to Beasley and Christofides (1989)
    - `ngroute::Bool = true`: whether to solve the *ng*-route relaxation [(Baldacci et al 2011)](https://pubsonline.informs.org/doi/10.1287/opre.1110.0975)
    - `neighborhoods::Union{Nothing, BitMatrix} = nothing` if solving the *ng*-route relaxation, the choice of *ng*-neighborhoods at each node: `neighborhoods[i][j] = 1` iff node `i` is in the *ng*-neighborhood of node `j` 
    - `ngroute_neighborhood_size::Int = Int(ceil(sqrt(graph.n_customers)))`: the number of customer nodes to include in the *ng*-set of each node
    - `ngroute_neighborhood_depots_size::String` and `ngroute_neighborhood_charging_size::String` determines the size of *ng*-sets at depots and charging stations respectively, can be one of `"small"`, `"medium"`, `"large"`. Recommended value `"small"`
    - `use_adaptive_ngroute::Bool = true`: whether to adaptively grow *ng*-relaxations according to [Martinelli et al (2014)](https://www.sciencedirect.com/science/article/abs/pii/S0377221714004123)
- Cutting planes
    - `use_SR3_cuts::Bool = true`: whether to use subset-row inequalities (includes *lm*-SRIs), first introduced in [Jepsen et al (2008)](https://pubsonline.informs.org/doi/10.1287/opre.1070.0449)
    - `use_lmSR3_cuts::Bool = true`: whether to use *lm*-SRIs instead of SRIs, according to [Pecin et al (2017)](https://link.springer.com/article/10.1007/s12532-016-0108-8)
    - `max_SR3_cuts::Int = 5`: how many cuts to add after each cut separation procedure; recommended to be a small number since cuts increase the pricing problem complexity
    - `randomize_cuts::Bool = false`: if `true`, selects which cuts to add in a random manner (weighted by violation of incumbent solution); if `false`, adds cuts greedily according to violation.