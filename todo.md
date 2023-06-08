# TODO:

## Code
- [x] acceleration strategy: ensure that it gets the same objective as non-accelerated
    - note: it might not, because we generate paths but solve the subpath relaxation. 
- [x] investigate instances of inequalities not satisfied: $$\text{CG}_{\text{label}}(\mathcal{P}) \geq \text{CG}_{2,\text{label}}(\mathcal{S}) \geq \text{CG}_{1,\text{label}}(\mathcal{S})$$
- [ ] unify code in `check_customers = true` and `check_customers = false` in `generate_base_labels()`
- [x] Investigate why non-elementary paths themselves don't converge
- [ ] Move code to `path_formulation.jl`, and implement the column generation on path vars
- [ ] Move back to slower (non-Floyd-Warshall) code in `generate_base_labels()`, in order to generate non-elementary subpaths w/ 2-loops
    - note: maybe we don't care about 2-loops!
- [ ] Investigate if there exists faster stitching code (buckets idea)

## Research
- [ ] Think if there exists problem settings of VRPs w/o TW and w/o capacity
- [ ] Investigate if partial node labels / iteratively expanding partial set of node labels can be leveraged 
