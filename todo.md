# TODO:

## Code
- [x] modify code in `desaulniers_benchmark.jl`: add history of 1 for all code
- [ ] check that lower bounds after cuts never go down
    - [x] in `subpath_stitching.jl`, for `generate_base_labels_ngroute[_alt]_sigma`: requires "first_node" (only if the subpath starts at a CS) and history of 1 (due to interaction with Christofides)
- [ ] think about how to implement "smart Christofides": for each label, 2 paths:
    - the actual best path, and
    - a second path with those parameters which does not visit the second-to-last node of the best path
- [ ] investigate finite termination issue after first round of cuts
    - [ ] `add_paths_to_path_model!`: need to add the paths to the relevant violated WSR3 inequalities (obviously!)
- [ ] commit 
- [ ] implement adaptive neighborhoods idea
