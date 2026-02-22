using Pkg
Pkg.activate("$(@__DIR__)/../")

include("$(@__DIR__)/../src/utils.jl")
include("$(@__DIR__)/../src/read_evrptw_data.jl")

using Glob
using Plots

for instance_type in [
    "c1", "r1", "rc1", 
    "c2", "r2", "rc2",
]
    for n_customers in [25, 50, 100]
        for fp in Glob.glob(
            "$(instance_type)*_$n_customers.txt",
            "data/evrptw/Instances", 
        )
            println(fp)
            filename = basename(fp)[1:end-4]
            data = read_evrptw_instance(
                filename,
                5,
                1,
                ;
                data_dir = "data/evrptw",
            )
            graph = generate_graph_from_data(data)
            p = plot_instance(data; graph = graph)
            if instance_type in ["r1", "r2"]
                plot!(p, xlims = (-5, 75), ylims = (-5, 85))
            elseif instance_type in ["c1", "c2", "rc1", "rc2"]
                plot!(p, xlims = (-5, 105), ylims = (-5, 105))
            end
            Plots.plot!(
                p,
                xlabel = "x",
                ylabel = "y",
                legend = :outerbottom,
                margin = 2Plots.mm,
                size = (500, 500),
            )
            Plots.savefig(p, "$(@__DIR__)/../data/evrptw_plots/$(filename).png")
            
            sorted_TWs = hcat(sort(collect.(zip(data.α[data.N_customers], data.β[data.N_customers])))...) ./ data.T
            p = Plots.plot(
                data.N_customers,
                [
                    sorted_TWs[1, :],
                    sorted_TWs[2, :],
                ],
                # seriestype = :bar,
                # bins = 20,
                label = false,
                title = "Time windows for customers in $(filename)",
            )

            Plots.savefig(p, "$(@__DIR__)/../data/evrptw_plots/$(filename)_TW.png")
        end
    end
end
