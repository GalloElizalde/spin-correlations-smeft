# Cargar paquetes
using EFTfitter
using BAT
using IntervalSets
using Distributions
using Plots
using Printf

# ===== Input Directory =====
input_dir = "cms_inputs"


# look for files with name: "cms_800GeV_*.jl"
input_files = filter(f -> endswith(f, ".jl") && occursin("cms_800GeV", f), readdir(input_dir, join=true))


# Results Directory
mkpath("EFTfitter_results")


# Loop over each file in "input dir"
for input_path in input_files
    println("\n>>> Analizando: $(input_path)")

    #  Define parameters, measurements, correlations...
    include(input_path)


    # Make Model
    model = EFTfitterModel(parameters, measurements, correlations, nuisances = nuisance_correlations)
    posterior = PosteriorMeasure(model)


    # Sampling
    algorithm = MCMCSampling(mcalg = MetropolisHastings(), nsteps = 10^5, nchains = 4)
    samples = bat_sample(posterior, algorithm).result


    # Base name for saved results
    name_base = split(basename(input_path), ".")[1]


    # Save Posterior Prob. plots
    savefig(plot(samples), "EFTfitter_results/$(name_base)_posterior.pdf")
    savefig(plot(samples, 0.9), "EFTfitter_results/$(name_base)_intervals.pdf")


    # Confidence Interval 90% 
    open("EFTfitter_results/$(name_base)_intervals.txt", "w") do io
        for pname in keys(parameters)
            intervals = get_smallest_interval_edges(samples, pname, 0.9, bins=200, atol=0.01)
            println("Intervalo 90% para $(pname):")
            println("  lower: $(intervals.lower)")
            println("  upper: $(intervals.upper)")

            println(io, "90% interval for $(pname):")
            println(io, "  lower: $(intervals.lower)")
            println(io, "  upper: $(intervals.upper)")
            println(io)  # línea en blanco entre parámetros
        end
    end

    # 1D Interval 
    #pname1 = first(keys(parameters))
    #intervals = get_smallest_interval_edges(samples, pname1, 0.9, bins=200, atol=0.01)
    #println("Intervalo 90% para $(pname1):")
    #println("  lower: $(intervals.lower)")
    #println("  upper: $(intervals.upper)")

    # Save Intervals in .txt
    #open("EFTfitter_results/$(name_base)_intervals_costheta.txt", "w") do io
    #    println(io, "90% interval for $(pname1):")
    #    println(io, "  lower: $(intervals.lower)")
    #    println(io, "  upper: $(intervals.upper)")
    #end

    # Ranking
    ranking = EFTfitter.rank_measurements(model)
    savefig(plot(ranking, title = "Ranking of measurements"),
            "EFTfitter_results/$(name_base)_ranking.pdf")
end

println("\n✓ Finish. Check 'EFTfitter_results'.")
