using CSV
using DataFrames
using RCall


include("../src/model.jl")


function run_series(plot = true; maxsteps = 1000, rename = true, model_params...)
    
    susceptible(x) = isempty(x) ? 0.0 : count(i == Susceptible for i in x) 
    infected(x) = isempty(x) ? 0.0 : count(i == Infected for i in x) 
    # recovered(x) = isempty(x) ? 0.0 : count(i == Recovered for i in x) 

    # adata = [(:status, f) for f in (susceptible, infected, recovered)]
    adata = [(:status, f) for f in (susceptible, infected)]

    # adata = [:status]
    mdata = [:total_infected]

    stopfn(model, step) = count(
        agent.status == Infected for agent in collect(allagents(model))
    ) == 0

    m = virulence_evo_model(; model_params...); 

    agent_df, model_df = run!(m, agent_step!, model_step!, stopfn; adata, mdata);

    println("\nTotal infected: $(model_df.total_infected[end])\n")

    if rename
        rename!(agent_df, :susceptible_status => :susceptible,
                          :infected_status => :infected);
    end

    if plot
        plot_series(agent_df)
    end

    return agent_df, model_df
end


function plot_series(agent_df)

    # Create temporary directory for CSV and plot output.
    if !isdir("tmp")
        mkdir("tmp")
    end

    # Write temporary CSV for plotting.
    CSV.write("tmp/series.csv", agent_df)

    # Call R plotting routine.
    R"""
    source("scripts/plot.R")

    plot_series("tmp/series.csv")
    """
end
