using CSV
using DataFrames
using RCall


include("../src/model.jl")


function run_series(plot = true; maxsteps = 3000, rename = true, whensteps = 100,
                                 model_params...)
    
    # Define aggregation of infection status.
    # susceptible(status_vec) = isempty(status_vec) ? 0.0 : count(i == Susceptible 
    #                                                             for i in status_vec) / length(status_vec)
    # infected(status_vec) = isempty(status_vec) ? 0.0 : count(i == Infected 
    #                                                          for i in status_vec) / length(status_vec)
    # adata = [(a -> a.status, f) for f in (susceptible, infected)]
    susceptible(status_vec) = isempty(status_vec) ? 
        0.0 : count(i == Susceptible for i in status_vec) 

    # ...then infected.
    # infected(status_vec) = isempty(status_vec) ? 
    #     0.0 : 
    #     count(i == Infected for i in status_vec)
    function infected(status_vec)
        # println(status_vec)
        return isempty(status_vec) ?  0.0 : count(i == Infected for i in status_vec)
    end

    # Define aggregation for pathogen virulence.
    # virulence(agent) = agent.pathogen.virulence
    # filtermean(virulence_vec) = mean(filter(!isnan, collect(virulence_vec)))

    is_minority(x) = x.group == Minority

    function filtermeanvirulence(pathogen_vec)
        # virulence(pathogen) = pathogen.virulence
        # virulence_vec = [pathogen.virulence for pathogen in pathogen_vec]
        # println(pathogen_vec)
        # println(typeof(pathogen_vec))
        # println(length(collect(pathogen_vec)))
        virulence_vec = [p.virulence for p in pathogen_vec if !isnan(p.virulence)]
        if isempty(virulence_vec) 
            return 0.0
        else
            return mean(virulence_vec)
        end
    end

    # Put all agent data aggregation together.
    # append!(adata, [(virulence, filtermean)])
    # adata = [(:status, susceptible), (:status, infected), (virulence, filtermean)]
    adata = [
             (:status, susceptible), (:status, infected), 
             (:pathogen, filtermeanvirulence),

             (:status, susceptible, is_minority), 
             (:status, infected, is_minority), 
             (:pathogen, filtermeanvirulence, is_minority),

             (:status, susceptible, !is_minority), 
             (:status, infected, !is_minority), 
             (:pathogen, filtermeanvirulence, !is_minority) #, !is_minority)
             # (virulence, filtermean, !is_minority)
            ]
    # Track total infected over time.
    mdata = [:total_infected]

    # Stop when the pathogen has gone extinct.
    stopfn(model, step) = (count(
        agent.status == Infected for agent in collect(allagents(model))
       ) == 0) || (step == maxsteps)

    # Initialize model.
    m = virulence_evo_model(; model_params...); 

    # Run model, collecting agent and model dataframes.
    when(model, step) = step % whensteps == 0
    agent_df, model_df = run!(m, agent_step!, model_step!, stopfn; adata, mdata, when);

    println("\nTotal infected: $(model_df.total_infected[end])\n")

    println(names(agent_df))

    if rename
        rename!(agent_df, :susceptible_status => :susceptible,
                          :infected_status => :infected,
                          :filtermeanvirulence_pathogen => :mean_virulence);
    end

    if plot
        plot_series(agent_df)
    end

    return agent_df, model_df
end


function run_trials(nreplicates = 10; outputfilename, experiment_kwargs...)

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
