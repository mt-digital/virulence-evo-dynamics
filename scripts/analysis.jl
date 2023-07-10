using CSV
using DataFrames
using RCall


include("../src/model.jl")


function run_series(plot = true; metapop_size = 1000,
                    maxsteps = 10000, rename = true, whensteps = 10,
                    virulence_init = 0.3, 
                    initial_infected_frac = 0.05, 
                    mutation_rate = 0.8,
                    mutation_variance = 0.2,
                    virulence_mortality_coeff = 0.01,
                    virulence_transmission_coeff = 0.01,
                    virulence_transmission_denom_summand = 0.3, 
                    min_group_frac = 0.4,
                    min_start = true,
                    maj_start = true,
                    global_add_rate = 2, global_death_rate = 2, 
                    min_homophily = 0.0, maj_homophily = 0.0,
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
    m = virulence_evo_model(; 
                            metapop_size,
                            initial_infected_frac, 
                            virulence_init,
                            mutation_rate,
                            mutation_variance,
                            virulence_mortality_coeff,
                            virulence_transmission_coeff,
                            virulence_transmission_denom_summand, 
                            min_group_frac,
                            min_start,
                            maj_start,
                            global_add_rate, global_death_rate, 
                            min_homophily, maj_homophily,
                            model_params...); 

    # Run model, collecting agent and model dataframes.
    when(model, step) = step % whensteps == 0
    agent_df, model_df = run!(m, agent_step!, model_step!, stopfn; adata, mdata, when);

    println("\nTotal infected: $(model_df.total_infected[end])\n")

    println(names(agent_df))

    if rename
        # rename!(agent_df, :susceptible_status => :susceptible,
        #                   :infected_status => :infected,
        #                   :filtermeanvirulence_pathogen => :mean_virulence)
                          # "infected_status_#84_f=is_minority" => :infected_majority);
        rename!(agent_df, "susceptible_status" => "susceptible",
                          "infected_status" => "infected",
                          "filtermeanvirulence_pathogen" => "mean_virulence",

                          "susceptible_status_is_minority" => 
                              "susceptible_minority",
                          "infected_status_is_minority" => 
                              "infected_minority",
                          "filtermeanvirulence_pathogen_is_minority" => 
                              "mean_virulence_minority",

                          "susceptible_status_#84_f=is_minority" => 
                              "susceptible_majority",
                          "infected_status_#84_f=is_minority" => 
                              "infected_majority",
                          "filtermeanvirulence_pathogen_#84_f=is_minority" =>
                              "mean_virulence_majority"
                           );
    end

    sum_series_copy = copy(agent_df.susceptible + agent_df.infected) 
    agent_df.susceptible ./= sum_series_copy
    agent_df.infected ./= sum_series_copy

    sum_series_copy = copy(agent_df.susceptible_minority 
                           + agent_df.infected_minority) 
    agent_df.susceptible_minority ./= sum_series_copy 
    agent_df.infected_minority ./= sum_series_copy 

    sum_series_copy = copy(agent_df.susceptible_majority 
                           + agent_df.infected_majority) 
    agent_df.susceptible_majority ./= sum_series_copy 
    agent_df.infected_majority ./= sum_series_copy

    if plot
        println("here")
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
