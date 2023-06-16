using DataFrames

using Distributed
using StatsBase


# Set up multiprocessing.
try
    num_cores = parse(Int, ENV["SLURM_CPUS_PER_TASK"])
    addprocs(num_cores)
catch
    desired_nprocs = length(Sys.cpu_info())

    if length(procs()) != desired_nprocs
        addprocs(desired_nprocs - 1)
    end
end


@everywhere using DrWatson
@everywhere quickactivate("..")
@everywhere include("model.jl")


function virulence_evo_experiment(nreplicates = 10, record_series = false; 
                                    maxsteps = 20_000, metapop_size = 2_000,
                                    virulence_init = collect(0.01:0.1:0.91),
                                    initial_infected_frac = 0.05, 
                                    mutation_rate = [0.0, 0.8],
                                    mutation_variance = 0.2,
                                    virulence_mortality_coeff = 0.01,
                                    virulence_transmission_coeff = 0.01,
                                    virulence_transmission_denom_summand = 0.3, 
                                    min_group_frac = 0.2,
                                    min_start = true,
                                    maj_start = false,
                                    global_add_rate = 2, global_death_rate = 2, 
                                    min_homophily = 0.0, maj_homophily = 0.0,
                                    whensteps = 10
    )

    rep_idx = collect(1:nreplicates)

    params_list = dict_list(
        @dict virulence_init mutation_rate global_add_rate rep_idx min_homophily maj_homophily
    )

    models = [
        virulence_evo_model(; 
            metapop_size, initial_infected_frac, mutation_variance, 
            virulence_mortality_coeff, virulence_transmission_coeff,
            virulence_transmission_denom_summand, min_group_frac, min_start,
            maj_start, params...)
        for params in params_list
    ]

    # Define aggregation of infection status, first susceptible...
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
    # virulence(pathogen_vec) = mean(filter(
    # virulence(agent) = agent.pathogen.virulence
    # filtermean(virulence_vec) = mean(filter(!isnan, collect(virulence_vec)))
    # filtermean(x) = mean(filter(!isnan, x))
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

    # function filtermeanvirulence_minority(person_vec)
    #     virulence(person) = person.pathogen.virulence
    #     virulence_vec = [virulence(p) for p in person_vec if p.group == Minority]
    #     filter!(!isnan, virulence_vec)
    #     return mean(virulence_vec)
    # end

    # function filtermeanvirulence_majority(person_vec)
    #     virulence(person) = person.pathogen.virulence
    #     virulence_vec = [virulence(p) for p in person_vec if p.group == Majority]
    #     filter!(!isnan, virulence_vec)
    #     return mean(virulence_vec)
    # end

    # Define group filtering.
    is_minority(x) = x.group == Minority


    # Put all agent data aggregation together.
    # append!(adata, [(virulence, filtermean)])
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
    mdata = [:mutation_rate, :virulence_init, :total_infected, :total_minority_infected, :total_majority_infected]

    # Stop when the pathogen has gone extinct or maxsteps reached.
    stopfn(model, step) = (count(
        agent.status == Infected for agent in collect(allagents(model))
       ) == 0) || (step == maxsteps)

    # Don't want to record every step in the series for full experiments.
    if record_series
        when(model, step) = step % whensteps == 0
    else
        when = stopfn
    end

    adf, mdf = ensemblerun!(models, agent_step!, model_step!, stopfn;
                            adata, mdata, when, parallel = true, 
                            showprogress = true)
    
    return innerjoin(adf, mdf, on = [:step, :ensemble])
end
