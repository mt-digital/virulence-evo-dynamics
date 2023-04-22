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
                                    initial_infected_frac = 0.005, 
                                    mutation_rate = [0.0, 0.8],
                                    mutation_variance = 0.2,
                                    virulence_mortality_coeff = 0.01,
                                    virulence_transmission_coeff = 0.01,
                                    virulence_transmission_denom_summand = 0.3, 
                                    # global_add_rate = [0,2], global_death_rate = [0,2], 
                                    global_add_rate = 2, global_death_rate = 2, 
                                    whensteps = 10
    )

    rep_idx = collect(1:nreplicates)

    params_list = dict_list(
        @dict virulence_init mutation_rate global_add_rate rep_idx
    )

    models = [
        virulence_evo_model(; 
            metapop_size, initial_infected_frac, mutation_variance, 
            virulence_mortality_coeff, virulence_transmission_coeff,
            virulence_transmission_denom_summand, params...)
        for params in params_list
    ]

    # Define aggregation of infection status, first susceptible...
    susceptible(status_vec) = isempty(status_vec) ? 
        0.0 : 
        count(i == Susceptible for i in status_vec) / 
            count(i != Dead for i in status_vec)

    # ...then infected.
    infected(status_vec) = isempty(status_vec) ? 
        0.0 : 
        count(i == Infected for i in status_vec) / 
            count(i != Dead for i in status_vec)

    # Define aggregation for pathogen virulence.
    virulence(agent) = agent.status == Dead ? NaN : agent.pathogen.virulence
    filtermean(virulence_vec) = mean(filter(!isnan, collect(virulence_vec)))


    # Put all agent data aggregation together.
    # append!(adata, [(virulence, filtermean)])
    adata = [(:status, susceptible), (:status, infected), (virulence, filtermean)]

    # Track total infected over time.
    mdata = [:mutation_rate, :virulence_init, :total_infected]

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
